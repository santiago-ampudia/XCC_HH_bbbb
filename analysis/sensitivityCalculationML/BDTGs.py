import os
import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
from tqdm import tqdm
import argparse
import sys
import subprocess

np.random.seed(42)

def load_data(signal_file, background_files, test_size=0.25):
    signal_data = pd.read_csv(signal_file)
    signal_data['label'] = 1  # Add label column with value 1 (signal)
    
    print(f"\nSignal file: {os.path.basename(signal_file)}")
    print(f"Total signal events: {len(signal_data)}")
    
    feature_names = [col for col in signal_data.columns if col != 'label']
    
    #same signal events for all background types
    signal_train, signal_test = train_test_split(
        signal_data, test_size=test_size, random_state=42
    )
    
    print(f"  Signal events for training: {len(signal_train)}")
    print(f"  Signal events for testing: {len(signal_test)}")
    
    signal_train.to_csv('signal_train.csv', index=False)
    signal_test.to_csv('signal_test.csv', index=False)
    print(f"  Saved signal training events to signal_train.csv")
    print(f"  Saved signal testing events to signal_test.csv")
    
    background_dfs = {}
    
    print("\nBackground files:")
    for bg_file in background_files:
        bg_name = os.path.splitext(os.path.basename(bg_file))[0]
        try:
            bg = pd.read_csv(bg_file)
            bg['label'] = 0  
            
            background_dfs[bg_name] = bg
            
            print(f"  {bg_name}: {len(bg)} events")
        except FileNotFoundError:
            print(f"  Warning: {bg_name} file not found, skipping")
    
    datasets = {}
    
    for bg_name, bg_data in background_dfs.items():
        bg_train, bg_test = train_test_split(
            bg_data, test_size=test_size, random_state=42
        )
        
        bg_train.to_csv(f'{bg_name}_train.csv', index=False)
        bg_test.to_csv(f'{bg_name}_test.csv', index=False)
        print(f"  Saved {bg_name} training events to {bg_name}_train.csv")
        print(f"  Saved {bg_name} testing events to {bg_name}_test.csv")
        
        X_train = pd.concat([signal_train[feature_names], bg_train[feature_names]])
        y_train = pd.concat([signal_train['label'], bg_train['label']])
        
        X_test = pd.concat([signal_test[feature_names], bg_test[feature_names]])
        y_test = pd.concat([signal_test['label'], bg_test['label']])
        
        # feature normalization
        scaler = StandardScaler()
        
        X_train_scaled = pd.DataFrame(
            scaler.fit_transform(X_train),
            columns=X_train.columns,
            index=X_train.index
        )
        X_test_scaled = pd.DataFrame(
            scaler.transform(X_test),
            columns=X_test.columns,
            index=X_test.index
        )
        
        datasets[bg_name] = {
            'X_train': X_train_scaled,
            'y_train': y_train.values,
            'X_test': X_test_scaled,
            'y_test': y_test.values,
            'scaler': scaler,
            'feature_names': feature_names
        }
        
        train_df = X_train.copy()
        train_df['label'] = y_train.values
        test_df = X_test.copy()
        test_df['label'] = y_test.values
        
        train_df.to_csv(f'signal_vs_{bg_name}_train.csv', index=False)
        test_df.to_csv(f'signal_vs_{bg_name}_test.csv', index=False)
        print(f"  Saved combined signal vs {bg_name} training events to signal_vs_{bg_name}_train.csv")
        print(f"  Saved combined signal vs {bg_name} testing events to signal_vs_{bg_name}_test.csv")
        
        print(f"\nDataset for {bg_name} vs Signal:")
        print(f"  Training: {len(X_train_scaled)} events ({len(signal_train)} signal, {len(bg_train)} background)")
        print(f"  Testing: {len(X_test_scaled)} events ({len(signal_test)} signal, {len(bg_test)} background)")
    
    return datasets, signal_test, {bg_name: train_test_split(bg, test_size=test_size, random_state=42)[1] 
                                  for bg_name, bg in background_dfs.items()}

def train_bdtg_models(datasets, n_estimators=200, learning_rate=0.1, max_depth=5):
    models = {}
    
    print("\nTraining XGBoost models for each background type vs signal:")
    for bg_name, data in tqdm(datasets.items(), desc="Background types"):
        print(f"\n  Training XGBoost for {bg_name} vs Signal")
        
        X_train = data['X_train']
        y_train = data['y_train']
        X_test = data['X_test']
        y_test = data['y_test']
        
        # for early stopping
        X_train_main, X_valid, y_train_main, y_valid = train_test_split(
            X_train, y_train, test_size=0.2, random_state=42
        )
        
        feature_selector = xgb.XGBClassifier(
            n_estimators=50,  # CHANGE??? fewer trees for feature selection
            learning_rate=0.1,
            max_depth=3,
            random_state=42,
            n_jobs=-1,
            tree_method='hist',
            subsample=0.8,  
            colsample_bytree=0.8  
        )
        
        feature_selector.fit(X_train_main, y_train_main)
        
        importances = feature_selector.feature_importances_
        feature_names = X_train.columns.tolist()
        
        sorted_idx = np.argsort(importances)[::-1]
        
        cumulative_importance = 0
        important_features = []
        
        for idx in sorted_idx:
            important_features.append(feature_names[idx])
            cumulative_importance += importances[idx]
            if cumulative_importance >= 0.8:  
                break
        
        print(f"  Selected {len(important_features)} out of {len(feature_names)} features")
        
        X_train_main_selected = X_train_main[important_features]
        X_valid_selected = X_valid[important_features]
        X_test_selected = X_test[important_features]
        
        try:
            import cupy
            gpu_available = True
            tree_method = 'gpu_hist'
            predictor = 'gpu_predictor'
            print("  Using GPU acceleration")
        except ImportError:
            gpu_available = False
            tree_method = 'hist'
            predictor = 'auto'
            print("  Using CPU (GPU not available)")
        
        scale_pos_weight = np.sum(y_train_main == 0) / np.sum(y_train_main == 1)
        
        model = xgb.XGBClassifier(
            n_estimators=n_estimators,
            learning_rate=learning_rate,
            max_depth=max_depth,
            min_child_weight=3,
            gamma=0.1,
            subsample=0.8,
            colsample_bytree=0.8,
            colsample_bylevel=0.8,
            colsample_bynode=0.8,
            reg_alpha=0.01,
            reg_lambda=1,
            scale_pos_weight=scale_pos_weight,
            random_state=42,
            tree_method=tree_method,
            predictor=predictor,
            n_jobs=-1,
            grow_policy='lossguide',
            max_bin=256,
            early_stopping_rounds=50,
            eval_metric=['auc', 'logloss']
        )
        
        # for early stopping
        eval_set = [(X_train_main_selected, y_train_main), (X_valid_selected, y_valid)]
        
        print("  Training model with early stopping...")
        model.fit(
            X_train_main_selected, 
            y_train_main,
            eval_set=eval_set,
            verbose=100  # Print progress every 100 iterations
        )
        
        best_iteration = model.best_iteration
        if best_iteration is not None:
            print(f"  Best iteration: {best_iteration}")
        else:
            print("  Early stopping not triggered, used all iterations")
        
        y_pred = model.predict_proba(X_test_selected)[:, 1]
        auc = roc_auc_score(y_test, y_pred)
        
        print(f"  AUC for {bg_name}: {auc:.4f}")
        
        models[bg_name] = {
            'model': model,
            'auc': auc,
            'features': important_features
        }
    
    return models

def apply_models(models, datasets, signal_test, background_tests):
   
    print("\nApplying models to all event topologies and creating prediction files")
    
    first_bg = next(iter(datasets))
    feature_names = datasets[first_bg]['feature_names']
    
    predictions = {}
    
    signal_predictions = {}
    for bg_name, model_info in models.items():
        model = model_info['model']
        
        scaler = datasets[bg_name]['scaler']
        
        selected_features = model_info['features']
        
        X_signal = signal_test[feature_names]
        X_signal_scaled = pd.DataFrame(
            scaler.transform(X_signal),
            columns=X_signal.columns,
            index=X_signal.index
        )
        
        X_signal_selected = X_signal_scaled[selected_features]
        
        preds = model.predict_proba(X_signal_selected)[:, 1]
        
        signal_predictions[bg_name] = preds
    
    # ONLY the XGB predictions (no original features)
    signal_df = pd.DataFrame()
    for bg_name, preds in signal_predictions.items():
        signal_df[f'BDTG_{bg_name}'] = preds
    
    signal_df.to_csv('signal_predictions.csv', index=False)
    predictions['signal'] = signal_df
    print(f"  Created prediction file for signal")
    
    for bg_name_test, bg_test in background_tests.items():
        bg_predictions = {}
        
        for bg_name_model, model_info in models.items():
            model = model_info['model']
            
            scaler = datasets[bg_name_model]['scaler']
            
            selected_features = model_info['features']
            
            X_bg = bg_test[feature_names]
            X_bg_scaled = pd.DataFrame(
                scaler.transform(X_bg),
                columns=X_bg.columns,
                index=X_bg.index
            )
            
            X_bg_selected = X_bg_scaled[selected_features]
            
            preds = model.predict_proba(X_bg_selected)[:, 1]
            
            bg_predictions[bg_name_model] = preds
        
        bg_df = pd.DataFrame()
        for bg_name_model, preds in bg_predictions.items():
            bg_df[f'BDTG_{bg_name_model}'] = preds
        
        bg_df.to_csv(f'{bg_name_test}_predictions.csv', index=False)
        predictions[bg_name_test] = bg_df
        
        print(f"  Created prediction file for {bg_name_test}")
    
    print("All prediction files created")
    return predictions

def calculate_significance(signal_count, background_count):
    if signal_count + background_count <= 0:
        return 0
    return signal_count / np.sqrt(signal_count + background_count)

def main(run_ga=False):
    signal_file = "signal.csv"
    background_files = [
        "Bqq.csv", "Btt.csv", "BZZ.csv", "BWW.csv", 
        "BqqX.csv", "BqqqqX.csv", "BqqHX.csv", "BZH.csv", 
        "Bpebb.csv", "Bpebbqq.csv", "BpeqqH.csv", "Bpett.csv"
    ]
    
    weights = {
        'signal': 0.0015552*1.155,  # weightHH
        'Bqq': 0.0349,              # weightqq
        'Btt': 0.503,               # weightttbar
        'BZZ': 0.17088*1.155,       # weightZZ
        'BWW': 0.5149,              # weightWW
        'BqqX': 0.04347826,         # weightqqX
        'BqqqqX': 0.04,             # weightqqqqX
        'BqqHX': 0.001,             # weightqqHX
        'BZH': 0.00207445*1.155,    # weightZH
        'Bpebb': 0.7536,            # weightpebb
        'Bpebbqq': 0.1522,          # weightpebbqq
        'BpeqqH': 0.1237,           # weightpeqqH
        'Bpett': 0.0570             # weightpett
    }
    
    # load and prepare data
    print("Loading and preparing data...")
    datasets, signal_test, background_tests = load_data(signal_file, background_files)
    
    # train XGB models (originally BDTG so naming matches that)
    print("\nTraining XGB models...")
    models = train_bdtg_models(datasets)
    
    # apply to all event topologies and create prediction files
    print("\nApplying models and creating prediction files...")
    predictions = apply_models(models, datasets, signal_test, background_tests)
    
    print("\nXGB analysis complete! All prediction files have been saved.")
    
    if run_ga:
        run_genetic_algorithm()

def run_genetic_algorithm():

    print("\nStarting Genetic Algorithm optimization...")
    try:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        ga_script = os.path.join(current_dir, "GA.py")
        
        subprocess.run([sys.executable, ga_script], check=True)
        print("\nGenetic Algorithm optimization completed successfully!")
    except subprocess.CalledProcessError as e:
        print(f"\nError running Genetic Algorithm: {e}")
    except Exception as e:
        print(f"\nUnexpected error: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run BDTG analysis with optional GA optimization')
    parser.add_argument('--run-ga', action='store_true', 
                        help='Run GA.py after generating BDTG predictions')
    
    args = parser.parse_args()
    
    main(run_ga=args.run_ga)