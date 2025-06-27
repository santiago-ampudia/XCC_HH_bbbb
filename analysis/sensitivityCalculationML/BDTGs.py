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

# Set random seeds for reproducibility
np.random.seed(42)

def load_data(signal_file, background_files, test_size=0.25):
    """
    Load and prepare data from CSV files for training and testing
    
    Args:
        signal_file (str): Path to the signal CSV file
        background_files (list): List of paths to background CSV files
        test_size (float): Proportion of data to use for testing (default: 0.25)
        
    Returns:
        dict: Dictionary containing datasets for each background type
    """
    # Load signal data from CSV file
    signal_data = pd.read_csv(signal_file)
    signal_data['label'] = 1  # Add label column with value 1 (signal)
    
    print(f"\nSignal file: {os.path.basename(signal_file)}")
    print(f"  Total signal events: {len(signal_data)}")
    
    # Get feature names (all columns except label)
    feature_names = [col for col in signal_data.columns if col != 'label']
    
    # Split signal data into training and testing sets with stratification
    # Use stratified sampling to maintain class distribution
    signal_train, signal_test = train_test_split(
        signal_data, test_size=test_size, random_state=42, stratify=signal_data['label']
    )
    
    print(f"  Signal events for training: {len(signal_train)}")
    print(f"  Signal events for testing: {len(signal_test)}")
    
    # Save signal training and testing sets to CSV files
    signal_train.to_csv('signal_train.csv', index=False)
    signal_test.to_csv('signal_test.csv', index=False)
    print(f"  Saved signal training events to signal_train.csv")
    print(f"  Saved signal testing events to signal_test.csv")
    
    # Dictionary to store individual background dataframes
    background_dfs = {}
    
    # Load each background file
    print("\nBackground files:")
    for bg_file in background_files:
        bg_name = os.path.splitext(os.path.basename(bg_file))[0]
        try:
            bg = pd.read_csv(bg_file)
            bg['label'] = 0  # Add label column with value 0 (background)
            
            # Store in the dictionary
            background_dfs[bg_name] = bg
            
            # Print background data information
            print(f"  {bg_name}: {len(bg)} events")
        except FileNotFoundError:
            print(f"  Warning: {bg_name} file not found, skipping")
    
    # Create a dictionary to store datasets for each background type
    datasets = {}
    
    # For each background type, create a dataset with signal vs that background
    for bg_name, bg_data in background_dfs.items():
        # Split background data into training and testing sets with stratification
        bg_train, bg_test = train_test_split(
            bg_data, test_size=test_size, random_state=42, stratify=bg_data['label']
        )
        
        # Save background training and testing sets to CSV files
        bg_train.to_csv(f'{bg_name}_train.csv', index=False)
        bg_test.to_csv(f'{bg_name}_test.csv', index=False)
        print(f"  Saved {bg_name} training events to {bg_name}_train.csv")
        print(f"  Saved {bg_name} testing events to {bg_name}_test.csv")
        
        # Combine with the same signal training and testing sets
        X_train = pd.concat([signal_train[feature_names], bg_train[feature_names]])
        y_train = pd.concat([signal_train['label'], bg_train['label']])
        
        X_test = pd.concat([signal_test[feature_names], bg_test[feature_names]])
        y_test = pd.concat([signal_test['label'], bg_test['label']])
        
        # Calculate and log class distribution for this dataset
        n_signal_train = len(signal_train)
        n_background_train = len(bg_train)
        n_signal_test = len(signal_test)
        n_background_test = len(bg_test)
        
        # Calculate class imbalance ratio
        imbalance_ratio = n_background_train / n_signal_train if n_signal_train > 0 else float('inf')
        
        print(f"\n  Class distribution for {bg_name} vs Signal:")
        print(f"    Training - Signal: {n_signal_train}, Background: {n_background_train}")
        print(f"    Testing - Signal: {n_signal_test}, Background: {n_background_test}")
        print(f"    Imbalance ratio (Background/Signal): {imbalance_ratio:.2f}")
        
        # Create a StandardScaler for feature normalization
        scaler = StandardScaler()
        
        # Fit the scaler on training data and transform both training and test data
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
        
        # Store in datasets dictionary with additional class imbalance information
        datasets[bg_name] = {
            'X_train': X_train_scaled,
            'y_train': y_train.values,
            'X_test': X_test_scaled,
            'y_test': y_test.values,
            'scaler': scaler,
            'feature_names': feature_names,
            'imbalance_ratio': imbalance_ratio,
            'n_signal_train': n_signal_train,
            'n_background_train': n_background_train
        }
        
        # Save combined (signal vs background) training and testing sets to CSV files
        # Create DataFrames with features and labels
        train_df = X_train.copy()
        train_df['label'] = y_train.values
        test_df = X_test.copy()
        test_df['label'] = y_test.values
        
        # Save to CSV files
        train_df.to_csv(f'signal_vs_{bg_name}_train.csv', index=False)
        test_df.to_csv(f'signal_vs_{bg_name}_test.csv', index=False)
        print(f"  Saved combined signal vs {bg_name} training events to signal_vs_{bg_name}_train.csv")
        print(f"  Saved combined signal vs {bg_name} testing events to signal_vs_{bg_name}_test.csv")
        
        # Print information about this dataset
        print(f"\nDataset for {bg_name} vs Signal:")
        print(f"  Training: {len(X_train_scaled)} events ({len(signal_train)} signal, {len(bg_train)} background)")
        print(f"  Testing: {len(X_test_scaled)} events ({len(signal_test)} signal, {len(bg_test)} background)")
    
    return datasets, signal_test, {bg_name: train_test_split(bg, test_size=test_size, random_state=42)[1] 
                                  for bg_name, bg in background_dfs.items()}

def train_bdtg_models(datasets, n_estimators=200, learning_rate=0.1, max_depth=5):
    """
    Train XGBoost models for each background vs signal dataset with optimized parameters
    
    Args:
        datasets (dict): Dictionary containing datasets for each background type
        n_estimators (int): Maximum number of boosting stages/trees
        learning_rate (float): Learning rate shrinks the contribution of each tree
        max_depth (int): Maximum depth of the individual regression estimators
        
    Returns:
        dict: Dictionary containing trained models for each background type
    """
    # Dictionary to store trained models
    models = {}
    
    # Train a separate model for each background type
    print("\nTraining XGBoost models for each background type vs signal:")
    for bg_name, data in tqdm(datasets.items(), desc="Background types"):
        print(f"\n  Training XGBoost for {bg_name} vs Signal")
        
        # Get the data
        X_train = data['X_train']
        y_train = data['y_train']
        X_test = data['X_test']
        y_test = data['y_test']
        
        # Get class imbalance information
        imbalance_ratio = data['imbalance_ratio']
        n_signal = data['n_signal_train']
        n_background = data['n_background_train']
        
        print(f"  Dataset info: {n_signal} signal, {n_background} background events")
        print(f"  Class imbalance ratio: {imbalance_ratio:.2f}")
        
        # Create validation set for early stopping with stratified sampling
        X_train_main, X_valid, y_train_main, y_valid = train_test_split(
            X_train, y_train, test_size=0.2, random_state=42, stratify=y_train
        )
        
        # Feature importance-based selection with optimized parameters
        feature_selector = xgb.XGBClassifier(
            n_estimators=50,  # Use fewer trees for feature selection
            learning_rate=0.1,
            max_depth=3,
            random_state=42,
            n_jobs=-1,
            tree_method='hist',
            subsample=0.8,  # Use 80% of data per tree for better generalization
            colsample_bytree=0.8  # Use 80% of features per tree for better generalization
        )
        
        # Train a quick model to get feature importances
        feature_selector.fit(X_train_main, y_train_main)
        
        # Get feature importances
        importances = feature_selector.feature_importances_
        feature_names = X_train.columns.tolist()
        
        # Sort features by importance
        sorted_idx = np.argsort(importances)[::-1]
        
        # Select features more intelligently - use top features that contribute to 80% importance
        cumulative_importance = 0
        important_features = []
        
        for idx in sorted_idx:
            important_features.append(feature_names[idx])
            cumulative_importance += importances[idx]
            if cumulative_importance >= 0.8:  # Keep features that contribute to 80% of importance
                break
        
        print(f"  Selected {len(important_features)} out of {len(feature_names)} features")
        
        # Print features in order of importance with their importance scores
        print(f"\n  Features for {bg_name} model in order of importance:")
        for i, idx in enumerate(sorted_idx[:len(important_features)]):
            feature = feature_names[idx]
            importance = importances[idx]
            print(f"    {i+1}. {feature}: {importance:.4f}")
        print()  # Add an empty line for better readability
        
        # Use only important features
        X_train_main_selected = X_train_main[important_features]
        X_valid_selected = X_valid[important_features]
        X_test_selected = X_test[important_features]
        
        # Check if GPU is available
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
        
        # Calculate scale_pos_weight for class imbalance handling
        # This is the ratio of negative (background) to positive (signal) instances
        scale_pos_weight = np.sum(y_train_main == 0) / np.sum(y_train_main == 1)
        
        print(f"  Calculated scale_pos_weight for class imbalance: {scale_pos_weight:.4f}")
        print(f"  This means background events are weighted {scale_pos_weight:.2f}x relative to signal events")
        
        # For extremely imbalanced datasets (common in rare signal analysis), 
        # we might want to cap the scale_pos_weight to prevent over-weighting
        max_scale_weight = 100.0  # Cap to prevent extreme weighting
        if scale_pos_weight > max_scale_weight:
            print(f"  Warning: Very high imbalance detected. Capping scale_pos_weight at {max_scale_weight}")
            scale_pos_weight = max_scale_weight
        
        # Initialize XGBoost model with highly optimized parameters for rare signal analysis
        model = xgb.XGBClassifier(
            n_estimators=n_estimators,
            learning_rate=learning_rate,
            max_depth=max_depth,
            min_child_weight=3,
            gamma=0.1,
            subsample=0.8,
            colsample_bytree=0.8,
            colsample_bylevel=0.8,  # Sample columns at each level for better performance
            colsample_bynode=0.8,   # Sample columns at each node for better performance
            reg_alpha=0.01,         # L1 regularization
            reg_lambda=1,           # L2 regularization
            scale_pos_weight=scale_pos_weight,  # Handle class imbalance - critical for rare signal analysis
            random_state=42,
            tree_method=tree_method,
            predictor=predictor,    # Use GPU predictor if available
            n_jobs=-1,              # Use all CPU cores
            grow_policy='lossguide', # More efficient tree growth policy
            max_bin=256,            # Optimize histogram bins
            early_stopping_rounds=50, # More patient early stopping
            eval_metric=['auc', 'logloss']  # Track multiple metrics
        )
        
        # Create evaluation set for early stopping
        eval_set = [(X_train_main_selected, y_train_main), (X_valid_selected, y_valid)]
        
        # Train the model with early stopping
        print("  Training model with early stopping and class imbalance handling...")
        model.fit(
            X_train_main_selected, 
            y_train_main,
            eval_set=eval_set,
            verbose=100  # Print progress every 100 iterations
        )
        
        # Get the best iteration
        best_iteration = model.best_iteration
        if best_iteration is not None:
            print(f"  Best iteration: {best_iteration}")
        else:
            print("  Early stopping not triggered, used all iterations")
        
        # Evaluate on test set
        y_pred = model.predict_proba(X_test_selected)[:, 1]
        auc = roc_auc_score(y_test, y_pred)
        
        print(f"  AUC for {bg_name}: {auc:.4f}")
        
        # Store model and selected features in dictionary
        models[bg_name] = {
            'model': model,
            'auc': auc,
            'features': important_features,
            'scale_pos_weight': scale_pos_weight
        }
    
    return models

def apply_models(models, datasets, signal_test, background_tests):
    """
    Apply trained models to all event topologies and create prediction files
    
    Args:
        models (dict): Dictionary containing trained models for each background type
        datasets (dict): Dictionary containing datasets for each background type
        signal_test (DataFrame): Signal test data
        background_tests (dict): Dictionary containing background test data for each type
        
    Returns:
        dict: Dictionary containing prediction DataFrames for each event topology
    """
    print("\nApplying models to all event topologies and creating prediction files...")
    
    # Get feature names from the first dataset
    first_bg = next(iter(datasets))
    feature_names = datasets[first_bg]['feature_names']
    
    # Create a dictionary to store predictions for each event topology
    predictions = {}
    
    # Apply models to signal test data
    signal_predictions = {}
    for bg_name, model_info in models.items():
        model = model_info['model']
        
        # Get the scaler for this background
        scaler = datasets[bg_name]['scaler']
        
        # Get the selected features for this model
        selected_features = model_info['features']
        
        # Scale signal test data
        X_signal = signal_test[feature_names]
        X_signal_scaled = pd.DataFrame(
            scaler.transform(X_signal),
            columns=X_signal.columns,
            index=X_signal.index
        )
        
        # Select only the features used by this model
        X_signal_selected = X_signal_scaled[selected_features]
        
        # Make predictions
        preds = model.predict_proba(X_signal_selected)[:, 1]
        
        # Store predictions
        signal_predictions[bg_name] = preds
    
    # Create a DataFrame with ONLY the BDTG predictions (no original features)
    signal_df = pd.DataFrame()
    for bg_name, preds in signal_predictions.items():
        signal_df[f'BDTG_{bg_name}'] = preds
    
    # Save signal predictions
    signal_df.to_csv('signal_predictions.csv', index=False)
    predictions['signal'] = signal_df
    print(f"  Created prediction file for signal")
    
    # Apply models to all background test data
    for bg_name_test, bg_test in background_tests.items():
        bg_predictions = {}
        
        for bg_name_model, model_info in models.items():
            model = model_info['model']
            
            # Get the scaler for this model
            scaler = datasets[bg_name_model]['scaler']
            
            # Get the selected features for this model
            selected_features = model_info['features']
            
            # Scale background test data
            X_bg = bg_test[feature_names]
            X_bg_scaled = pd.DataFrame(
                scaler.transform(X_bg),
                columns=X_bg.columns,
                index=X_bg.index
            )
            
            # Select only the features used by this model
            X_bg_selected = X_bg_scaled[selected_features]
            
            # Make predictions
            preds = model.predict_proba(X_bg_selected)[:, 1]
            
            # Store predictions
            bg_predictions[bg_name_model] = preds
        
        # Create a DataFrame with ONLY the BDTG predictions (no original features)
        bg_df = pd.DataFrame()
        for bg_name_model, preds in bg_predictions.items():
            bg_df[f'BDTG_{bg_name_model}'] = preds
        
        # Save background predictions
        bg_df.to_csv(f'{bg_name_test}_predictions.csv', index=False)
        predictions[bg_name_test] = bg_df
        
        print(f"  Created prediction file for {bg_name_test}")
    
    # Create plots for each model showing signal vs its associated background
    print("\nCreating output distribution plots for each model...")
    
    # Define color mapping similar to ROOT colors
    color_map = {
        'signal': 'blue',           # kBlue
        'Bqq': '#CC0000',           # kRed+2
        'Btt': '#FF0099',           # kPink+6
        'BZZ': '#33CC00',           # kGreen+1
        'BWW': '#CCCC00',           # kYellow-2
        'BqqX': '#006666',          # kTeal-6
        'BqqqqX': '#CC00CC',        # kMagenta+2
        'BqqHX': '#FF6600',         # kOrange+7
        'BZH': '#FFAA00',           # kOrange
        'Bpebb': '#660099',         # kViolet-4
        'Bpebbqq': '#009999',       # kTeal-5
        'BpeqqH': '#99CC00',        # kSpring-8
        'Bpett': '#FF9933'          # kOrange+1
    }
    
    # Define process labels using matplotlib's math rendering
    process_labels = {
        'signal': r'$\gamma\gamma \rightarrow HH \rightarrow b\bar{b}b\bar{b}$',
        'Bqq': r'$\gamma\gamma \rightarrow q\bar{q}$',
        'Btt': r'$\gamma\gamma \rightarrow t\bar{t}$',
        'BZZ': r'$\gamma\gamma \rightarrow ZZ$',
        'BWW': r'$\gamma\gamma \rightarrow W^+W^-$',
        'BqqX': r'$e\gamma \rightarrow q\bar{q}$',
        'BqqqqX': r'$e\gamma \rightarrow q\bar{q}q\bar{q}$',
        'BqqHX': r'$e\gamma \rightarrow q\bar{q}H$',
        'BZH': r'$\gamma\gamma \rightarrow ZH$',
        'Bpebb': r'$e^+e^- \rightarrow b\bar{b}$',
        'Bpebbqq': r'$e^+e^- \rightarrow b\bar{b}q\bar{q}$',
        'BpeqqH': r'$e^+e^- \rightarrow ZH$',
        'Bpett': r'$e^+e^- \rightarrow t\bar{t}$'
    }
    
    # Use matplotlib's built-in math rendering instead of LaTeX
    plt.rcParams['text.usetex'] = False
    plt.rcParams['mathtext.fontset'] = 'dejavusans'  # Use a modern math font
    
    for bg_name, model_info in models.items():
        # Get predictions for signal and this specific background
        signal_preds = signal_predictions[bg_name]
        bg_preds = background_tests[bg_name][f'BDTG_{bg_name}'] if f'BDTG_{bg_name}' in background_tests[bg_name] else predictions[bg_name][f'BDTG_{bg_name}']
        
        # Create figure
        plt.figure(figsize=(10, 6))
        
        # Plot histograms
        bins = np.linspace(0, 1, 50)
        plt.hist(signal_preds, bins=bins, alpha=0.7, label=process_labels['signal'], color=color_map['signal'], density=True)
        plt.hist(bg_preds, bins=bins, alpha=0.7, label=process_labels[bg_name], color=color_map[bg_name], density=True)
        
        # Add labels with adjusted positions and sizes
        plt.xlabel('XGBoost Output Score', fontsize=10, loc='right')
        plt.ylabel('Normalized Density', fontsize=10, loc='top')
        
        # No title
        # plt.title(f'XGBoost Output Distribution: Signal vs {bg_name}', fontsize=16)
        
        plt.legend(fontsize=12)
        plt.grid(alpha=0.3)
        
        # Save figure
        plot_filename = f'xgboost_output_signal_vs_{bg_name}.png'
        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved plot to {plot_filename}")
    
    print("All prediction files and plots created successfully!")
    return predictions

def calculate_significance(signal_count, background_count):
    """
    Calculate signal significance as signal/sqrt(signal+background)
    
    Args:
        signal_count (float): Number of signal events
        background_count (float): Number of background events
        
    Returns:
        float: Signal significance
    """
    if signal_count + background_count <= 0:
        return 0
    return signal_count / np.sqrt(signal_count + background_count)

def main(run_ga=False):
    """
    Main function to run the entire analysis pipeline
    
    Args:
        run_ga (bool): Whether to run GA.py after generating predictions
    """
    # Define file paths
    signal_file = "signal.csv"
    background_files = [
        "Bqq.csv", "Btt.csv", "BZZ.csv", "BWW.csv", 
        "BqqX.csv", "BqqqqX.csv", "BqqHX.csv", "BZH.csv", 
        "Bpebb.csv", "Bpebbqq.csv", "BpeqqH.csv", "Bpett.csv"
    ]
    
    # Initialize weights 
    weights = {
        'signal': 0.001552*1.155,  # weightHH
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
    
    # Load and prepare data
    print("Loading and preparing data...")
    datasets, signal_test, background_tests = load_data(signal_file, background_files)
    
    # Train BDTG models
    print("\nTraining BDTG models...")
    models = train_bdtg_models(datasets)
    
    # Apply models to all event topologies and create prediction files
    print("\nApplying models and creating prediction files...")
    predictions = apply_models(models, datasets, signal_test, background_tests)
    
    print("\nBDTG analysis complete! All prediction files have been saved.")
    
    # Run GA.py if requested
    if run_ga:
        run_genetic_algorithm()

def run_genetic_algorithm():
    """
    Run the Genetic Algorithm optimization using GA.py
    """
    print("\nStarting Genetic Algorithm optimization...")
    try:
        # Get the current directory
        current_dir = os.path.dirname(os.path.abspath(__file__))
        ga_script = os.path.join(current_dir, "GA.py")
        
        # Run GA.py as a subprocess
        subprocess.run([sys.executable, ga_script], check=True)
        print("\nGenetic Algorithm optimization completed successfully!")
    except subprocess.CalledProcessError as e:
        print(f"\nError running Genetic Algorithm: {e}")
    except Exception as e:
        print(f"\nUnexpected error: {e}")

if __name__ == "__main__":
    # Set up command line argument parser
    parser = argparse.ArgumentParser(description='Run BDTG analysis with optional GA optimization')
    parser.add_argument('--run-ga', action='store_true', 
                        help='Run GA.py after generating BDTG predictions')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Run main function with parsed arguments
    main(run_ga=args.run_ga)