import os
import sys
import subprocess
import pandas as pd
import numpy as np
import time
import argparse
from datetime import datetime

def run_bdtg_training():
    print("\n=== Running BDTG Training ===")
    try:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        bdtg_script = os.path.join(current_dir, "BDTGs.py")
        
        subprocess.run([sys.executable, bdtg_script], check=True)
        print("BDTG training completed successfully!")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running BDTG training: {e}")
        return False
    except Exception as e:
        print(f"Unexpected error during BDTG training: {e}")
        return False

def run_ga_optimization():
    print("\n=== Running Genetic Algorithm Optimization ===")
    try:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        ga_script = os.path.join(current_dir, "GA.py")
        
        subprocess.run([sys.executable, ga_script], check=True)
        
        results_file = os.path.join(current_dir, "ga_results.txt")
        
        if not os.path.exists(results_file):
            print("GA results file not found!")
            return None
        
        with open(results_file, 'r') as f:
            content = f.read()
        
        significance_line = content.split('\n')[0]
        significance = float(significance_line.split(': ')[1])
        
        cross_section_data = {}
        if "=== CROSS SECTION MEASUREMENT ===" in content:
            cross_section_section = content.split("=== CROSS SECTION MEASUREMENT ===")[1].strip()
            lines = cross_section_section.split('\n')
            
            for line in lines:
                if ': ' in line:
                    key, value = line.split(': ', 1)
                    key = key.strip().lower().replace(' ', '_')
                    try:
                        cross_section_data[key] = float(value.split()[0])
                    except ValueError:
                        cross_section_data[key] = value
        
        thresholds_file = os.path.join(current_dir, "optimal_thresholds.csv")
        thresholds_df = pd.read_csv(thresholds_file)
        
        results = {
            'significance': significance,
            'thresholds': thresholds_df.iloc[0].to_dict(),
            'timestamp': datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        
        if cross_section_data:
            results['cross_section'] = cross_section_data
            print(f"GA optimization completed with significance: {significance:.6f} and cross-section: {cross_section_data.get('cross_section', 'N/A'):.6f} fb")
        else:
            print(f"GA optimization completed with significance: {significance:.6f}")
        
        return results
    
    except subprocess.CalledProcessError as e:
        print(f"Error running GA optimization: {e}")
        return None
    except Exception as e:
        print(f"Unexpected error during GA optimization: {e}")
        return None

def calculate_statistics(results_list):
    if not results_list:
        return None
    
    significances = [result['significance'] for result in results_list]
    
    stats = {
        'significance_mean': np.mean(significances),
        'significance_median': np.median(significances),
        'significance_std': np.std(significances),
        'significance_min': np.min(significances),
        'significance_max': np.max(significances),
        'significance_variance': np.var(significances),
        'count': len(significances)
    }
    
    if any('cross_section' in result for result in results_list):
        xs_params = set()
        for result in results_list:
            if 'cross_section' in result:
                xs_params.update(result['cross_section'].keys())
        
        for param in xs_params:
            param_values = [result['cross_section'].get(param, 0) 
                           for result in results_list 
                           if 'cross_section' in result and param in result['cross_section']]
            
            if param_values:
                param_prefix = f'xs_{param}_'
                stats.update({
                    f'{param_prefix}mean': np.mean(param_values),
                    f'{param_prefix}median': np.median(param_values),
                    f'{param_prefix}std': np.std(param_values),
                    f'{param_prefix}min': np.min(param_values),
                    f'{param_prefix}max': np.max(param_values),
                    f'{param_prefix}variance': np.var(param_values)
                })
    
    return stats

def save_results_to_csv(results_list, stats, filename="analysisResults.csv"):
    runs_data = []
    for i, result in enumerate(results_list):
        run_data = {
            'run_number': i + 1,
            'significance': result['significance'],
            'timestamp': result['timestamp']
        }
        for feature, value in result['thresholds'].items():
            run_data[f'threshold_{feature}'] = value
        
        if 'cross_section' in result:
            for key, value in result['cross_section'].items():
                run_data[f'xs_{key}'] = value
        
        runs_data.append(run_data)
    
    runs_df = pd.DataFrame(runs_data)
    
    stats_data = [{
        'statistic': stat_name,
        'value': stat_value
    } for stat_name, stat_value in stats.items()]
    
    stats_df = pd.DataFrame(stats_data)
    
    with open(filename, 'w') as f:
        f.write("=== INDIVIDUAL RUNS ===\n")
        f.write(runs_df.to_csv(index=False))
        f.write("\n=== STATISTICS ===\n")
        f.write(stats_df.to_csv(index=False))
    
    print(f"Results saved to {filename}")

def main():
    parser = argparse.ArgumentParser(description='Run BDTG training and GA optimization multiple times')
    parser.add_argument('--runs', type=int, default=5, help='Number of GA optimization runs (default: 5)')
    args = parser.parse_args()
    
    num_runs = args.runs
    
    print(f"=== Starting Analysis with {num_runs} GA optimization runs ===")
    start_time = time.time()
    
    if not run_bdtg_training():
        print("BDTG training failed. Exiting.")
        return
    
    results_list = []
    for i in range(num_runs):
        print(f"\n--- Starting GA Run {i+1}/{num_runs} ---")
        result = run_ga_optimization()
        if result:
            results_list.append(result)
    
    if results_list:
        stats = calculate_statistics(results_list)
        
        print("\n=== STATISTICS ===")
        for stat_name, stat_value in stats.items():
            print(f"{stat_name}: {stat_value:.6f}" if isinstance(stat_value, float) else f"{stat_name}: {stat_value}")
        
        save_results_to_csv(results_list, stats)
    else:
        print("No valid results obtained from GA runs.")
    
    execution_time = time.time() - start_time
    hours, remainder = divmod(execution_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    print(f"\nTotal execution time: {int(hours)}h {int(minutes)}m {seconds:.2f}s")

if __name__ == "__main__":
    main()