from numba import jit
import numpy as np

@jit(nopython=True)
def fast_significance_calculation(thresholds, signal_array, bg_arrays, bg_weights, signal_weight):

    n_signal = signal_array.shape[0]
    n_features = signal_array.shape[1]
    signal_mask = np.ones(n_signal, dtype=np.bool_)
    
    for i in range(n_features):
        signal_mask = signal_mask & (signal_array[:, i] > thresholds[i])
    
    surviving_signal = np.sum(signal_mask) * signal_weight
    
    total_surviving_background = 0.0
    
    for b in range(len(bg_arrays)):
        bg_data = bg_arrays[b]
        n_bg = bg_data.shape[0]
        bg_mask = np.ones(n_bg, dtype=np.bool_)
        
        for i in range(n_features):
            bg_mask = bg_mask & (bg_data[:, i] > thresholds[i])
        
        surviving_bg = np.sum(bg_mask) * bg_weights[b]
        total_surviving_background += surviving_bg
    
    if surviving_signal + total_surviving_background > 0:
        significance = surviving_signal / np.sqrt(surviving_signal + total_surviving_background)
    else:
        significance = 0.0
    
    return significance

_cached_data = None

def prepare_data_for_numba(signal_df, background_dfs, weights):
   
    global _cached_data
    
    if _cached_data is None:
        signal_array = signal_df.values
        
        bg_arrays = []
        bg_weights = []
        
        for bg_type, bg_df in background_dfs.items():
            bg_arrays.append(bg_df.values)
            bg_weights.append(weights[bg_type])
        
        signal_weight = weights['signal']
        
        _cached_data = (signal_array, bg_arrays, np.array(bg_weights), signal_weight)
        print("Data prepared for optimized calculation")
    
    return _cached_data

def optimized_evaluate(individual, signal_df, background_dfs, weights):
    signal_array, bg_arrays, bg_weights, signal_weight = prepare_data_for_numba(signal_df, background_dfs, weights)
    
    significance = fast_significance_calculation(
        np.array(individual, dtype=np.float64), 
        signal_array, 
        bg_arrays, 
        bg_weights, 
        signal_weight
    )
    
    return (significance,)
