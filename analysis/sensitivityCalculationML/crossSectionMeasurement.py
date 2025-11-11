import numpy as np
import matplotlib.pyplot as plt
import math
import os
import argparse

# fixed luminosity value in fb^-1 for XCC at 380 GeV, 10-years runtime
LUMINOSITY = 4900.0

def find_cross_section_HHbbbb(HH_remaining, back_remaining, luminosity=LUMINOSITY):
    
    total_remaining = HH_remaining + back_remaining
    
    x_min = 50.0
    x_max = 4000.0
    
    branching_ratio_Hbb = 0.5824  
    total_BR = branching_ratio_Hbb * branching_ratio_Hbb
    total_effi = HH_remaining / (1462612 * 0.001225 * total_BR)
    
    b = back_remaining
    n = total_remaining
    
    print(f"b: {b}    n: {n}    n-b: {n-b}")
    print(f"Using luminosity: {luminosity} fb^-1")
    
    print("Starting chi-squared calculation")
    x_values = np.arange(x_min, x_max, 0.01)
    s_values = x_values * total_effi * total_BR
    
    # Calculate ln(n!) using Stirling's approximation for large n, exact for small n
    if n < 100:
        ln_n_factorial = sum(np.log(i) for i in range(1, int(n) + 1))
    else:
        # Stirling's approximation: ln(n!) ≈ n*ln(n) - n
        ln_n_factorial = n * np.log(n) - n
    
    # New chi-squared formula: χ² = 2*(s+b+ln(n!)-n*ln(s+b))
    chi_squared_values_raw = 2 * (s_values + b + ln_n_factorial - n * np.log(s_values + b))
    
    # Find minimum before normalization
    min_idx = np.argmin(chi_squared_values_raw)
    min_chi_squared_raw = chi_squared_values_raw[min_idx]
    min_x = x_values[min_idx]
    
    # Normalize chi-squared so minimum is at 0
    chi_squared_values = chi_squared_values_raw - min_chi_squared_raw
    min_chi_squared = 0.0
    
    cross_section = min_x / luminosity
    
    print(f"Signal efficiency: {total_effi:.6f}")
    print(f"Branching ratio: {total_BR:.6f}")
    print(f"min ChiSquared: {min_chi_squared:.6f} for x: {min_x:.6f}")
    print(f"Cross section: {cross_section:.6f} fb")
    
    # Use ATLAS plotting style
    plt.style.use('classic')
    fig, ax = plt.subplots(figsize=(10, 8))
    
    plot_step = max(1, len(x_values) // 10000)
    ax.plot(x_values[::plot_step], chi_squared_values[::plot_step], 'b-', linewidth=1.5)
    
    ax.axhline(y=1, color='r', linestyle='--', alpha=0.7, linewidth=1.5)
    
    target_chi = 1.0
    
    left_mask = (x_values < min_x) & (np.abs(chi_squared_values - target_chi) < 0.1)
    if np.any(left_mask):
        left_candidates = x_values[left_mask]
        left_chi_values = chi_squared_values[left_mask]
        left_distances = np.abs(left_chi_values - target_chi)
        left_idx = np.argmin(left_distances)
        x_error_left = left_candidates[left_idx]
        distance_min_left = left_distances[left_idx]
    else:
        x_error_left, distance_min_left = find_error_point(min_x, min_chi_squared_raw, True, total_effi, total_BR, n, b, ln_n_factorial)
    
    right_mask = (x_values > min_x) & (np.abs(chi_squared_values - target_chi) < 0.1)
    if np.any(right_mask):
        right_candidates = x_values[right_mask]
        right_chi_values = chi_squared_values[right_mask]
        right_distances = np.abs(right_chi_values - target_chi)
        right_idx = np.argmin(right_distances)
        x_error_right = right_candidates[right_idx]
        distance_min_right = right_distances[right_idx]
    else:
        x_error_right, distance_min_right = find_error_point(min_x, min_chi_squared_raw, False, total_effi, total_BR, n, b, ln_n_factorial)
    
    error_left = abs((x_error_left - min_x) / (luminosity * cross_section))
    error_right = abs((x_error_right - min_x) / (luminosity * cross_section))
    
    print(f"errorLeft: {error_left:.6f} for x: {x_error_left:.6f} with distance: {distance_min_left:.6f}")
    print(f"errorRight: {error_right:.6f} for x: {x_error_right:.6f} with distance: {distance_min_right:.6f}")
    
    ax.axvline(x=x_error_left, color='r', linestyle='--', alpha=0.7, linewidth=1.5)
    ax.axvline(x=x_error_right, color='r', linestyle='--', alpha=0.7, linewidth=1.5)
    
    # Bold axis labels with proper formatting matching other plots
    ax.set_xlabel(r'$\mathbf{\sigma_{HH}} \times \mathbf{\mathcal{L}}$ $\mathbf{(fb^{-1})}$', fontsize=14, fontweight='bold', loc='right')
    ax.set_ylabel(r'$\mathbf{\chi^{2} = -2ln(L_{s+b})}$', fontsize=14, fontweight='bold', loc='top')
    
    # Set axis ranges
    ax.set_xlim(0, 3000)
    
    y_range = np.percentile(chi_squared_values[np.isfinite(chi_squared_values)], [5, 95])
    y_min = 0.0
    y_max = min(10, y_range[1])
    ax.set_ylim(y_min, y_max)
    
    # ATLAS style: increase tick label size
    ax.tick_params(axis='both', which='major', labelsize=14)
    
    # Save plot in the same directory as this script file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    plot_path = os.path.join(script_dir, 'chiSquared.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {plot_path}")
    
    # Display the plot on screen
    plt.show()
    
    return {
        'cross_section': float(cross_section),
        'error_top': float(error_right),
        'error_bottom': float(error_left),
        'min_chi_squared': float(min_chi_squared),
        'min_x': float(min_x)
    }

def find_error_point(min_x, min_chi_squared_raw, is_left, total_effi, total_BR, n, b, ln_n_factorial):
   
    target_chi = 1.0  # After normalization, we look for chi^2 = 1
    x_error = min_x
    distance_min = float('inf')
    
    step = 0.01
    
    if is_left:
        search_range = np.arange(min_x - 500, min_x, step)
    else:
        search_range = np.arange(min_x, min_x + 500, step)
    
    for x in search_range:
        s = x * total_effi * total_BR
        # New chi-squared formula: χ² = 2*(s+b+ln(n!)-n*ln(s+b))
        chi_squared_raw = 2 * (s + b + ln_n_factorial - n * math.log(s + b))
        # Normalize by subtracting minimum
        chi_squared = chi_squared_raw - min_chi_squared_raw
        
        if abs(chi_squared - target_chi) < distance_min:
            distance_min = abs(chi_squared - target_chi)
            x_error = x
            
            if distance_min < 0.001:
                break
    
    return x_error, distance_min

def parse_arguments():
    parser = argparse.ArgumentParser(description='Calculate cross section for HH→bbbb process')
    parser.add_argument('--hh', type=float, required=True, 
                        help='Number of HH signal events remaining after selection')
    parser.add_argument('--back', type=float, required=True,
                        help='Number of background events remaining after selection')
    parser.add_argument('--output', type=str, default='chiSquared.png',
                        help='Output file path for the chi-squared plot (default: chiSquared.png)')
    return parser.parse_args()

# example usage:
if __name__ == "__main__":
    args = parse_arguments()
    
    result = find_cross_section_HHbbbb(args.hh, args.back)
    
    print(f"\nFinal Result:")
    print(f"Cross section: {result['cross_section']:.4f} ± {result['error_top']:.4f}/{result['error_bottom']:.4f} fb")
