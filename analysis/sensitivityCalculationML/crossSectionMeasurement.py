import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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
    
    chi_squared_values = 2 * ((s_values + b) - n * np.log(s_values + b) + n * np.log(b))
    
    min_idx = np.argmin(chi_squared_values)
    min_chi_squared = chi_squared_values[min_idx]
    min_x = x_values[min_idx]
    
    cross_section = min_x / luminosity
    
    print(f"Signal efficiency: {total_effi:.6f}")
    print(f"Branching ratio: {total_BR:.6f}")
    print(f"min ChiSquared: {min_chi_squared:.6f} for x: {min_x:.6f}")
    print(f"Cross section: {cross_section:.6f} fb")
    
    plt.figure(figsize=(10, 8))
    
    plot_step = max(1, len(x_values) // 10000)
    plt.plot(x_values[::plot_step], chi_squared_values[::plot_step], 'b-')
    
    plt.axhline(y=min_chi_squared + 1, color='r', linestyle='-', linewidth=2)
    
    target_chi = min_chi_squared + 1
    
    left_mask = (x_values < min_x) & (np.abs(chi_squared_values - target_chi) < 0.1)
    if np.any(left_mask):
        left_candidates = x_values[left_mask]
        left_chi_values = chi_squared_values[left_mask]
        left_distances = np.abs(left_chi_values - target_chi)
        left_idx = np.argmin(left_distances)
        x_error_left = left_candidates[left_idx]
        distance_min_left = left_distances[left_idx]
    else:
        x_error_left, distance_min_left = find_error_point(min_x, min_chi_squared, True, total_effi, total_BR, n, b)
    
    right_mask = (x_values > min_x) & (np.abs(chi_squared_values - target_chi) < 0.1)
    if np.any(right_mask):
        right_candidates = x_values[right_mask]
        right_chi_values = chi_squared_values[right_mask]
        right_distances = np.abs(right_chi_values - target_chi)
        right_idx = np.argmin(right_distances)
        x_error_right = right_candidates[right_idx]
        distance_min_right = right_distances[right_idx]
    else:
        x_error_right, distance_min_right = find_error_point(min_x, min_chi_squared, False, total_effi, total_BR, n, b)
    
    error_left = abs((x_error_left - min_x) / (luminosity * cross_section))
    error_right = abs((x_error_right - min_x) / (luminosity * cross_section))
    
    print(f"errorLeft: {error_left:.6f} for x: {x_error_left:.6f} with distance: {distance_min_left:.6f}")
    print(f"errorRight: {error_right:.6f} for x: {x_error_right:.6f} with distance: {distance_min_right:.6f}")
    
    plt.axvline(x=x_error_left, color='r', linestyle='--', alpha=0.7)
    plt.axvline(x=x_error_right, color='r', linestyle='--', alpha=0.7)
    plt.axvline(x=min_x, color='g', linestyle='-', alpha=0.7)
    
    plt.xlabel('cross-section*luminosity', fontsize=12)
    plt.ylabel('$\\chi^{2}$', fontsize=14)
    plt.title('Chi-Squared vs Cross Section', fontsize=14)
    
    y_range = np.percentile(chi_squared_values[np.isfinite(chi_squared_values)], [5, 95])
    y_min = max(min_chi_squared - 2, y_range[0])
    y_max = min(min_chi_squared + 10, y_range[1])
    plt.ylim(y_min, y_max)
    
    legend_text = [
        "$\\chi^{2} \\equiv -2\\ln\\frac{L_{s+b}}{L_{b}}$",
        "$L_{s+b} = \\prod_{i} \\frac{e^{-(s_{i}+b_{i})} (s_{i} + b_{i})^{n_{i}}}{n_{i}!}$",
        "$L_{b} = \\prod_{i} \\frac{e^{-b_{i}} b_{i}^{n_{i}}}{n_{i}!}$",
        "uncertainty"
    ]
    
    legend_elements = [
        Line2D([0], [0], color='white', lw=0, label=legend_text[0]),
        Line2D([0], [0], color='white', lw=0, label=legend_text[1]),
        Line2D([0], [0], color='white', lw=0, label=legend_text[2]),
        Line2D([0], [0], color='red', lw=2, label=legend_text[3])
    ]
    
    plt.legend(handles=legend_elements, loc='upper right', fontsize=10, 
              frameon=False, handlelength=0, handletextpad=0)
    
    result_text = f"Cross section: {cross_section:.2f} $\\pm$ {error_right:.2f}/{error_left:.2f} fb"
    plt.annotate(result_text, xy=(0.05, 0.95), xycoords='axes fraction',
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.8),
                 fontsize=12)
    
    os.makedirs('analysis', exist_ok=True)
    plt.savefig('analysis/chiSquared.png', dpi=300, bbox_inches='tight')
    print(f"Plot saved to: analysis/chiSquared.png")
    plt.close()
    
    return {
        'cross_section': float(cross_section),
        'error_top': float(error_right),
        'error_bottom': float(error_left),
        'min_chi_squared': float(min_chi_squared),
        'min_x': float(min_x)
    }

def find_error_point(min_x, min_chi_squared, is_left, total_effi, total_BR, n, b):
   
    target_chi = min_chi_squared + 1
    x_error = min_x
    distance_min = float('inf')
    
    step = 0.01
    
    if is_left:
        search_range = np.arange(min_x - 500, min_x, step)
    else:
        search_range = np.arange(min_x, min_x + 500, step)
    
    for x in search_range:
        s = x * total_effi * total_BR
        chi_squared = 2 * ((s + b) - n * math.log(s + b) + n * math.log(b))
        
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
    parser.add_argument('--output', type=str, default='analysis/chiSquared.png',
                        help='Output file path for the chi-squared plot (default: analysis/chiSquared.png)')
    return parser.parse_args()

# example usage:
if __name__ == "__main__":
    args = parse_arguments()
    
    result = find_cross_section_HHbbbb(args.hh, args.back)
    
    print(f"\nFinal Result:")
    print(f"Cross section: {result['cross_section']:.4f} ± {result['error_top']:.4f}/{result['error_bottom']:.4f} fb")