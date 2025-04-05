import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
from deap import base, creator, tools, algorithms
import multiprocessing
import time
import math
import copy
from numba_optimization import optimized_evaluate, prepare_data_for_numba
import threading
from concurrent.futures import ThreadPoolExecutor

weights = {
    'signal': 0.0015552*1.155 * 4,  
    'Bqq': 0.0349 * 4,              
    'Btt': 0.503 * 4,               
    'BZZ': 0.17088*1.155 * 4,       
    'BWW': 0.5149 * 4,              
    'BqqX': 0.04347826 * 4,         
    'BqqqqX': 0.04 * 4,             
    'BqqHX': 0.001 * 4,             
    'BZH': 0.00207445*1.155 * 4,    
    'Bpebb': 0.7536 * 4,            
    'Bpebbqq': 0.1522 * 4,          
    'BpeqqH': 0.1237 * 4,           
    'Bpett': 0.0570 * 4             
}

def load_data():
    signal_df = pd.read_csv('signal_predictions.csv')
    background_dfs = {}
    background_types = ['Bqq', 'Btt', 'BZZ', 'BWW', 'BqqX', 'BqqqqX', 'BqqHX', 
                       'BZH', 'Bpebb', 'Bpebbqq', 'BpeqqH', 'Bpett']
    
    for bg_type in background_types:
        bg_file = f'{bg_type}_predictions.csv'
        background_dfs[bg_type] = pd.read_csv(bg_file)
    
    return signal_df, background_dfs

def calculate_significance(thresholds, signal_df, background_dfs):
    feature_names = signal_df.columns
    
    signal_mask = np.ones(len(signal_df), dtype=bool)
    for i, feature in enumerate(feature_names):
        signal_mask &= (signal_df[feature] > thresholds[i])
    
    surviving_signal = np.sum(signal_mask) * weights['signal']
    
    total_surviving_background = 0
    for bg_type, bg_df in background_dfs.items():
        bg_mask = np.ones(len(bg_df), dtype=bool)
        for i, feature in enumerate(feature_names):
            bg_mask &= (bg_df[feature] > thresholds[i])
        
        surviving_bg = np.sum(bg_mask) * weights[bg_type]
        total_surviving_background += surviving_bg
    
    if surviving_signal + total_surviving_background > 0:
        significance = surviving_signal / np.sqrt(surviving_signal + total_surviving_background)
    else:
        significance = 0
    
    return significance

def evaluate(individual, signal_df, background_dfs):
    return optimized_evaluate(individual, signal_df, background_dfs, weights)

def create_educated_guesses(signal_df, background_dfs, n_guesses=5):
    feature_names = signal_df.columns
    n_features = len(feature_names)
    guesses = []
    
    feature_stats = {}
    for feature in feature_names:
        signal_mean = signal_df[feature].mean()
        signal_std = signal_df[feature].std()
        bg_means = [bg_df[feature].mean() for bg_type, bg_df in background_dfs.items()]
        
        feature_stats[feature] = {
            'signal_mean': signal_mean,
            'signal_std': signal_std,
            'bg_means': bg_means
        }
    
    for percentile in [50, 60, 70, 80, 90]:
        thresholds = []
        for feature in feature_names:
            threshold = np.percentile(signal_df[feature], percentile)
            thresholds.append(threshold)
        guesses.append(thresholds)
    
    for k in [0.5, 1.0]:
        thresholds = []
        for feature in feature_names:
            stats = feature_stats[feature]
            threshold = max(0, min(1, stats['signal_mean'] - k * stats['signal_std']))
            thresholds.append(threshold)
        guesses.append(thresholds)
    
    thresholds = []
    for feature in feature_names:
        stats = feature_stats[feature]
        avg_bg_mean = np.mean(stats['bg_means'])
        midpoint = (stats['signal_mean'] + avg_bg_mean) / 2
        threshold = max(0, min(1, midpoint))
        thresholds.append(threshold)
    guesses.append(thresholds)
    
    if len(guesses) > n_guesses:
        fitnesses = [calculate_significance(guess, signal_df, background_dfs) for guess in guesses]
        sorted_indices = np.argsort(fitnesses)[::-1]
        guesses = [guesses[i] for i in sorted_indices[:n_guesses]]
    
    return guesses

def adaptive_mutation(individual, mu, sigma, indpb, gen, max_gen, fitness_improvement):
    progress = gen / max_gen
    
    if fitness_improvement > 0.05:
        adaptive_sigma = sigma * (0.5 - 0.3 * progress)
        adaptive_indpb = indpb * (0.5 - 0.3 * progress)
    else:
        adaptive_sigma = sigma * (1.0 - 0.5 * progress)
        adaptive_indpb = indpb * (1.0 - 0.5 * progress)
    
    adaptive_sigma = max(0.05, adaptive_sigma)
    adaptive_indpb = max(0.05, adaptive_indpb)
    
    for i in range(len(individual)):
        if random.random() < adaptive_indpb:
            individual[i] += random.gauss(mu, adaptive_sigma)
            individual[i] = max(0, min(1, individual[i]))
    
    return individual,

def calculate_shared_fitness(population, sigma_share=0.1):
    n = len(population)
    shared_fitnesses = []
    
    for i in range(n):
        ind_i = population[i]
        raw_fitness = ind_i.fitness.values[0]
        
        sharing_factor = 0
        for j in range(n):
            if i == j:
                continue
                
            ind_j = population[j]
            
            distance = math.sqrt(sum((a - b) ** 2 for a, b in zip(ind_i, ind_j)))
            
            if distance < sigma_share:
                sharing_component = 1 - (distance / sigma_share) ** 2
                sharing_factor += sharing_component
        
        shared_fitness = raw_fitness / (1 + sharing_factor)
        shared_fitnesses.append(shared_fitness)
    
    return shared_fitnesses

def setup_genetic_algorithm(signal_df, background_dfs):
    creator.create("FitnessMax", base.Fitness, weights=(1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMax)
    
    toolbox = base.Toolbox()
    
    pool = multiprocessing.Pool()
    toolbox.register("map", pool.map)
    
    toolbox.register("attr_float", random.uniform, 0, 1)
    
    n_features = len(signal_df.columns)
    toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_float, n=n_features)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    
    toolbox.register("evaluate", evaluate, signal_df=signal_df, background_dfs=background_dfs)
    
    toolbox.register("mate", tools.cxBlend, alpha=0.5)
    toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=0.4, indpb=0.4)
    toolbox.register("select", tools.selTournament, tournsize=7)
    
    return toolbox, pool

def eaAdvanced(population, toolbox, cxpb, mutpb, ngen, stats=None,
               halloffame=None, verbose=__debug__, signal_df=None, background_dfs=None):
    num_elites = 15
    stagnation_limit = 15
    restart_diversity_factor = 0.7
    
    logbook = tools.Logbook()
    logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])
    
    invalid_ind = [ind for ind in population if not ind.fitness.valid]
    fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
    for ind, fit in zip(invalid_ind, fitnesses):
        ind.fitness.values = fit
    
    if halloffame is not None:
        halloffame.update(population)
    
    record = stats.compile(population) if stats else {}
    logbook.record(gen=0, nevals=len(invalid_ind), **record)
    if verbose:
        print(logbook.stream)
    
    best_fitness_history = [max([ind.fitness.values[0] for ind in population])]
    best_individual_history = [copy.deepcopy(tools.selBest(population, 1)[0])]
    
    early_stop_threshold = 0.001
    patience = 20
    
    for gen in range(1, ngen + 1):
        if gen >= 5:
            fitness_improvement = (best_fitness_history[-1] - best_fitness_history[-5]) / best_fitness_history[-5]
        else:
            fitness_improvement = 1.0
        
        offspring = toolbox.select(population, len(population) - num_elites)
        
        elites = tools.selBest(population, num_elites)
        
        offspring = list(map(toolbox.clone, offspring))
        
        for i in range(1, len(offspring), 2):
            if random.random() < cxpb:
                offspring[i-1], offspring[i] = toolbox.mate(offspring[i-1], offspring[i])
                del offspring[i-1].fitness.values
                del offspring[i].fitness.values
        
        for i in range(len(offspring)):
            if random.random() < mutpb:
                offspring[i], = adaptive_mutation(offspring[i], 0, 0.4, 0.4, gen, ngen, fitness_improvement)
                del offspring[i].fitness.values
        
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit
        
        offspring.extend(elites)
        
        if gen % 5 == 0:
            shared_fitnesses = calculate_shared_fitness(offspring)
            for ind, shared_fit in zip(offspring, shared_fitnesses):
                ind.fitness.values = (shared_fit,)
        
        if halloffame is not None:
            halloffame.update(offspring)
        
        population[:] = offspring
        
        current_best_fitness = max([ind.fitness.values[0] for ind in population])
        best_fitness_history.append(current_best_fitness)
        best_individual_history.append(copy.deepcopy(tools.selBest(population, 1)[0]))
        
        record = stats.compile(population) if stats else {}
        logbook.record(gen=gen, nevals=len(invalid_ind), **record)
        if verbose:
            print(logbook.stream)
        
        if gen > stagnation_limit:
            recent_best = max(best_fitness_history[-stagnation_limit:])
            if best_fitness_history[-1] <= recent_best and best_fitness_history[-1] == best_fitness_history[-stagnation_limit]:
                print(f"\nGeneration {gen}: Stagnation detected! Performing partial restart...")
                
                keep_count = int(len(population) * (1 - restart_diversity_factor))
                best_individuals = tools.selBest(population, keep_count)
                
                new_individuals = [toolbox.individual() for _ in range(len(population) - keep_count)]
                
                fitnesses = toolbox.map(toolbox.evaluate, new_individuals)
                for ind, fit in zip(new_individuals, fitnesses):
                    ind.fitness.values = fit
                
                population[:] = best_individuals + new_individuals
                
                mutpb = 0.5
                print(f"Restart complete. New mutation rate: {mutpb}")
        
        if gen > patience:
            improvement = (best_fitness_history[-1] - best_fitness_history[-patience]) / best_fitness_history[-patience]
            if improvement < early_stop_threshold:
                print(f"\nEarly stopping at generation {gen}: Improvement below threshold ({improvement:.6f} < {early_stop_threshold})")
                break
    
    best_gen_idx = best_fitness_history.index(max(best_fitness_history))
    best_individual = best_individual_history[best_gen_idx]
    
    if best_gen_idx < len(best_individual_history) - 1:
        print(f"\nBest individual found at generation {best_gen_idx}, not in final population. Restoring best individual.")
        worst_idx = min(range(len(population)), key=lambda i: population[i].fitness.values[0])
        population[worst_idx] = toolbox.clone(best_individual)
    
    return population, logbook

def island_model(signal_df, background_dfs, n_islands=5, n_migrations=5, 
                 island_size=100, n_gen_per_migration=40, total_gen=100):
    print("Starting Parallel Island Model optimization...")
    
    base_toolbox, pool = setup_genetic_algorithm(signal_df, background_dfs)
    
    print("Generating educated initial guesses...")
    educated_guesses = create_educated_guesses(signal_df, background_dfs, n_guesses=min(5, island_size//2))
    
    islands = []
    for i in range(n_islands):
        random_pop = base_toolbox.population(n=island_size)
        
        if educated_guesses and i == 0:
            for j, guess in enumerate(educated_guesses):
                if j < len(random_pop):
                    for k, val in enumerate(guess):
                        random_pop[j][k] = val
        
        islands.append(random_pop)
    
    island_hofs = [tools.HallOfFame(1) for _ in range(n_islands)]
    
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("std", np.std)
    stats.register("min", np.min)
    stats.register("max", np.max)
    
    global_hof = tools.HallOfFame(1)
    
    for migration in range(n_migrations):
        print(f"\nMigration cycle {migration+1}/{n_migrations}")
        
        with ThreadPoolExecutor(max_workers=min(n_islands, multiprocessing.cpu_count())) as executor:
            futures = []
            for i in range(n_islands):
                print(f"Submitting Island {i+1}/{n_islands} for evolution")
                future = executor.submit(
                    eaAdvanced,
                    islands[i],
                    base_toolbox,
                    0.7,  # cxpb
                    0.3,  # mutpb
                    n_gen_per_migration,
                    stats,
                    island_hofs[i],
                    True,  # verbose
                    signal_df,
                    background_dfs
                )
                futures.append((i, future))
            
            for i, future in futures:
                islands[i], _ = future.result()
                global_hof.update(islands[i])
                print(f"Island {i+1}/{n_islands} evolution completed")
        
        if migration < n_migrations - 1:
            print("\nPerforming migration between islands...")
            
            diversity_measures = []
            for island in islands:
                total_dist = 0
                count = 0
                for i in range(len(island)):
                    for j in range(i+1, len(island)):
                        dist = math.sqrt(sum((a - b) ** 2 for a, b in zip(island[i], island[j])))
                        total_dist += dist
                        count += 1
                avg_dist = total_dist / count if count > 0 else 0
                diversity_measures.append(avg_dist)
            
            base_migration_rate = 0.1
            migration_rates = []
            for div in diversity_measures:
                rate = base_migration_rate * (1 + (1 - div))
                migration_rates.append(min(0.3, max(0.05, rate)))
            
            for i in range(n_islands):
                source = i
                dest = (i + 1) % n_islands
                
                n_migrants = max(1, int(island_size * migration_rates[source]))
                print(f"  Island {source+1} → Island {dest+1}: {n_migrants} migrants (rate: {migration_rates[source]:.2f})")
                
                migrants = tools.selBest(islands[source], n_migrants)
                
                worst_indices = sorted(range(len(islands[dest])), 
                                       key=lambda j: islands[dest][j].fitness.values[0])[:n_migrants]
                
                for j, worst_idx in enumerate(worst_indices):
                    islands[dest][worst_idx] = base_toolbox.clone(migrants[j])
    
    combined_population = []
    for island in islands:
        combined_population.extend(island)
    
    global_hof.update(combined_population)
    
    pool.close()
    
    best_individual = global_hof[0]
    best_fitness = best_individual.fitness.values[0]
    
    print(f"\nBest fitness across all islands: {best_fitness:.6f}")
    
    return best_individual, best_fitness, None

def calculate_event_statistics(best_thresholds, signal_df, background_dfs):
    feature_names = signal_df.columns
    results = {}
    
    signal_mask = np.ones(len(signal_df), dtype=bool)
    for i, feature in enumerate(feature_names):
        signal_mask &= (signal_df[feature] > best_thresholds[i])
    
    results['signal'] = {
        'initial_events': len(signal_df),
        'surviving_events': np.sum(signal_mask),
        'initial_weighted': len(signal_df) * weights['signal'],
        'surviving_weighted': np.sum(signal_mask) * weights['signal']
    }
    
    for bg_type, bg_df in background_dfs.items():
        bg_mask = np.ones(len(bg_df), dtype=bool)
        for i, feature in enumerate(feature_names):
            bg_mask &= (bg_df[feature] > best_thresholds[i])
        
        results[bg_type] = {
            'initial_events': len(bg_df),
            'surviving_events': np.sum(bg_mask),
            'initial_weighted': len(bg_df) * weights[bg_type],
            'surviving_weighted': np.sum(bg_mask) * weights[bg_type]
        }
    
    return results

def plot_evolution(logbook):
    gen = logbook.select("gen")
    fit_max = logbook.select("max")
    fit_avg = logbook.select("avg")
    
    plt.figure(figsize=(12, 8))
    plt.plot(gen, fit_max, 'b-', label='Maximum Fitness')
    plt.plot(gen, fit_avg, 'r-', label='Average Fitness')
    plt.xlabel('Generation')
    plt.ylabel('Fitness (Significance)')
    plt.title('Evolution of Fitness')
    plt.legend(loc='best')
    plt.grid(True)
    plt.savefig('fitness_evolution.png')
    plt.close()

def main():
    start_time = time.time()
    
    print("Loading data...")
    signal_df, background_dfs = load_data()
    
    print("\nRunning Island Model Genetic Algorithm...")
    best_individual, best_significance, _ = island_model(
        signal_df, background_dfs, 
        n_islands=5,
        n_migrations=5,
        island_size=100,
        n_gen_per_migration=40,
        total_gen=100
    )
    
    print("\nCalculating event statistics...")
    event_stats = calculate_event_statistics(best_individual, signal_df, background_dfs)
    
    print("\n=== GENETIC ALGORITHM RESULTS ===")
    print(f"Best Significance: {best_significance:.6f}")
    print("\nOptimal Thresholds:")
    for i, feature in enumerate(signal_df.columns):
        print(f"{feature}: {best_individual[i]:.6f}")
    
    print("\nEvent Statistics:")
    for topology, stats in event_stats.items():
        print(f"\n{topology}:")
        print(f"  Initial events: {stats['initial_events']}")
        print(f"  Surviving events: {stats['surviving_events']} ({stats['surviving_events']/stats['initial_events']*100:.2f}%)")
        print(f"  Initial weighted: {stats['initial_weighted']:.6f}")
        print(f"  Surviving weighted: {stats['surviving_weighted']:.6f} ({stats['surviving_weighted']/stats['initial_weighted']*100:.2f}%)")
    
    total_initial_bg = sum(stats['initial_weighted'] for topology, stats in event_stats.items() if topology != 'signal')
    total_surviving_bg = sum(stats['surviving_weighted'] for topology, stats in event_stats.items() if topology != 'signal')
    
    print("\nTotal Background:")
    print(f"  Initial weighted: {total_initial_bg:.6f}")
    print(f"  Surviving weighted: {total_surviving_bg:.6f} ({total_surviving_bg/total_initial_bg*100:.2f}%)")
    
    execution_time = time.time() - start_time
    hours, remainder = divmod(execution_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    print(f"\nExecution time: {int(hours)}h {int(minutes)}m {seconds:.2f}s")
    
    with open('ga_results.txt', 'w') as f:
        f.write(f"Best Significance: {best_significance:.6f}\n\n")
        f.write("Optimal Thresholds:\n")
        for i, feature in enumerate(signal_df.columns):
            f.write(f"{feature}: {best_individual[i]:.6f}\n")
        
        f.write("\nEvent Statistics:\n")
        for topology, stats in event_stats.items():
            f.write(f"\n{topology}:\n")
            f.write(f"  Initial events: {stats['initial_events']}\n")
            f.write(f"  Surviving events: {stats['surviving_events']} ({stats['surviving_events']/stats['initial_events']*100:.2f}%)\n")
            f.write(f"  Initial weighted: {stats['initial_weighted']:.6f}\n")
            f.write(f"  Surviving weighted: {stats['surviving_weighted']:.6f} ({stats['surviving_weighted']/stats['initial_weighted']*100:.2f}%)\n")
        
        f.write("\nTotal Background:\n")
        f.write(f"  Initial weighted: {total_initial_bg:.6f}\n")
        f.write(f"  Surviving weighted: {total_surviving_bg:.6f} ({total_surviving_bg/total_initial_bg*100:.2f}%)\n")
        f.write(f"\nExecution time: {int(hours)}h {int(minutes)}m {seconds:.2f}s\n")
    
    with open('optimal_thresholds.csv', 'w') as f:
        f.write(','.join(signal_df.columns) + '\n')
        f.write(','.join(f"{threshold:.6f}" for threshold in best_individual))
    
    print("\nOptimal thresholds saved to 'optimal_thresholds.csv'")
    
    try:
        print("\n=== CROSS SECTION MEASUREMENT ===")
        signal_surviving = event_stats['signal']['surviving_weighted']
        background_surviving = total_surviving_bg
        
        print(f"Running cross-section measurement with:")
        print(f"  Signal events: {signal_surviving:.6f}")
        print(f"  Background events: {background_surviving:.6f}")
        
        import crossSectionMeasurement
        
        cross_section_result = crossSectionMeasurement.find_cross_section_HHbbbb(
            signal_surviving, background_surviving)
        
        with open('ga_results.txt', 'a') as f:
            f.write("\n\n=== CROSS SECTION MEASUREMENT ===\n")
            f.write(f"Signal events: {signal_surviving:.6f}\n")
            f.write(f"Background events: {background_surviving:.6f}\n")
            f.write(f"Cross section: {cross_section_result['cross_section']:.6f} fb\n")
            f.write(f"Error (top): {cross_section_result['error_top']:.6f}\n")
            f.write(f"Error (bottom): {cross_section_result['error_bottom']:.6f}\n")
        
        print("\nCross-section measurement completed and results saved to 'ga_results.txt'")
        print(f"Cross-section plot saved to 'analysis/chiSquared.png'")
        
    except Exception as e:
        print(f"\nError running cross-section measurement: {e}")
        print("Please ensure crossSectionMeasurement.py is in the same directory.")

if __name__ == "__main__":
    main()