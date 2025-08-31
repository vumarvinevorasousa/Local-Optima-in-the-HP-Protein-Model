import os
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import stats
import seaborn as sns
from scipy.stats import norm
from fold_graph import graph_fold
from fold_graph_3D import graph_fold_3D

def plot_distribution(folder, protein, runs, gens, _3D):
    best_stability = 0
    best_folds = []
    protein_length = len(protein)

    for run in range(runs):
        best = []

        folder_run = f'/runs/run{run}'
        if os.path.isdir(f'{folder}/{protein}/{folder_run}'):
            with open(f'{folder}/{protein}/{folder_run}/folds.txt', 'r') as f:
                best = f.readlines()[-1]
            # best = best.splitlines()[0].replace(" ", "").replace("[", "").replace("]", "")
            best = best.splitlines()[0].split("[")[1].replace("]", "").replace(" ", "")

            index = 0
            fold = []
            while index < len(protein) - 1:
                fold.append(best[index])
                index += 1
            free_energy = best[len(protein) - 1:]

            best_folds.append([fold, int(free_energy)])
            if int(free_energy) < best_stability:
                best_stability = int(free_energy)

    # Create a bar plot that shows the distribution of the different free energy values that have been recorded
    energy_bucket = [0 for i in range(abs(best_stability - 1))]
    for fold in best_folds:
        energy_bucket[abs(fold[1])] += 1
    xaxis_energy_bucket = [-abs(i) for i in range(len(energy_bucket))]

    data = []
    for i in range(len(xaxis_energy_bucket)):
        data += [xaxis_energy_bucket[i]]*energy_bucket[i]
    mu, std = norm.fit(data)

    fig, ax = plt.subplots()
    ax = sns.histplot(data=data, bins=len(energy_bucket), discrete=True, label="smth")

    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu, std)
    plt.plot(x, [x*runs for x in p], 'k', linewidth=2)

    plt.gca().invert_xaxis()
    if _3D:
        plt.title(f'[3D] Stability Distribution, Best: {best_stability} (L={len(protein)}, R={runs})')
    else:
        plt.title(f'[2D] Stability Distribution, Best: {best_stability} (L={len(protein)}, R={runs})')
    plt.xlabel("Stability")
    plt.legend(labels=[f'mu={round(mu, 1)}, std={round(std, 2)}'])
    
    plt.savefig(f'{folder}/{protein}/results/bar_energy_and_curve.png')
    plt.close()

def graph_best(folder, protein, runs, gens, _3D):
    best_stability = 0
    best_folds = []
    protein_length = len(protein)

    for run in range(runs):
        best = []

        folder_run = f'/runs/run{run}'
        if os.path.isdir(f'{folder}/{protein}/{folder_run}'):
            with open(f'{folder}/{protein}/{folder_run}/folds.txt', 'r') as f:
                best = f.readlines()[-1]
            # best = best.splitlines()[0].replace(" ", "").replace("[", "").replace("]", "")
            best = best.splitlines()[0].split("[")[1].replace("]", "").replace(" ", "")

            index = 0
            fold = []
            while index < len(protein) - 1:
                fold.append(best[index])
                index += 1
            free_energy = best[len(protein) - 1:]

            best_folds.append([fold, int(free_energy)])
            if int(free_energy) < best_stability:
                best_stability = int(free_energy)

    # Sort the folds according to their free energy and keep track of duplicates; and store these in a file
    fold_bucket = []
    for fold in best_folds:
        if not fold_bucket:
            fold_bucket.append([fold, 1])
        else:
            found = False
            for f in fold_bucket:
                if f[0][0] == fold[0]:
                    found = True
                    f[1] += 1
                    break

            if not found:
                fold_bucket.append([fold, 1])

    fold_energy_bucket = [[] for i in range(abs(best_stability - 1))]
    for fold in fold_bucket:
        fold_energy_bucket[abs(fold[0][1])].append([fold[0][0], fold[1]])

    index = 0
    with open(f'{folder}/{protein}/results/best/best_folds.txt', "w") as f:
        for bucket in fold_energy_bucket:
            f.write(str(index) + "\n")
            for fold in bucket:
                f.write("".join(fold[0]) + " " + str(fold[1]) + "\n")
            
            index -= 1

    # Create a 2D/3D graphical representation for each of the folds
    index = 0
    os.mkdir(f'{folder}/{protein}/results/best/figures')
    for bucket in fold_energy_bucket:
        folder_bucket = f'{folder}/{protein}/results/best/figures/{index}/'
        os.mkdir(folder_bucket)

        file_index = 1
        for fold in bucket:
            if _3D:
                graph_fold_3D(protein, fold[0], folder_bucket + str(file_index) + "_#" + str(fold[1]) + ".png")
            else:
                graph_fold(protein, fold[0], folder_bucket + str(file_index) + "_#" + str(fold[1]) + ".png")
            file_index += 1

        index -= 1
