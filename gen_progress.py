import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np 
import math
import datetime
import gc
from fold_graph_3D import get_coords_3D

def get_coords(fold):
    pos = [[0,0]]
    for dir in fold:
        match dir:
            case "u":
                pos.append([pos[-1][0], pos[-1][1] + 1])
            case "d":
                pos.append([pos[-1][0], pos[-1][1] - 1])
            case "l":
                pos.append([pos[-1][0] - 1, pos[-1][1]])
            case "r":
                pos.append([pos[-1][0] + 1, pos[-1][1]])
    
    return pos

def free_energy_edges(protein, coords, _3D):
    edges = []
    i = 0
    while i < len(protein):
        if protein[i] == "H":
            j = i + 1
            while j < len(protein):
                if protein[j] == "H":
                    if _3D:
                        if (((coords[i][0] + 1 == coords[j][0] or coords[i][0] - 1 == coords[j][0]) and coords[i][1] == coords[j][1] and coords[i][2] == coords[j][2]
                        or (coords[i][1] + 1 == coords[j][1] or coords[i][1] - 1 == coords[j][1]) and coords[i][0] == coords[j][0] and coords[i][2] == coords[j][2]
                        or (coords[i][2] + 1 == coords[j][2] or coords[i][2] - 1 == coords[j][2]) and coords[i][1] == coords[j][1] and coords[i][0] == coords[j][0])
                        and not (abs(i - j) == 1)):
                            edges.append([i,j])
                    else:
                        if (((coords[i][0] + 1 == coords[j][0] or coords[i][0] - 1 == coords[j][0]) and coords[i][1] == coords[j][1]
                        or (coords[i][1] + 1 == coords[j][1] or coords[i][1] - 1 == coords[j][1]) and coords[i][0] == coords[j][0])
                        and not (abs(i - j) == 1)):
                            edges.append([i,j])
                j += 1
        i += 1
    
    return edges

def gen_progress(path, protein, runs, gens, _3D):
    gen_progress = []
    for run in range(runs):
        gen_prog = []
        h_bonds = []
        highest_stability = 0

        folder = f'runs/run{run}'
        folds = []
        with open(f'{path}/{protein}/{folder}/folds.txt', 'r') as f:
            folds = f.readlines()

        last_edges = []
        tmp_max_stb = 0
        for entry in folds:
            stb = int(entry.split("[")[1].split("]")[1].replace(" ", "").replace("\n", ""))
            if stb < tmp_max_stb:
                fold = entry.split("[")[1].split("]")[0].replace(" ", "")
                gen = int(entry.split("[")[0].replace(" ", ""))

                new_edges = []
                if _3D:
                    new_edges = free_energy_edges(protein, get_coords_3D(fold), True)
                else:
                    new_edges = free_energy_edges(protein, get_coords(fold), False)

                if last_edges != new_edges:
                    last_edges = new_edges

                tmp_max_stb = stb

        for entry in folds:
            stb = int(entry.split("]")[1].replace(" ", "").replace("\n", ""))
            if stb < highest_stability:
                fold = entry.split("[")[1].split("]")[0].replace(" ", "")
                gen = int(entry.split("[")[0].replace(" ", ""))

                new_edges = []
                if _3D:
                    new_edges = free_energy_edges(protein, get_coords_3D(fold), True)
                else:
                    new_edges = free_energy_edges(protein, get_coords(fold), False)

                for h_bond in h_bonds:
                    if h_bond not in new_edges:
                        h_bonds.remove(h_bond)

                for edge in new_edges:
                    if (not edge in h_bonds) and (edge in last_edges):
                        h_bonds.append(edge)
                        gen_prog.append([gen, edge, stb])

                highest_stability = stb

        gen_progress.append([highest_stability, gen_prog])
        # del folds
        # gc.collect()
    
    highest_total_stability = 0
    for progress in gen_progress:
        if progress[0] < highest_total_stability:
            highest_total_stability = progress[0]

    energy_bucket = [[-i,[]] for i in range(abs(highest_total_stability - 1))]
    energy = 0
    while energy < len(energy_bucket):
        for progress in gen_progress:
            if progress[0] == -energy:
                energy_bucket[energy][1].append(progress[1])

        energy += 1
    # for bucket in energy_bucket:
    #     print(bucket)

    # del gen_progress
    # gc.collect()

    if not os.path.isdir(f'{path}/{protein}/results/generational_progress/ind'):
        os.makedirs(f'{path}/{protein}/results/generational_progress/ind')
    for bucket in energy_bucket:
        # print(bucket)
        # print("")
        # per energy
        if bucket[1] != [] and bucket[1] != [[]]:
            max_y = 0
            x_labels = []
            y_points = []
            avg_gen_for_max_stb = 0
            # per run
            for run in bucket[1]:
                avg_gen_for_max_stb += run[-1][0]
                if max_y < run[-1][0]:
                    max_y = run[-1][0]
                for bond in run:
                    x_labels.append(bond[1][0])
                    y_points.append(bond[0])
            x, y = (list(p) for p in zip(*sorted(zip([int(i) for i in x_labels], y_points))))
            x = [str(p) for p in x]

            mpl.rcParams["font.size"] = 15
            fig = plt.figure(figsize=(25, 8))
            ax = fig.add_subplot(111)
            ax.scatter(x, y, color="#3e559c", label=str(bucket[0]))
            plt.axhline(y=avg_gen_for_max_stb/len(bucket[1]), color='blue', label="Average generation of highest achieved stability")
            
            plt.legend()
            if _3D:
                plt.title(f'[3D] Generational progress of stability {bucket[0]} (L={len(protein)}, G={gens})', fontsize=28)
            else:
                plt.title(f'[2D] Generational progress of stability {bucket[0]} (L={len(protein)}, G={gens})', fontsize=28)
            plt.xlabel("H-amino index", fontsize=20)
            plt.ylabel("Generation", fontsize=20)
            plt.yscale('log')

            plt.savefig(f'{path}/{protein}/results/generational_progress/ind/' + str(bucket[0]) + "_log.png")
            plt.close()

            # del x_labels
            # del y_points
            # gc.collect()

def gen_progress_two_stabilities(path, protein, runs, gens, stability1, stability2, _3D):
    gen_progress = []
    for run in range(runs):
        gen_prog = []
        h_bonds = []
        highest_stability = 0

        folder = f'runs/run{run}'
        folds = []
        with open(f'{path}/{protein}/{folder}/folds.txt', 'r') as f:
            folds = f.readlines()

        last_edges = []
        tmp_max_stb = 0
        for entry in folds:
            stb = int(entry.split("[")[1].split("]")[1].replace(" ", "").replace("\n", ""))
            if stb < tmp_max_stb:
                fold = entry.split("[")[1].split("]")[0].replace(" ", "")
                gen = int(entry.split("[")[0].replace(" ", ""))

                new_edges = []
                if _3D:
                    new_edges = free_energy_edges(protein, get_coords_3D(fold), True)
                else:
                    new_edges = free_energy_edges(protein, get_coords(fold), False)

                if last_edges != new_edges:
                    last_edges = new_edges

                tmp_max_stb = stb

        for entry in folds:
            stb = int(entry.split("]")[1].replace(" ", "").replace("\n", ""))
            if stb < highest_stability:
                fold = entry.split("[")[1].split("]")[0].replace(" ", "")
                gen = int(entry.split("[")[0].replace(" ", ""))

                new_edges = []
                if _3D:
                    new_edges = free_energy_edges(protein, get_coords_3D(fold), True)
                else:
                    new_edges = free_energy_edges(protein, get_coords(fold), False)

                for h_bond in h_bonds:
                    if h_bond not in new_edges:
                        h_bonds.remove(h_bond)

                for edge in new_edges:
                    if (not edge in h_bonds) and (edge in last_edges):
                        h_bonds.append(edge)
                        gen_prog.append([gen, edge, stb])

                highest_stability = stb

        gen_progress.append([highest_stability, gen_prog])
        del folds
        gc.collect()

    energy_bucket = [[-i,[]] for i in [abs(stability1), abs(stability2)]]
    energy = 0
    while energy < len(energy_bucket):
        for progress in gen_progress:
            if progress[0] == stability1:
                energy_bucket[0][1].append(progress[1])
            if progress[0] == stability2:
                energy_bucket[1][1].append(progress[1])

        energy += 1

    del gen_progress
    gc.collect()

    if not os.path.isdir(f'{path}/{protein}/results/generational_progress/comp'):
        os.makedirs(f'{path}/{protein}/results/generational_progress/comp')
    x_labels = []
    y_points = []
    y_avg = []
    for bucket in energy_bucket:
        # per energy
        max_y = 0
        x_label = []
        y_point = []
        avg_gen_for_max_stb = 0
        # per run
        for run_ in bucket[1]:
            avg_gen_for_max_stb += run_[-1][0]
            if max_y < run_[-1][0]:
                max_y = run_[-1][0]
            for bond in run_:
                x_label.append(str(bond[1][0]))# + "\n-\n" + str(bond[1][1]))
                y_point.append(bond[0])

        x_labels.append(x_label)
        y_points.append(y_point)
        y_avg.append(avg_gen_for_max_stb/len(bucket[1]))

    fig = plt.figure(figsize=(30, 8))
    ax = fig.add_subplot(111)

    plot_x = []
    plot_y = []
    for x_label, y_point in zip(x_labels, y_points):
        x, y = (list(p) for p in zip(*sorted(zip([int(i) for i in x_label], y_point))))
        plot_x.append([str(p) for p in x])
        plot_y.append(y)

    color_wheel = [["#d43338", "red"], ["#3e559c", "blue"]]
    stbs = [stability1, stability2]
    trans = [0.4, 0.2]
    i = 0
    while i < len(x_labels):
        ax.scatter(plot_x[i], plot_y[i], color=color_wheel[i][0], alpha=trans[i], label=f'{stbs[i]}')
        plt.axhline(y=y_avg[i], color=color_wheel[i][1], label=f'Average generation of highest achieved stability of {stbs[i]}')
        i += 1

    plt.legend()
    # plt.title(f'Generational progress of the H-bond forming process for stability {stability1} and {stability2}', fontsize=28)
    if _3D:
        plt.title(f'[3D] Generational progress of stability {stability1} and {stability2} (L={len(protein)}, G={gens})', fontsize=28)
    else:
        plt.title(f'[2D] Generational progress of stability {stability1} and {stability2} (L={len(protein)}, G={gens})', fontsize=28)
    plt.xlabel("H-amino index", fontsize=20)
    plt.ylabel("Generation", fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.yscale('log')

    plt.savefig(f'{path}/{protein}/results/generational_progress/comp/{stability1}~{stability2}_log.png')
    plt.close()

# def setup():
#     print("start: " + str(datetime.datetime.now()))
#     runs = 1000
#     length = 100
#     proteins = ["HPPPPPPPHPPPPPPPPHPHPHHPHPHHHPPPPPPPHHPPPPHPPHHHPPPHHPPHPPPHHHPPPPPHPPPPPPPHHPPHPHHPHPPPHPHPHPHHPHPH", "PPHHHPPHHPHPHHPHHHHHHHPHPPPPHPPPPHHPHHHHHPHHHPPHHPHPPPPHHPHHPHPHHHPPPPHPHHHPHPHPHPHPHHHPPPHHPPPHHHHH", "PPHHPPHPPHHHHPHPHPHHPPPHPHPPPPHHPHHHPPPHHHPPPHHHHPPPPHHHHPPHHPPHPPPHHPPHHPPPPPHHHPHPPPPHPPPHHHPPHHPP"]
#     stbs = [[[-12, -18], [-7, -21]], [[-18, -25], [-11, -30]], [[-16, -22], [-6, -26]]]
#     stb_index = 0
#     for protein in proteins:
#         print("     " + protein)
#         path = f'1000_sample/length_{length}'
#         os.makedirs(f'{path}/{protein}/results/generational_progress/ind')
#         os.makedirs(f'{path}/{protein}/results/generational_progress/comp')
        
#         print("     - individual")
#         gen_progress(path, protein, runs)
#         print(f'     - compare {stbs[stb_index][0][0]} and {stbs[stb_index][0][1]}')
#         gen_progress_two_stabilities(path, protein, runs, stbs[stb_index][0][0], stbs[stb_index][0][1])
#         print(f'     - compare {stbs[stb_index][1][0]} and {stbs[stb_index][1][1]}')
#         gen_progress_two_stabilities(path, protein, runs, stbs[stb_index][1][0], stbs[stb_index][1][1])
        
#         stb_index += 1
#         print(" next: " + str(datetime.datetime.now()))

# def setup():
#     runs = 1000
#     lengths = [100]
#     for length in lengths:
#         path = f'1000_sample/length_{length}'
#         protein = "PPHHPPHPPHHHHPHPHPHHPPPHPHPPPPHHPHHHPPPHHHPPPHHHHPPPPHHHHPPHHPPHPPPHHPPHHPPPPPHHHPHPPPPHPPPHHHPPHHPP"
#         # os.makedirs(f'{path}/{protein}/results/generational_progress/ind')
#         # os.makedirs(f'{path}/{protein}/results/generational_progress/comp')
#         gen_progress_two_stabilities(path, protein, runs, -13, -24)
#         gen_progress_two_stabilities(path, protein, runs, -9, -26)

# setup()