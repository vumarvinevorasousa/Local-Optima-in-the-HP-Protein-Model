import os
import numpy as np
import matplotlib.pyplot as plt
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

def plot_similarity_protein(path, protein, _3D):
    f = open(f'{path}/{protein}/results/best/best_folds.txt')
    best_folds = f.read().splitlines()

    edge_matrices = []
    energy_index = -1
    for fold_line in best_folds:
        if fold_line.lstrip("-").isdigit():
            energy_index += 1

            edge_matrix = [[0 for i in range(len(protein))] for j in range(len(protein))]
            edge_matrices.append(edge_matrix)
        else:
            fold = fold_line.split()[0]
            edges = []
            if _3D:
                edges = free_energy_edges(protein, get_coords_3D(fold), True)
            else:
                edges = free_energy_edges(protein, get_coords(fold), False)
                
            for edge in edges:
                edge_matrices[energy_index][edge[0]][edge[1]] += 1

    if not os.path.isdir(f'{path}/{protein}/results/similarity'):
        os.mkdir(f'{path}/{protein}/results/similarity')

    energy = 0
    for matrix in edge_matrices:
        count = 0
        dist = 0
        i = 1
        while i < len(protein):
            j = i + 1
            while j < len(protein):
                if matrix[i][j] > 0:
                    count += 1
                    dist += j-i
                j += 1
            i += 1

        diag = 0
        if count > 0:
            diag = round(dist/count, 2)

        if _3D:
            plt.title(f'[3D] Stability {energy}, D={diag} (L={len(protein)})', fontsize=18)
        else:
            plt.title(f'[2D] Stability {energy}, D={diag} (L={len(protein)})', fontsize=18)
        plt.imshow(matrix, cmap='Blues', interpolation='nearest', origin='lower')
        plt.colorbar()
        plt.xlabel("Amino index")
        plt.ylabel("Amino index")
        plt.savefig(f'{path}/{protein}/results/similarity/{energy}.png')
        plt.close()

        energy -= 1
        