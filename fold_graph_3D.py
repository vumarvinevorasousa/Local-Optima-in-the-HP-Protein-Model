import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
import operator
import numpy as np

op = {
    "+": operator.add,
    "-": operator.sub
}

def free_energy_edges(protein, coords):
    edges = []
    i = 0
    while i < len(protein):
        if protein[i] == "H":
            coords_run = []
            j = i + 1
            while j < len(protein):
                if protein[j] == "H":
                    if (((coords[i][0] + 1 == coords[j][0] or coords[i][0] - 1 == coords[j][0]) and coords[i][1] == coords[j][1] and coords[i][2] == coords[j][2]
                        or (coords[i][1] + 1 == coords[j][1] or coords[i][1] - 1 == coords[j][1]) and coords[i][0] == coords[j][0] and coords[i][2] == coords[j][2]
                        or (coords[i][2] + 1 == coords[j][2] or coords[i][2] - 1 == coords[j][2]) and coords[i][1] == coords[j][1] and coords[i][0] == coords[j][0])
                        and not (abs(i - j) == 1)):
                        coords_run.append(coords[i])
                        coords_run.append(coords[j])
                j += 1
            if coords_run:
                edges.append(coords_run)
        i += 1
    
    return edges

def opposite_operator(op):
    if op == "+":
        return "-"
    else:
        return "+"

def check_orientation(current, last):
    if current[0] != last[0]:
        if current[0] - last[0] == 1:
            return "+x"
        else:
            return "-x"
    elif current[1] != last[1]:
        if current[1] - last[1] == 1:
            return "+y"
        else:
            return "-y"
    else:
        if current[2] - last[2] == 1:
            return "+z"
        else:
            return "-z"

def match_orientation(orientation, d, coord):
    coords = []
    oper = orientation[0]
    if oper == "-":
            oper = opposite_operator(oper)
    
    match d:
        case "f":
            coords = [coord[0], op[oper](coord[1], 1), coord[2]]
        case "b":
            coords = [coord[0], op[opposite_operator(oper)](coord[1], 1), coord[2]]
        case "l":
            coords = [op[opposite_operator(oper)](coord[0], 1), coord[1], coord[2]]
        case "r":
            coords = [op[oper](coord[0], 1), coord[1], coord[2]]
    
    return coords

def get_coords_3D(fold):
    coords = [[0,0,0], [0,1,0]]
    orientation = "+y"
    for d in fold[1:]:
        current = coords[-1]
        match d:
            case "u":
                coords.append([current[0], current[1], current[2] + 1])
            case "d":
                coords.append([current[0], current[1], current[2] - 1])
            case _:
                coords.append(match_orientation(orientation, d, current))
                orientation = check_orientation(coords[-1], current)
    
    return coords

def graph_fold_3D(protein, fold, save_at):
    coords = get_coords_3D(fold)

    fig = plt.figure()
    ax = plt.axes(projection='3d')

    color_map = []
    for amino in protein:
        if amino == "H":
            color_map.append("#3b5dc4")
        else:
            color_map.append("#f0a326")

    xdata = []
    ydata = []
    zdata = []
    for i in range(len(coords)):
        xdata.append(coords[i][0])
        ydata.append(coords[i][1])
        zdata.append(coords[i][2])

    ax.scatter3D(xdata, ydata, zdata, c=color_map)

    energy_edges_collection = free_energy_edges(protein, coords)
    xedges = [[] for i in range(len(energy_edges_collection))]
    yedges = [[] for i in range(len(energy_edges_collection))]
    zedges = [[] for i in range(len(energy_edges_collection))]
    for i in range(len(energy_edges_collection)):
        for j in range(len(energy_edges_collection[i])):
            xedges[i].append(energy_edges_collection[i][j][0])
            yedges[i].append(energy_edges_collection[i][j][1])
            zedges[i].append(energy_edges_collection[i][j][2])
    
    ax.plot3D(xdata, ydata, zdata, c='black')
    for i in range(len(xedges)):
        ax.plot3D(xedges[i], yedges[i], zedges[i], linestyle='dashed', c='black')

    ax.set_zticks([])
    ax.set_yticks([])
    ax.set_xticks([])

    plt.savefig(save_at)
    plt.show()
    plt.close()
    