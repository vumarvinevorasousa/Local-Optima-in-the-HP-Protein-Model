import operator
import os
import random
from fold_graph_3D import graph_fold_3D
import datetime

op = {
    "+": operator.add,
    "-": operator.sub
}

_2D = "fblr"
_3D = "ud"
 
dir_dict_2D = {
    "f": ["r", "b", "l"],
    "b": ["l", "f", "r"],
    "l": ["f", "r", "b"],
    "r": ["b", "l", "f"]
}

dir_dict_3D_z1 = {
    "u": ["r", "d", "l"],
    "d": ["l", "u", "r"],
    "l": ["u", "r", "d"],
    "r": ["d", "l", "u"]
}

dir_dict_3D_z2 = {
    "u": ["f", "d", "b"],
    "d": ["b", "u", "f"],
    "f": ["d", "b", "u"],
    "b": ["u", "f", "d"]
}

opposite_dict_2D = {
    "f": "b",
    "b": "f",
    "l": "r",
    "r": "l"
}

def calc_free_energy_with_coords_3D(protein, coords):
    free_energy = 0
    i = 0
    while i < len(protein):
        if protein[i] == "H":
            j = i + 1
            while j < len(protein):
                if protein[j] == "H":
                    if (((coords[i][0] + 1 == coords[j][0] or coords[i][0] - 1 == coords[j][0]) and coords[i][1] == coords[j][1] and coords[i][2] == coords[j][2]
                        or (coords[i][1] + 1 == coords[j][1] or coords[i][1] - 1 == coords[j][1]) and coords[i][0] == coords[j][0] and coords[i][2] == coords[j][2]
                        or (coords[i][2] + 1 == coords[j][2] or coords[i][2] - 1 == coords[j][2]) and coords[i][1] == coords[j][1] and coords[i][0] == coords[j][0])
                        and not (abs(i - j) == 1)):
                        free_energy -= 1
                j += 1
        i += 1
    
    return free_energy

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

def get_coords_with_start_3D(fold, coord, orientation):
    coords = [coord]
    orien = orientation
    for d in fold:
        current = coords[-1]
        match d:
            case "u":
                coords.append([current[0], current[1], current[2] + 1])
            case "d":
                coords.append([current[0], current[1], current[2] - 1])
            case _:
                coords.append(match_orientation(orien, d, current))
                orien = check_orientation(coords[-1], current)
    
    return coords

def check_collision_3D(d_pre, d_post):
    coll = False
    if len(d_pre) > len(d_post):
        for d in d_post:
            if d in d_pre:
                coll = True
                break
    else:
        for d in d_pre:
            if d in d_post:
                coll = True
                break

    return coll

def rotate_3D(fold, new_dir):
    new_fold = [new_dir]

    if (new_dir in "lr" and fold[0] in _3D) or (new_dir in _3D and fold[0] in "lr"):
        new_fold = [new_dir]
        options = dir_dict_3D_z1.get(fold[0])
        rotate_index = -1
        for option in options:
            rotate_index += 1
            if option == new_dir:
                break

        for dir in fold[1:]:
            if dir in "fb":
                new_fold.append(dir)
            else:
                new_fold.append(dir_dict_3D_z1[dir][rotate_index])

    elif new_dir in _2D and fold[0] in _2D:
        new_fold = [new_dir]
        options = dir_dict_2D.get(fold[0])
        rotate_index = -1
        for option in options:
            rotate_index += 1
            if option == new_dir:
                break

        for dir in fold[1:]:
            if dir in "ud":
                new_fold.append(dir)
            else:
                new_fold.append(dir_dict_2D[dir][rotate_index])

    elif (new_dir in "fb" and fold[0] in _3D) or (new_dir in _3D and fold[0] in "fb") or (new_dir in _3D and fold[0] in _3D):
        new_fold = [new_dir]
        options = dir_dict_3D_z2.get(fold[0])
        rotate_index = -1
        for option in options:
            rotate_index += 1
            if option == new_dir:
                break

        for dir in fold[1:]:
            if dir in "lr":
                new_fold.append(dir)
            else:
                new_fold.append(dir_dict_3D_z2[dir][rotate_index])
    
    return new_fold

def hill_climb_3D(protein, save_at_folder, runs, gens):
    os.makedirs(save_at_folder)
    for run in range(runs):
        muts_and_colls = []

        # Start with a straight line (all forward)
        fold_forward = ["r" for i in range(len(protein) - 1)]
        coords = [[0,x,0] for x in range(len(protein) - 1)]
        fold_length = len(fold_forward)

        folds = [[0, fold_forward, 0]]
        for gen in range(gens):
            prev_fold = folds[-1]

            # Randomly choose an index in the list of directions to mutate
            index = random.choice(range(1, fold_length))

            # Randomly choose a new direction (only of the available ones)
            dirs = ["f", "b", "l", "r", "u", "d"]
            old_dir = prev_fold[1][index]
            dirs.remove(old_dir)
            if prev_fold[1][index - 1] in _2D:
                if old_dir != opposite_dict_2D[prev_fold[1][index - 1]]:
                    dirs.remove(opposite_dict_2D[prev_fold[1][index - 1]])
            new_dir = random.choice(dirs)

            # Rotate the remainder of the fold past the index of the mutation to find the appropiate result after mutating
            fold_pre = prev_fold[1][:index]
            fold_post = rotate_3D(prev_fold[1][index:], new_dir)
            coords_post = get_coords_with_start_3D(fold_post, coords[index], check_orientation(coords[index], coords[index - 1]))

            # Check for collision and whether the mutated fold has a lower free energy value compared to the previous generation
            if not check_collision_3D(coords[:index], coords_post):
                coords_temp = coords[:index] + coords_post
                new_energy = calc_free_energy_with_coords_3D(protein, coords_temp)
                if new_energy <= prev_fold[2]:
                    muts_and_colls.append("10")

                    coords = coords_temp
                    folds.append([gen, fold_pre + fold_post, new_energy])
                else:
                    muts_and_colls.append("00")
            else:
                muts_and_colls.append("01")

        # For all runs and each generation, we store the generation, the fold and whether a mutation or collision has taken place
        folder = save_at_folder + f'/run{run}'
        os.mkdir(folder)
        with open(folder + "/folds.txt", "w") as f:
            for step in folds:
                f.write(str(step[0]) + " [" + " ".join(step[1]) + "] " + str(step[2]) + "\n")

        with open(folder + "/muts_and_colls.txt", "w") as f:
            for change in muts_and_colls:
                f.write(change + "\n")
