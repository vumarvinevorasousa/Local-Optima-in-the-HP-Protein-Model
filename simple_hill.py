import random
import os

rotate_dict = {
    "u": ["r", "d", "l"],
    "d": ["l", "u", "r"],
    "l": ["u", "r", "d"],
    "r": ["d", "l", "u"]
}

def get_coords_with_start(ds, start):
    coords = [start]
    for dir in ds:
        match dir:
            case "u":
                coords.append([coords[-1][0], coords[-1][1] + 1])
            case "d":
                coords.append([coords[-1][0], coords[-1][1] - 1])
            case "l":
                coords.append([coords[-1][0] - 1, coords[-1][1]])
            case "r":
                coords.append([coords[-1][0] + 1, coords[-1][1]])
    
    return coords

def check_collision(d_pre, d_post):
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

def rotate(fold, new_dir):
    new_fold = [new_dir]
    
    options = rotate_dict.get(fold[0])
    rotate_index = -1
    for option in options:
        rotate_index += 1
        if option == new_dir:
            break

    for dir in fold[1:]:
        new_fold.append(rotate_dict[dir][rotate_index])

    return new_fold

def calc_free_energy_with_coords(protein, coords):
    free_energy = 0
    i = 0
    while i < len(protein):
        if protein[i] == "H":
            j = i + 1
            while j < len(protein):
                if protein[j] == "H":
                    if (((coords[i][0] + 1 == coords[j][0] or coords[i][0] - 1 == coords[j][0]) and coords[i][1] == coords[j][1]
                        or (coords[i][1] + 1 == coords[j][1] or coords[i][1] - 1 == coords[j][1]) and coords[i][0] == coords[j][0])
                        and not (abs(i - j) == 1)):
                        free_energy -= 1
                j += 1
        i += 1
    
    return free_energy

def hill_climb(protein, save_at_folder, runs, gens):
    os.makedirs(save_at_folder)
    for run in range(runs):
        muts_and_colls = []

        # Start with a straight line (all to the right)
        fold_right = ["r" for i in range(len(protein) - 1)]
        coords = [[x,0] for x in range(len(protein) - 1)]
        fold_length = len(fold_right)
        
        folds = [[0, fold_right, 0]]
        for gen in range(gens):
            prev_fold = folds[-1]
            
            # Randomly choose an index in the list of directions to mutate
            index = random.choice(range(1, fold_length))
            
            # Randomly choose a new direction (other than the one it already is)
            dirs = ["u", "d", "l", "r"]
            old_dir = prev_fold[1][index]
            dirs.remove(old_dir)
            new_dir = random.choice(dirs)
            
            # Rotate the remainder of the fold past the index of the mutation to find the appropiate result after mutating
            fold_pre = prev_fold[1][:index]
            fold_post = rotate(prev_fold[1][index:], new_dir)
            coords_post = get_coords_with_start(fold_post, coords[index])

            # Check for collisions and whether the mutated fold has a lower free energy value than the previous generation
            if not check_collision(coords[:index], coords_post):
                coords_temp = coords[:index] + coords_post
                new_energy = calc_free_energy_with_coords(protein, coords_temp)
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
                