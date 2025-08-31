import matplotlib.pyplot as plt
import matplotlib as mpl

def plot_average_accepted_mutations_and_stability_from_file(folder, protein, runs, gens, _3D):
    mutations = [[] for i in range(gens)]
    stabilities = [[] for i in range(gens)]

    for i in range(runs):
        f1 = open(f'{folder}/{protein}/runs/run{i}/muts_and_colls.txt')
        f2 = open(f'{folder}/{protein}/runs/run{i}/folds.txt')
        muts = f1.read().splitlines()
        stbs = f2.read().splitlines()

        mut_index = 0
        stb_index = 0
        while mut_index < len(muts):
            mutations[mut_index] += muts[mut_index][0]
            if muts[mut_index][0] == '1':
                stb_index += 1
            stabilities[mut_index].append(stbs[stb_index].split("]")[1].replace(" ", ""))

            mut_index += 1

    avg_muts = []
    avg_stbs = []
    for muts, stbs in zip(mutations, stabilities):
        if not muts == []:
            avg_muts.append(sum(map(int, muts))/len(muts))
            avg_stbs.append(sum(map(int, stbs))/len(stbs))

    mpl.rcParams["font.size"] = 13
    fig, ax1 = plt.subplots(figsize=(10,5))
    ax2 = ax1.twinx()

    ax1.plot(avg_muts, label="mutations")
    ax2.plot(avg_stbs, label="stability", color='red')
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1+h2, l1+l2)

    if _3D:
        fig.suptitle(f'[3D] Accepted Mutations (L={len(protein)}, R={runs}, G={gens})', fontsize=18)
    else:
        fig.suptitle(f'[2D] Accepted Mutations (L={len(protein)}, R={runs}, G={gens})', fontsize=18)
    ax1.set_xlabel("Generation", fontsize=14)
    ax1.set_ylabel("Average of accepted mutations", fontsize=14)
    ax2.set_ylabel("Average of stabilities", fontsize=14)
    ax1.set_xscale('log')

    plt.savefig(f'{folder}/{protein}/results/avg_muts_stbs.png')
    plt.close()
    