import networkx as nx
from matplotlib import pyplot as plt

def free_energy_edges(protein, coords):
    edges = []
    i = 0
    while i < len(protein):
        if protein[i] == "H":
            j = i + 1
            while j < len(protein):
                if protein[j] == "H":
                    if (((coords[i][0] + 1 == coords[j][0] or coords[i][0] - 1 == coords[j][0]) and coords[i][1] == coords[j][1]
                        or (coords[i][1] + 1 == coords[j][1] or coords[i][1] - 1 == coords[j][1]) and coords[i][0] == coords[j][0])
                        and not (abs(i - j) == 1)):
                        edges.append([i,j])
                j += 1
        i += 1
    
    return edges

def calc_figsize(fold):
    coords = [[0,0]]
    for d in fold:
        match d:
            case "u":
                coords.append([coords[-1][0], coords[-1][1] + 1])
            case "d":
                coords.append([coords[-1][0], coords[-1][1] - 1])
            case "l":
                coords.append([coords[-1][0] - 1, coords[-1][1]])
            case "r":
                coords.append([coords[-1][0] + 1, coords[-1][1]])

    lim_x = [0,0]
    lim_y = [0,0]
    for coord in coords:
        lim_x = [min(lim_x[0], coord[0]), max(lim_x[1], coord[0])]
        lim_y = [min(lim_y[0], coord[1]), max(lim_y[1], coord[1])]

    return [len(range(lim_x[0], lim_x[1])) + 1, len(range(lim_y[0], lim_y[1])) + 1]

def graph_fold(protein, fold, save_at):
    # Creates coordinates for each amino in a fold
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

    # Add the aminos in the fold as a node in the graph with their corresponding coordinates and color
    G = nx.Graph()
    for i in range(len(protein)):
        color = ""
        if protein[i] == "H":
            color = "#3b5dc4"
        else:
            color = "#f0a326"
        G.add_node(i, col=color, pos=(pos[i][0], pos[i][1]), edge=color)

    # Add the connections in the fold as edges in the graph
    for i in range(len(fold)):
        G.add_edge(i, i + 1)

    # Add the hydrophilic bonds in the fold as edges in the graph
    h_bond_edges = free_energy_edges(protein, pos)
    for edge in h_bond_edges:
        G.add_edge(edge[0], edge[1])

    edge_color = []
    style_edge = []
    for edge in G.edges():
        if edge[1] - edge[0] == 1:
            edge_color.append("black")
            style_edge.append("-")
        else:
            edge_color.append("black")
            style_edge.append("--")

    pos = nx.get_node_attributes(G, 'pos')
    color_nodes = [n[1] for n in G.nodes(data="col")]
    node_edge = [n[1] for n in G.nodes(data="edge")]

    lims = calc_figsize(fold)
    plt.figure(figsize=(lims[0], lims[1]), dpi=150)
    nx.draw(G, pos=pos, node_color=color_nodes, edgecolors=node_edge, linewidths=3.5, node_size=450)
    nx.draw_networkx_edges(G, style=style_edge, edge_color=edge_color, width=3.0, pos=pos)
    plt.savefig(save_at)
    plt.close()
    