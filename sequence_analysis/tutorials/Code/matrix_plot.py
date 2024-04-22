# Plot transition matrices
# Use with conda environment bnp-mrf

import matplotlib.pyplot as plt
import networkx, os
import numpy as np

import os
import numpy as np

print(os.getcwd())

virtual = True
nb_states = 5

if not(virtual):
    home_dir = "D:\\"
else:
    home_dir = "/mnt/transfert"
base_path = os.path.join(home_dir, "devlp_shared", "AppleCultivarsTreatments")
data_path = base_path + os.sep + "Data"
os.chdir(base_path)


matrix_file = base_path + os.sep + "Models" + os.sep + "seq1v_" + \
    str(nb_states) + "s_L.mat"

TM = np.fromfile(matrix_file)
TM = TM.reshape(nb_states, nb_states)

G = networkx.from_numpy_matrix(TM, create_using=networkx.DiGraph)
blues = cm = plt.get_cmap('Blues')

f = plt.figure(1, figsize=(9, 9))
ax = f.add_subplot(1,1,1)

#networkx.draw_circular(G, node_size=150, node_color='#1f78b4',
                       #with_labels=True, arrows=True, 
                       #arrow_size=140,
                       #width=4,
                       #edge_cmap = plt.cm.Blues)

# G = networkx.star_graph(81)
initial_nodes = np.array(list(G.nodes()))
final_nodes = initial_nodes + 1

networkx.relabel_nodes(G, dict(zip(initial_nodes, final_nodes)), copy=False)

colors = [G.get_edge_data(a, b)['weight'] for (a,b) in G.edges()]    

# networkx.draw_circular(G, edge_color = range(81), edge_cmap = plt.cm.Blues)

options = {
    "node_color": "#A0CBE2",
    "edge_color": colors,
    "width": 3,
    "edge_cmap": plt.cm.Blues,
    "with_labels": True,
    "node_size": 150, 
    "with_labels": True,
}

                       
networkx.draw_circular(G, **options) 
plt.show() 
