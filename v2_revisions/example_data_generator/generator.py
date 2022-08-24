import networkx as nx

G = nx.ring_of_cliques(5, 10)
cliques = nx.find_cliques(G)
for i, c in enumerate(cliques):
    for n in c:
        print(i, n)
print("=====")

nx.write_edgelist(G, "ring_of_cliques.txt", data=False)
