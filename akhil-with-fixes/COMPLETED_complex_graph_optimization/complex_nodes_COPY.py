import time 
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import pandas as pd
"""
Given a simple undirected graph with edges between nodes, transform into a complex graph
and find edges between complex nodes.
"""


def build_complex_graph_dict(simple_graph: dict):
    complex_edges_graph = dict()
    edges_added = set()  # to avoid duplicate additions into the graph

    for n1, n2_set in simple_graph.items():
        for n2 in n2_set:
            n3_set = simple_graph[n2]
            for n3 in n3_set:
                cn1 = tuple(sorted([n1, n2]))  # first complex node
                cn2 = tuple(sorted([n2, n3]))  # second complex node

                if cn1 == cn2:  # check if both are same complex nodes (e.g. (1,8) and (8,1))
                    # node cannot point itself
                    continue

                possible_same_edges = {(cn1, cn2), (cn2, cn1)}
                if len(possible_same_edges & edges_added) > 0:  # check if edge already in the graph
                    continue

                # add complex edge
                edge = cn1, cn2
                if cn1 not in complex_edges_graph:
                    complex_edges_graph[cn1] = set()
                complex_edges_graph[cn1].add(cn2)
                edges_added.add(edge)

    return complex_edges_graph

def build_complex_graph_dict_v2(simple_graph):

  complex_edges_graph = defaultdict(list)

  for n1, n1_neighbors in simple_graph.items():
    for i in range(len(n1_neighbors)):
      for j in range(i+1, len(n1_neighbors)):
        edge = ((n1, n1_neighbors[i]), (n1, n1_neighbors[j]))
        complex_edges_graph[edge[0]].append(edge[1])

  return complex_edges_graph

if __name__ == "__main__":
    num_edges = []
    time_list = []
    typ_list = []
    input_graph_file_txt = '../peelback_code/peelback_tsv/peelback' + str(1990) + '.tsv'
    plt.figure(figsize=(8, 6), dpi=80)
    for _ in range(5):
      for edge_count in range(20):
        for typ in ['regular', 'optimized']:
          num_edges.append(2**edge_count)
          start_time = time.time()
          
          # STEP 1: Read input file, build simple graph, and print edges
          simple_graph = dict()  # input graph in dictionary form
          input_graph_num_edges = 0
          #print('INPUT EDGES')
          with open(input_graph_file_txt) as graph_data:
            lines = graph_data.readlines()
            for edge in range(2**edge_count):
                  line = lines[edge]
                  n1, n2 = line.split('\t', 1)
                  n1, n2 = n1.strip(), n2.strip()

                  # simple graph has edge for both directions; this is useful in generating complex edges later
                  if n1 not in simple_graph:
                      simple_graph[n1] = list()
                  simple_graph[n1].append(n2)
                  if n2 not in simple_graph:
                      simple_graph[n2] = list()
                  simple_graph[n2].append(n1)
                  #print(f'{n1} -> {n2}')
                  input_graph_num_edges += 1
          #print(f'Num Undirected Edges in Input Graph: {input_graph_num_edges}')

          #print('\n')

          # STEP 2: Transform simple graph to complex graph and find edges
          if typ == 'regular':
            complex_edges_graph = build_complex_graph_dict(simple_graph)  # complex graph in dictionary form
          else:
            complex_edges_graph = build_complex_graph_dict_v2(simple_graph)  # complex graph in dictionary form

          # STEP 3: Print complex graph edges
          complex_graph_num_edges = 0
          #print('COMPLEX EDGES')
          for cn1, cn2_set in complex_edges_graph.items():
              for cn2 in cn2_set:
                  #print(f'{cn1} -> {cn2}')
                  complex_graph_num_edges += 1
          time_list.append(time.time()-start_time)
          typ_list.append(typ)
        #print(f'Num Undirected Edges in Complex Graph: {complex_graph_num_edges}')
    
    data = {'size': num_edges, 'time': time_list, 'type': typ_list}
    df = pd.DataFrame.from_dict(data)
    sns.barplot(x='size', y='time', hue='type', data=df, palette="Blues_d")
    plt.xlabel("Edges in Graph (# of Edges)")
    plt.ylabel("Time to Create Complex Edge Graph (s)")
    plt.title("Scalability of Complex Edge Graph Code")
    plt.yscale("log")
    locs, labels = plt.xticks()
    plt.setp(labels, rotation=45)
    plt.savefig('../plots/scalability_complex_graph_code.png')
