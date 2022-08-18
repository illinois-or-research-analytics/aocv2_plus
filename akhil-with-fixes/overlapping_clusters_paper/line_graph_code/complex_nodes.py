"""
Given a simple undirected graph with edges between nodes, transform into a complex graph
and find edges between complex nodes.
"""
from collections import defaultdict


def build_complex_graph_dict(simple_graph: dict):
    complex_edges_graph = dict()
    edges_added = set()  # to avoid duplicate additions into the graph

    for n1, n2_list in simple_graph.items():
        for n2 in n2_list:
            n3_list = simple_graph[n2]
            for n3 in n3_list:
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

    input_graph_file_txt = "exosome_k53_jc250.tsv"
    output_graph_file_txt = 'exosome_k53_jc250_line_graph.tsv'

    # STEP 1: Read input file, build simple graph, and print edges
    simple_graph = dict()  # input graph in dictionary form
    input_graph_num_edges = 0
    print('INPUT EDGES')
    with open(input_graph_file_txt) as graph_data:
        lines = graph_data.readlines()
        for line in lines:
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
    print(f'Num Undirected Edges in Input Graph: {input_graph_num_edges}')

    print('\n')

    # STEP 2: Transform simple graph to complex graph and find edges
    complex_edges_graph = build_complex_graph_dict_v2(simple_graph)  # complex graph in dictionary form

    # STEP 3: Print complex graph edges
    complex_graph_num_edges = 0
    #print('COMPLEX EDGES')
    
    output_writer = open(output_graph_file_txt, 'w')
    
    for cn1, cn2_set in complex_edges_graph.items():
        for cn2 in cn2_set:
            #print(f'{cn1} -> {cn2}')
            output_writer.write(f'{cn1}\t{cn2}\n')
            complex_graph_num_edges += 1

    print(f'Num Undirected Edges in Complex Graph: {complex_graph_num_edges}')
