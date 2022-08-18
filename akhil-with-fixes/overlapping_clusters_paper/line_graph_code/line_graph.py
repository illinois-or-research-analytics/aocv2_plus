import networkx as nx
import pickle
from collections import defaultdict
import pickle

# From peelback_code/overlapping_clsuter_sets.py

def generate_graph(network_file):
  return nx.read_edgelist(network_file, delimiter='\t', create_using=nx.DiGraph)


def convert_line_graph_readable(network_file):
  output_file = 'exosome_k53_jc250_line_graph_readable.tsv'
  pickle_file = open('ids_no_a1.pickle', 'wb')
  output_file_writer = open(output_file, 'w')
  ids = defaultdict()
  counter = 0
  with open(network_file, 'r') as network_file_reader:
      line = network_file_reader.readline()
      i = 0
      while len(line) > 0:
        if i % int(1371745288/100) == 0:
          print(int(100*i/1371745288))
        v1, v2 = line.split('\t')
        v1, v2 = v1.strip(), v2.strip()
        if v1 in ids:
          new_v1 = ids[v1]
        else:
          ids[v1] = counter
          new_v1 = ids[v1]
          counter += 1

        if v2 in ids:
          new_v2 = ids[v2]
        else:
          ids[v2] = counter
          new_v2 = ids[v2]
          counter += 1
      
        output_file_writer.write(str(new_v1) + '\t' + str(new_v2) + '\n')
        line = network_file_reader.readline()
        i += 1


  pickle.dump(ids, pickle_file)

def main():
  convert_line_graph_readable('./exosome_k53_jc250_line_graph.tsv')
  #G = nx.read_edgelist('exosome_k53_a1_jc250_line_graph_readable.tsv', delimiter = '\t')
  #print(G.number_of_nodes())
  #print(G.number_of_edges())




if __name__ == "__main__":
  main()
