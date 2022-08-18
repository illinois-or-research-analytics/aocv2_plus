import pstats

if __name__ == '__main__':
  p = pstats.Stats('complex_nodes_profile.txt')
  p.sort_stats('tottime').print_stats(25)
