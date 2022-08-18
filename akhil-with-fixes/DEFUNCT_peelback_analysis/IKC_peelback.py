import os

IKC_directory =  '/home/akhilrj2/ERNIE_Plus/Illinois/clustering/eleanor/code'


def main():
  os.chdir(IKC_directory)
  tsv_directory = '/home/akhilrj2/spring_2022_research/akhil/peelback_tsv/peelback'
  clustering_directory = '/home/akhilrj2/spring_2022_research/akhil/peelback_clusters/peelback'
 
  for i in range(0, 31, 5):
    os.system('pipenv run python3 IKC.py -e ' + tsv_directory + str(1990+i) + '.tsv -o ' + clustering_directory + str(1990+i) + '.clustering -k 56')



if __name__ == '__main__':
  main()
