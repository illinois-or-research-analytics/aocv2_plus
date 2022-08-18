import os




def run_ikc(network_file, k=5):
  os.chdir('/home/akhilrj2/ERNIE_Plus/Illinois/clustering/eleanor/code')
  print('pipenv run python3 IKC.py -e ' + '/home/akhilrj2/spring_2022_research/akhil/jc_analysis/jc_tsv/' +  str(network_file) + '.tsv -o ../../../../../spring_2022_research/' + str(network_file) + 'ikc.clustering -k ' + str(k))
  os.system('pipenv run python3 IKC.py -e ' + '/home/akhilrj2/spring_2022_research/akhil/jc_analysis/jc_tsv/'+ str(network_file) + '.tsv -o ../../../../../spring_2022_research/' + str(network_file) + 'ikc.clustering -k ' + str(k))


def main():
  for jc in [100, 250, 500, 1000]:
    for year in range(0, 31, 5):
      run_ikc('jc'+str(jc)+'year'+str(1990+year))



if __name__ == '__main__':
  main()
