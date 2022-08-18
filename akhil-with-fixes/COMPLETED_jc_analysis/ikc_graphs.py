import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import matplotlib.pyplot as plt
import seaborn as sns
import os
from collections import defaultdict
import pandas as pd
import plotly.figure_factory as ff
import matplotlib


def generate_graphs():
  df = pd.read_csv('clusterings.csv')
  for i in range(0, 31, 5):
    df1 = df.loc[(df['k'] >= 5) & (df['cluster_size'] > 1) & (df['year'] == i)]
    fin_df = df1.groupby(['jc'])['cluster_size', 'k'].describe()
    print(fin_df.to_string())

  '''
  df1 = df.loc[(df['k'] >= 5) &  (df['k'] <= 20) & (df['cluster_size'] > 1) & (df['cluster_size'] < 100) & (df['year'] == 30)]
  df2 = df.loc[(df['k'] >=20) & (df['cluster_size'] > 1) & (df['cluster_size'] < 100) & (df['year'] == 30)]
  df3 = df.loc[(df['k'] >= 5) &  (df['k'] <= 20) & (df['cluster_size'] > 100) & (df['year'] == 30)]
  df4 = df.loc[(df['k'] >= 20) & (df['cluster_size'] > 100) & (df['year'] == 30)]
  dfs = [df1, df2, df3, df4]
  #print('Summary Statistics for Exosome Citation Network Peelback Year ', 1990+year)
  
  #print('\n')
  #print('Summary Statistics for k-values of clusters grouped by JC')
  for df in dfs:
    fin_df = df.groupby(['jc'])['cluster_size', 'k'].describe()
    print(fin_df.to_string())
  '''
 
  
  '''
  cmap = sns.cubehelix_palette(rot=-.2, as_cmap=True)
  g = sns.relplot(
      data=df,
      x="k", y="cluster_size",
      hue="year", size="jc",
      palette=cmap, sizes=(10, 200),
  )
  g.set(xscale="linear", yscale="log")
  g.ax.xaxis.grid(True, "minor", linewidth=.25)
  g.ax.yaxis.grid(True, "minor", linewidth=.25)
  g.despine(left=True, bottom=True)
  
  plt.savefig('../plots/jc_year_relplot.png')
  #print('Summary Statistics for cluster size of clusters grouped by JC')
  #print(df_year.groupby(['jc'])['cluster_size'].describe())

  '''
  '''
  plt.figure()
  sns.scatterplot(data=df_year, x='k', y='cluster_size', hue='jc')
  plt.yscale('symlog')
  plt.xlim(0,50)
  plt.savefig('../plots/jc_scatterplot' + str(1990+year) + '.png')
  print(1990+year)
  '''
def main():
  clusters_list = []
  for jc in [100, 250, 500, 1000]:
    for year in range(0, 31, 5):
      filepath = './jc_clusterings/jc'+str(jc)+'year'+str(1990+year)+'ikc.clustering'
      print(filepath)
      f = open(filepath)
      clusters = defaultdict(list)
      for line in f.readlines():
        line = [i.strip() for i in line.split(',')]
        cluster_index = str(line[1]) + '-' + str(line[2])
        clusters[cluster_index].append(line[0])
        
      clusters_list.append((jc, year, clusters))

  to_df = defaultdict(list)
  for item in clusters_list:
    for key, value in item[2].items():
      k = int(key.split('-')[1])
      to_df['k'].append(k)
      to_df['year'].append(item[1])
      to_df['jc'].append(item[0])
      to_df['cluster_size'].append(len(item[2][key]))
  
  df = pd.DataFrame.from_dict(to_df)
  df.to_csv('clusterings.csv')
  '''
  sns.scatterplot(data=df, x='k', y='cluster_size', hue='jc')
  plt.yscale()
  plt.savefig('../plots/test_jc.png')
  '''

if __name__ == '__main__':
  generate_graphs()
