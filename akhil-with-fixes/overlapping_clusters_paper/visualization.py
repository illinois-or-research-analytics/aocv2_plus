import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import statistics


'''
Method to analyze a dictionary of information in the form of a scatterplot

Input:
  information {dict} - dictionary of information to analyze
  x_col str - column name to graph on x-axis
  x_scale str - scaling factor for x-axis
  y_col str - column name to graph on y-axis
  y_scale str - scaling factor for y-axis
  hue_val str - hue value for third variable in scatterplot
  output_file_path str - file path to save output scatterplot figure

Output:
  None
'''
def scatterplot_analysis(information, x_col, x_scale, y_col, y_scale, hue_val, output_file_path):
  df = pd.DataFrame.from_dict(information)
  
  plt.figure()
  ax = sns.scatterplot(data=df, x=x_col, y=y_col, hue=hue_val)
  ax.set_xscale(x_scale)
  ax.set_yscale(y_scale)
  plt.savefig(output_file_path)

'''
Method to analyze a dictionary of information in the form of a scatterplot

Input:
  frequency_list {dict} - dictionary of  information to analyze
  output_file_path str - file path to save output scatterplot figure
  min_k_core int - minimum k value to parse
  inclusion criterion str - criteria to include an node into a cluster {k, mcd}

Output:
  None
'''
def histogram_analysis(frequency_list, output_file_path, min_k_core, inclusion_criterion):
  plt.figure()
  print(min(frequency_list), max(frequency_list), statistics.median(frequency_list))
  sns.histplot(frequency_list,  bins=100, log_scale = True) 

  plt.xlabel('Intersection Size')
  plt.ylabel('Frequency')
  plt.savefig(str(min_k_core) + '_' + str(inclusion_criterion) + '_'  + output_file_path)

