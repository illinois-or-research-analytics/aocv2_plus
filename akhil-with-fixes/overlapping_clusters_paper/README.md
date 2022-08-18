# Overlapping Clustering Paper Contents

## Overlapping KMP Pipeline

This program runs the step 5 of the OKMP method and saves it to your local directory. To run this program, 
run the following code. 

```console
python3 overlapping_kmp_pipeline.py --network-file [path to tsv network file] --clustering [path to input kmp-valid disjoint clustering csv file]
-- output-path [path to save outputted overlapping clustering] --min-k-core [min k of input clustering] --rank-type [percent, percentile] --rank-val [top n percent of candidate nodes to consider adding, opposite if percentile]
-- inclusion-criterion [k, mcd] --candidate-criterion [total_degree, indegree, random, seed] --candidate-file [file path to custom candidate list] --experiment-name [name of experiment being run] --experiment-num [experiment number being run] --config [bool, bool, bool, bool] (run overlapping?, display cluster stats?, include marker node analysis?, save outputs?)
```

Shown below are the command flags and what purpose each flag accomplishes.

- network-file  - str  - file path to network tsv file of the format node\_id[tab]node\_id
- clustering-file - str - file path to clustering tsv file of the format cluster\_id[tab]node\_id
- min-k-core  - int - integer defining the min-k-core by which to add candidate nodes if selecting k as inclusion criterion. It is still mandatory even if not using k, so a dummy value must be filled in
- rank-type  - [percent, percentile] - choose whether to use percent or percentile when calculating list of candidate nodes to generate
- rank-val  - float - the percent or percentile value used to generate the list of candidate nodes
- inclusion-criterion - [k, mcd] - utilize k or mcd of cluster as the inclusion criteria for adding a candidate node to a cluster
- candidate-criterion  - [total\_degree, indegree, random, seed] - criterion by which the rank type and rank val sort the list of nodes in the network in deciding the candidate nodes
- candidate-file  - str - only unrequired field that, if included, will override the other candidate selection values and simply use the nodes specified in the candidate file as candidate nodes. They must be in the format node\_id[new line]node\_id
- experiment-name - str - string to specify the name of the experiment in order to name the statistic csv files that are ouputted with the final clustering file
- experiment-num - int - must correspond to a directory in the environment where the program will be run with the format experiment\_{experiment-num}
- config - [bool, bool, bool, bool] - list of four boolean values that correspond to the following criteria. 1. run the overlapping step of the kmp pipeline. 2. Display the cluster statistics to the console. 3. Include an analysis of marker node coverage. 4. Save output clustering files and cluster statistics


## Visualization

Utility file to run histogram analysis or scatterplot analysis on clusterings


## Overlapping Clustering Stats

Utility file ot run other analysis on clustering files outputted from the Overlapping KMP Pipeline


## Utils

Various utility functions that may be needed in preprocessing data


