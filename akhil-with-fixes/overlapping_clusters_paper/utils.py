import csv

'''
Helper method to convert clustering file to core node file

Input:
  clustering str - clustering file to parse
  output_path str - output path to write core node file to
Output:
  None
'''
def generate_core_node_file(clustering, output_path):
  with open(clustering, 'r') as clustering_reader, open(output_path, 'w') as output_path_writer:
    for line in clustering_reader:
      output_path_writer.write(line.split('')[1] + ' ' + line.split(',')[0] + '\n')

'''
Helper method to convert csv file to tsv file

Input:
  input_csv str - input path to csv file to convert
  output_tsv str - output path to write converted tsv file
Output:
  None
'''
def csv_to_tsv(input_csv, output_tsv):
  with open(input_csv,'r') as csvin, open(output_tsv, 'w') as tsvout:
    csvin = csv.reader(csvin)
    tsvout = csv.writer(tsvout, delimiter='\t')

    for row in csvin:
      tsvout.writerow(row)
