import csv
import pandas as pd
import psycopg2 as pg

# Helper code to generate csv file of marker data from the original set of marker nodes from PostgresQL table

def main():
  conn = pg.connect(
      database='ernieplus',
  )
  cursor = conn.cursor()


  query = 'select doi, integer_id from ernieplus.public.marker_nodes_integer_pub'
  sql_query = pd.read_sql_query(query, conn)
  df = pd.DataFrame(sql_query)
  df.to_csv(r'~/spring_2022_research/akhil/overlapping_clusters_paper/original_markers.csv', index = False)

  if conn:
    cursor.close()
    conn.close()


if __name__ == '__main__':
  main()
