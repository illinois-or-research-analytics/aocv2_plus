import pandas as pd
import psycopg2 as pg



def main():

  conn = pg.connect(
      database = 'ernieplus'
  )

  cursor = conn.cursor()

  jc = 250
  query = 'select citing_integer_id, cited_integer_id from akhilrj2.jc_table where citing_outdegree <= ' + str(jc) + ' and cited_outdegree <= ' + str(jc)
  sql_query = pd.read_sql_query(query, conn)
  df = pd.DataFrame(sql_query)
  df.to_csv(r'/home/akhilrj2/spring_2022_research/akhil/overlapping_clusters_pipeline/exosome_1900_2010_sabpq_jc250.tsv', sep='\t', index=False, header=False)

  if conn:
    cursor.close()
    conn.close()


if __name__ == '__main__':
  main()
