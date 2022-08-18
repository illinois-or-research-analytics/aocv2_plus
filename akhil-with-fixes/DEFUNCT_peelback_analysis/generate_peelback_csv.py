import csv
import pandas as pd
import psycopg2 as pg
#from sqlalchemy import create_engine


def main():
  conn = pg.connect(
    database="ernieplus",
  )
  cursor = conn.cursor()

  for i in range(0, 31, 5):
    peelback_query = "select citing_integer_id, cited_integer_id from akhilrj2.peelback_table where peelback_year = " + str(1990 + i)
    sql_query = pd.read_sql_query(peelback_query, conn)
    peelback_csv = "peelback" + str(1990+i) + ".tsv"
    df = pd.DataFrame(sql_query)
    df.to_csv (r'/home/akhilrj2/spring_2022_research/akhil/peelback_tsv/' + peelback_csv, sep='\t', index = False, header=False)

  if conn:
    cursor.close()
    conn.close()


if __name__ == "__main__":
  main()
