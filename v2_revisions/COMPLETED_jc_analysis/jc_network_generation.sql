DROP TABLE IF EXISTS akhilrj2.jc_table;

CREATE TABLE IF NOT EXISTS akhilrj2.jc_table (
  citing_integer_id int,
  cited_integer_id int,
  citing_year int,
  cited_year int,
  citing_outdegree int,
  cited_outdegree int
);

/*DROP TABLE IF EXISTS akhilrj2.outdegree_table;*/

CREATE TABLE IF NOT EXISTS akhilrj2.outdegree_table (
  citing_integer_id int,
  outdegree int
);


CREATE OR REPLACE FUNCTION createOutDegreeTable()
RETURNS  VOID
LANGUAGE plpgsql
AS $$
BEGIN 
    INSERT INTO akhilrj2.outdegree_table(citing_integer_id, outdegree) SELECT citing_integer_id, COUNT(*) as outdegree from dimensions.exosome_1900_2010_sabpq_deduplicated group by citing_integer_id;
END;
$$;


CREATE OR REPLACE FUNCTION createJCTable()
RETURNS  VOID
LANGUAGE plpgsql
AS $$
BEGIN 
    INSERT INTO akhilrj2.jc_table(citing_integer_id, cited_integer_id, citing_year, cited_year, citing_outdegree, cited_outdegree) select d1.citing_integer_id, d1.cited_integer_id, d1.citing_year, d1.cited_year, d1.outdegree, o2.outdegree from (select d.citing_integer_id, d.cited_integer_id, d.citing_year, d.cited_year, o1.outdegree from dimensions.exosome_1900_2010_sabpq_deduplicated as d inner join akhilrj2.outdegree_table as o1 on o1.citing_integer_id = d.citing_integer_id) as d1 inner join akhilrj2.outdegree_table as o2 on o2.citing_integer_id = d1.cited_integer_id;
END;
$$;

SELECT createJCTable();

