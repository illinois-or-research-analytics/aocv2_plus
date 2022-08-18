/*DROP TABLE IF EXISTS akhilrj2.peelback_table;*/

CREATE TABLE IF NOT EXISTS akhilrj2.peelback_table (
  citing_integer_id int, 
  cited_integer_id int,
  peelback_year int
);

CREATE OR REPLACE FUNCTION createPeelbackNetwork()
RETURNS VOID 
LANGUAGE plpgsql
AS $$
BEGIN 
 FOR i IN 0..31 LOOP
  INSERT INTO akhilrj2.peelback_table(citing_integer_id, cited_integer_id, peelback_year) SELECT citing_integer_id,cited_integer_id, 1990+i from dimensions.exosome_1900_2010_sabpq_deduplicated where citing_year <= 1990+i and cited_year <= 1990+i;
 END LOOP;
END;
$$;

SELECT createPeelbackNetwork();

CREATE OR REPLACE FUNCTION testPeelbackNetwork()
RETURNS VOID 
LANGUAGE plpgsql
AS $$
DECLARE
    publication_count INT := 0; 
BEGIN 
 FOR i IN 0..31 LOOP
  SELECT COUNT(*) FROM akhilrj2.peelback_table WHERE peelback_year = 1990 + i INTO publication_count;
  raise notice 'year: %, publication_count: %', 1990 + i, publication_count;
 END LOOP;
END;
$$;

SELECT testPeelbackNetwork();

/*SELECT * FROM akhilrj2.peelback_table where peelback_year = 1991*/
