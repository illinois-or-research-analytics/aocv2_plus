-- integer_ids from exosome data that are in the retraction_doi list

SELECT distinct retraction_doi,edcn.integer_id,edcn.year FROM cleaned_rw cr
INNER JOIN dimensions.exosome_dimensions_complete_nodelist edcn
  ON edcn.doi=lower(cr.retraction_doi);

-- integer_ids from exosome data that are in the original_doi list
with cte as(
SELECT edcn.integer_id,edcn.year,
       original_paper_doi, orig_year,
       retraction_doi, retraction_year
FROM cleaned_rw cr
INNER JOIN dimensions.exosome_dimensions_complete_nodelist edcn
  ON edcn.doi=lower(cr.original_paper_doi))
select distinct retraction_year,orig_year,original_paper_doi from cte;

SELECT MAX(record_id),original_paper_doi, COUNT(DISTINCT record_id)
FROM cleaned_rw
WHERE original_paper_doi != 'Unavailable'
OR original_paper_doi != 'unavailable'
GROUP BY original_paper_doi
HAVING COUNT(DISTINCT record_id) > 1;

SELECT retraction_doi,original_paper_doi, max(record_id) OVER (PARTITION BY original_paper_doi) as record
FROM cleaned_rw where original_paper_doi='10.1001/archpediatrics.2012.999'

select count(record_id) from cleaned_rw where retraction_doi is not null;