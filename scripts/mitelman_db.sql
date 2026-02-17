# This is a bigquery script to download fusion events from mitelman database
# Run with google account and update dates accordingly, 4 releases per year
# Check `isb-cgc-bq.mitelman_versioned.INFORMATION_SCHEMA.TABLES`
SELECT g.Gene, COUNT(*) AS CNT
FROM `isb-cgc-bq.mitelman_versioned.MolBiolClinAssoc_2026_01` c
JOIN `isb-cgc-bq.mitelman_versioned.MolClinGene_2026_01` g
  ON (g.RefNo = c.RefNo 
      AND g.InvNo = c.InvNo 
      AND g.MolClin = c.MolClin)
WHERE g.Gene LIKE '%::%'
  AND g.MolClin = 'M'
  AND c.GeneShort NOT LIKE '%::%,%::%'
GROUP BY Gene
ORDER BY CNT DESC, Gene
