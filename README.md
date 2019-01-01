# litSearchTruthsetTesting
There are 2 truth set files (one for BRCA1, one for BRCA2) that contain LOVD literature data. They are the first 2 arguments for the repository's main module, litSearchVsTruthSet. BRCA1 truth set must be first argument (arg O), BRCA2 truth set must be second (arg 1). 

A 3rd file, labeled exampleLitSearchOutput can be used as a 3rd argument (arg 2) for the main module, litSearchVsTruthSet, to test LitSearch output against the truthset. Output: a pubs dictionary (in its own outPubs.py file) and performance stats (printed in terminal). Use this example output file to model any new LitSearch output to be used accurately in litSearchVsTruthSet, which can then be used as a new 3rd argument. 

Use the resulting pubs data structure to interrogate any PMID of your choice. Pubs data structure is provided below: 

pubs = {'PMID': { ('variantalias1', 'variantalias2', ...) : {'inMunchOutput': boolean, 'inTruthSet': boolean} }, 'PMID':{():{}}, ...  }
