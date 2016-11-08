clusters = [cluster, cluster, cluster]

cluster = [{gene_id:1, phage_id:1, translation:"proteinString", start:0, end:227, locus_tag:"name",blastp:[hits],clustulo:[diffHits]}, {gene}, {gene}]

hits = [{query_id:"same as above gene_id", subject_id:222, query_start:1, subject_start:4, e_value:float, percent_identity(ident):float}, {hit}, {hit}]

for cluster in clusters:
	# Build a dictionary of golden genes
	{1:[gene_ids], 2:[gene_ids]}

	for gene in cluster:
		if gene != golden:
			# find a one to one match within the golden phage dictionary starting at 1
			# one to one is the same start codon
			listOfTuples = [tuple{postAdjustmentDistance, goldenRank, %ident, ...more?}]


	for each goldenHit in gene[blastp]:
		relative_gene_start = hit["query_start"]
		relative_golden_start = hit["subject_start"]

		if gene_start != 1 and golden_start != 1:
			do nothing

		elif:

			gene going forward (absolute_gene_start > absolute_gene_end):
				tooLong(relative_golden_start == 1)
				tooShort(relative_gene_start == 1)

			gene going backward (absolute_gene_end > absolute_gene_start):