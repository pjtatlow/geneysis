clusters = [cluster, cluster, cluster]

cluster = [{gene_id:1, phage_id:1, translation:"proteinString", start:0, end:227, locus_tag:"name",blastp:[hits],clustulo:[diffHits]}, {gene}, {gene}]

hits = [{query_id:"same as above gene_id", subject_id:222, query_start:1, subject_start:4, e_value:float, percent_identity(ident):float}, {hit}, {hit}]

startCodons = ["ATG","GTG","TGT"]

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
        
        idealMoveDistance = abs(relative_gene_start - relative_golden_start)
        phageGenome = getPhageGenome(gene['phage_id'])

		if gene_start != 1 and golden_start != 1:
			do nothing

		elif:

            bestGeneStart = gene['start']

			if gene going forward (absolute_gene_start > absolute_gene_end):
				if tooLong(relative_golden_start == 1):
                    move -> along genome (increasing the start)

                    for i in range(1,ideal_move_distance)
                        currentStart = gene['start'] + (3 * i) # increase our start
                        codon = phageGenome[currentStart:currentStart+3] # codon is going forward
                        if codon in startCodons:
                            bestGeneStart = currentStart

				elif tooShort(relative_gene_start == 1):
                    move <- along genome (decreasing the start)

                    for i in range(1,ideal_move_distance)
                        currentStart = gene['start'] - (3 * i) # decrease our start
                        codon = phageGenome[currentStart:currentStart+3] # codon is going forward
                        if codon in startCodons:
                            bestGeneStart = currentStart

			elif gene going backward (absolute_gene_end > absolute_gene_start):
                if tooLong(relative_golden_start == 1):
                    move <- along genome (decreasing the start)

                    for i in range(1,ideal_move_distance)
                        currentStart = gene['start'] - (3 * i) # decrease our start
                        codon = phageGenome[currentStart:currentStart-3] # codon is going backward

                        codon = codon[::-1] # reverses codon because we're going backwards

                        if codon in startCodons:
                            bestGeneStart = currentStart

                elif tooShort(relative_gene_start == 1):
                    move -> along genome (increasing the start)

                    for i in range(1,ideal_move_distance)
                        currentStart = gene['start'] + (3 * i) #increase our start
                        codon = phageGenome[currentStart:currentStart-3] # codon is going backward

                        codon = codon[::-1] # reverses codon because we're going backwards

                        if codon in startCodons:
                            bestGeneStart = currentStart

            #by the time we get here bestGeneStart should have the best index of the genome for our new start.
            
