from functions import *
import argparse


def tooLongForward(db, gene, ideal_move_distance, start_codons, stop_codons):
    # Init bestGeneStart
    farBestGeneStart = None
    closeBestGeneStart = None
    # Run through all the potential starts
    for i in xrange(1,ideal_move_distance*2): # doubled to we have equal search space on both sides
        currentStart = gene['start'] + (3 * i) # increase our start 3 at a time
        phage = get_phage(db, gene['phage_id'])
        phageGenome = phage['seq']
        codon = phageGenome[currentStart:currentStart+3] # codon is going forward

        print codon

        if codon in stop_codons:
            print "Found stop codon at {}".format(currentStart)
            break
        elif codon in start_codons and i > ideal_move_distance:
            print "***far"
            farBestGeneStart = currentStart
            break
        elif codon in start_codons and i <= ideal_move_distance:
            print "***on or before"
            closeBestGeneStart = currentStart
    return closeBestGeneStart, farBestGeneStart


def tooLongBackward(db, gene, ideal_move_distance, start_codons, stop_codons):
    # Init bestGeneStart
    farBestGeneStart = None
    closeBestGeneStart = None
    # Run through all the potential starts
    for i in xrange(1,ideal_move_distance*2): # doubled to we have equal search space on both sides
        currentStart = gene['start'] - (3 * i) # decrease our start 3 at a time
        phage = get_phage(db, gene['phage_id'])
        phageGenome = phage['seq']
        codon = phageGenome[currentStart:currentStart+3] # codon is going backward
        codon = codon[::-1] # reverse the codon

        print codon

        if codon in stop_codons:
            print "Found stop codon at {}".format(currentStart)
            break
        elif codon in start_codons and i > ideal_move_distance:
            print "far"
            farBestGeneStart = currentStart
            break
        elif codon in start_codons and i <= ideal_move_distance:
            print "on or before"
            closeBestGeneStart = currentStart
    return closeBestGeneStart, farBestGeneStart


def tooLong(db, gene, ideal_move_distance, start_codons, stop_codons):
    if gene['start'] < gene['end']:
        print "Forward"
        return tooLongForward(db, gene, ideal_move_distance, start_codons, stop_codons)
    elif gene['start'] > gene['end']:
        print "Backward"
        return tooLongBackward(db, gene, ideal_move_distance, start_codons, stop_codons)



def tooShortForward(db, gene, ideal_move_distance, start_codons, stop_codons):
    # Init bestGeneStart
    farBestGeneStart = None
    closeBestGeneStart = None
    # Run through all the potential starts
    for i in xrange(1,ideal_move_distance*2): # doubled to we have equal search space on both sides
        currentStart = gene['start'] - (3 * i) # decrease our start 3 at a time
        phage = get_phage(db, gene['phage_id'])
        phageGenome = phage['seq']
        codon = phageGenome[currentStart:currentStart+3]

        print codon

        if codon in stop_codons:
            print "Found stop codon at {}".format(currentStart)
            break
        elif codon in start_codons and i > ideal_move_distance:
            print "far"
            farBestGeneStart = currentStart
            break
        elif codon in start_codons and i <= ideal_move_distance:
            print "on or before"
            closeBestGeneStart = currentStart
    return closeBestGeneStart, farBestGeneStart


def tooShortBackward(db, gene, ideal_move_distance, start_codons, stop_codons):
    # Init bestGeneStart
    farBestGeneStart = None
    closeBestGeneStart = None
    # Run through all the potential starts
    for i in xrange(1,ideal_move_distance*2): # doubled to we have equal search space on both sides
        currentStart = gene['start'] + (3 * i) # increase our start 3 at a time
        phage = get_phage(db, gene['phage_id'])
        phageGenome = phage['seq']
        codon = phageGenome[currentStart:currentStart+3]
        codon = codon[::-1] # reverse the codon

        print codon

        if codon in stop_codons:
            print "Found stop codon at {}".format(currentStart)
            break
        elif codon in start_codons and i > ideal_move_distance:
            print "far"
            farBestGeneStart = currentStart
            break
        elif codon in start_codons and i <= ideal_move_distance:
            print "on or before"
            closeBestGeneStart = currentStart
    return closeBestGeneStart, farBestGeneStart


def tooShort(db, gene, ideal_move_distance, start_codons, stop_codons):
    if gene['start'] < gene['end']:
        print "Forward"
        return tooShortForward(db, gene, ideal_move_distance, start_codons, stop_codons)
    elif gene['start'] > gene['end']:
        print "Backward"
        return tooShortBackward(db, gene, ideal_move_distance, start_codons, stop_codons)



# makes the best possible adjustments for a given cluster, aligning any genes that do not belong to 
def adjust_cluster(db,cluster,golden_phage_id,start_codons,stop_codons):
    #first we need to make a list of all the golden phage proteins that are in this create_cluster
    golden_phages = get_golden_phages(db)
    golden_genes = get_golden_genes(golden_phages, cluster)

    print golden_genes.values()
    print golden_phages

    #print golden_genes
    
    if len(golden_genes) == 0: #make sure there is at least one gene from the golden phage in the cluster
        return
    
    for gene in cluster:
        if gene['phage_id'] not in golden_phages:
            closeStart = [] # Because if we have muliple golds we will have many different starts to try?
            farStart = [] # Because if we have muliple golds we will have many different starts to try?
            print
            print "New Gene"            
            possible_adjustments = []
            for index, gold_ids in enumerate(golden_genes.values()):
                print index, gold_ids
                for gold_id in gold_ids:
                    blastp_hit = None
                    for hit in gene['hits']: #find hit in gene for that gold_id
                        if hit['subject_id'] == gold_id:
                            blastp_hit = hit
                            break
                    if blastp_hit is not None: # if we found a hit for that gene, continue
                        #print "Gene", gene['id'], "and gene", gold_id, "has hit", blastp_hit
                        golden_start = blastp_hit['subject_start']
                        gene_start = blastp_hit['query_start']

                        # our gene is too short and we need to move the start upstream
                        if gene_start == 1 and golden_start == 1:
                            print "They are already perfectly aligned!"
                        elif gene_start == 1:
                            print "Too Short"
                            ideal_move_distance = golden_start - gene_start
                            newCloseStart, newFarStart = tooShort(db, gene, ideal_move_distance, start_codons, stop_codons)
                            closeStart.append(newCloseStart)
                            farStart.append(newFarStart)
                        # our gene is too long and we need to trim it down
                        elif golden_start == 1:
                            print "Too Long"
                            ideal_move_distance = gene_start - golden_start
                            newCloseStart, newFarStart = tooLong(db, gene, ideal_move_distance, start_codons, stop_codons)

                            ######## Just to debug #####################
                            print "Starts:", newCloseStart, newFarStart
                            print "Original:", gene['start']
                            if closeStart == None:
                                print "closeStart is None"
                            if farStart == None:
                                print "closeStart is None"
                            ###################################


                            closeStart.append(newCloseStart)
                            farStart.append(newFarStart)
                        # right now we do nothing...
                        else:
                            print "Neither one starts at 1..."


                    else:
                        print "Gene", gene['id'], "has no blastp hit for golden gene", gold_id, gene['hits']







        #     pass
        # if not gene['adjusted'] and gene['phage_id'] != golden_phage_id: # only adjust non-golden genes that have not been adjusted
        #     golden_hit = None
        #     gold_id = None
        #     # find the closest golden gene blastp his
        #     for hit in gene['hits']:
                
        #         # make sure we're dealing with a blastp hit is to a golden gene
        #         if hit['type'] == "blastp" and hit['subject_id'] in golden_genes:
        #             # first golden gene blastp hit
        #             if golden_hit is None:
        #                 gold_id = hit['subject_id']
        #                 golden_hit = hit
        #             # if the evalue of this hit is smaller than the best one so far, then this is out new best one
        #             elif hit['e_value'] < golden_hit['e_value']:
        #                 gold_id = hit['subject_id']
        #                 golden_hit = hit


        #     print "Gold ID:", gold_id

        #     if gold_id == None:
        #         print "There was no matches in blastp to the golden genes"
        #     else:
        #         # IN A BLASTP CLUSTER, THERE IS A SUBJECT_START AND A QUERY_START
        #         golden_start = golden_hit['subject_start']
        #         gene_start = golden_hit['query_start']

        #         print golden_start
        #         print gene_start

        #          #our gene is too short and we need to move the start upstream
        #         if gene_start == 1 and golden_start == 1:
        #             print "They are already perfectly aligned!"
        #         elif gene_start == 1:
        #             print "Option 1"
        #         #our gene is too long and we need to trim it down
        #         elif golden_start == 1:
        #             print "Option 2"
        #         # right now we do nothing...
        #         else:
        #             print "Neither one starts at 1..."

###############################################################################################################


parser = argparse.ArgumentParser()
parser.add_argument("--clustalo_cutoff", help="Minimum percent identity when clustering", default=32.5)
parser.add_argument("--blastp_cutoff", help="Minimum e-value when clustering", default=1e-35)
args = parser.parse_args()

# Hard coded for dev
start_codons = ['ATG','GTG','TTG'] #CCG
stop_codons = ['TAG', 'TAA', 'TGA']
cluster_id = 17
db = connect_db("geneysis.db")
cluster = get_cluster(db, cluster_id, args)
golden_phage_id = 5
print cluster[4]['phage_id']
#print get_gene(db, 1579)

adjust_cluster(db,cluster,golden_phage_id,start_codons,stop_codons)
'''
db = connect_db("geneysis.db")
phage = get_phage(db, 3)
print phage['seq']
'''