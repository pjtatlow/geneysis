from functions import *
import argparse

def updateStart(db, gene_id, newStart, rev_comp):
    cur = db.cursor()
    if !rev_comp:
        cur.execute("UPDATE gene SET start = " + str(newStart) + " WHERE id = " + str(gene_id))
    else:
        cur.execute("UPDATE gene SET end = " + str(newStart) + " WHERE id = " + str(gene_id))
    db.commit()


def findBestStart(db, gene, potentialStarts, ideal_move_distance):
    print potentialStarts
    if gene['rev_comp']:
        check = gene['end']
    else:
        check = gene['start']
    imd = ideal_move_distance * 3

    diffs = [np.abs(s - imd) for s in potentialStarts]
    return potentialStarts[np.argmax(diffs)]


def tooLongForwardNormal(db, gene, ideal_move_distance, start_codons, stop_codons):
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


def tooLongBackwardNormal(db, gene, ideal_move_distance, start_codons, stop_codons):
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


def tooLongForwardRevComp(db, gene, ideal_move_distance, start_codons, stop_codons):
    # Init bestGeneStart
    farBestGeneStart = None
    closeBestGeneStart = None
    # Run through all the potential starts
    for i in xrange(1,ideal_move_distance*2): # doubled to we have equal search space on both sides
        currentStart = gene['end'] + (3 * i) - 2 # increase our start 3 at a time
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


def tooLongBackwardRevComp(db, gene, ideal_move_distance, start_codons, stop_codons):
    # Init bestGeneStart
    farBestGeneStart = None
    closeBestGeneStart = None
    # Run through all the potential starts
    for i in xrange(1,ideal_move_distance*2): # doubled to we have equal search space on both sides
        currentStart = gene['end'] - (3 * i) - 2 # decrease our start 3 at a time
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


def tooLong(db, gene, ideal_move_distance, start_codons, stop_codons,revcomp_start_codons,revcomp_stop_codons):
    if gene['rev_comp'] == False:
        if gene['start'] < gene['end']:
            print "Forward"
            return tooLongForwardNormal(db, gene, ideal_move_distance, start_codons, stop_codons)
        elif gene['start'] > gene['end']:
            print "Backward"
            return tooLongBackwardNormal(db, gene, ideal_move_distance, start_codons, stop_codons)
    else:
        if gene['start'] < gene['end']:
            print "Forward"
            return tooLongForwardRevComp(db, gene, ideal_move_distance, revcomp_start_codons, revcomp_stop_codons)
        elif gene['start'] > gene['end']:
            print "Backward"
            return tooLongBackwardRevComp(db, gene, ideal_move_distance, revcomp_start_codons, revcomp_stop_codons)



def tooShortForwardNormal(db, gene, ideal_move_distance, start_codons, stop_codons):
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


def tooShortBackwardNormal(db, gene, ideal_move_distance, start_codons, stop_codons):
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


def tooShortForwardRevComp(db, gene, ideal_move_distance, start_codons, stop_codons):
    # Init bestGeneStart
    farBestGeneStart = None
    closeBestGeneStart = None
    # Run through all the potential starts
    for i in xrange(1,ideal_move_distance*2): # doubled to we have equal search space on both sides
        currentStart = gene['end'] - (3 * i) - 2 # decrease our start 3 at a time
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


def tooShortBackwardRevComp(db, gene, ideal_move_distance, start_codons, stop_codons):
    # Init bestGeneStart
    farBestGeneStart = None
    closeBestGeneStart = None
    # Run through all the potential starts
    for i in xrange(1,ideal_move_distance*2): # doubled to we have equal search space on both sides
        currentStart = gene['end'] + (3 * i) - 2 # increase our start 3 at a time
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
    if gene['rev_comp'] == False:
        if gene['start'] < gene['end']:
            print "Forward"
            return tooShortForwardNormal(db, gene, ideal_move_distance, start_codons, stop_codons)
        elif gene['start'] > gene['end']:
            print "Backward"
            return tooShortBackwardNormal(db, gene, ideal_move_distance, start_codons, stop_codons)
    else:
        if gene['start'] < gene['end']:
            print "Forward"
            return tooShortForwardRevComp(db, gene, ideal_move_distance, revcomp_start_codons, revcomp_stop_codons)
        elif gene['start'] > gene['end']:
            print "Backward"
            return tooShortBackwardRevComp(db, gene, ideal_move_distance, revcomp_start_codons, revcomp_stop_codons)



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
            potentialStarts = set() # Because if we have muliple golds we will have many different starts to try?
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
                            newCloseStart, newFarStart = tooShort(db, gene, ideal_move_distance, start_codons, stop_codons,revcomp_start_codons,revcomp_stop_codons)
                            if newCloseStart != None:
                                potentialStarts.add(newCloseStart)
                            if newFarStart != None:
                                potentialStarts.add(newFarStart)
                        # our gene is too long and we need to trim it down
                        elif golden_start == 1:
                            print "Too Long"
                            ideal_move_distance = gene_start - golden_start
                            newCloseStart, newFarStart = tooLong(db, gene, ideal_move_distance, start_codons, stop_codons,revcomp_start_codons,revcomp_stop_codons)

                            ######## Just to debug #####################
                            print "Starts:", newCloseStart, newFarStart
                            print "Original:", gene['start']
                            ###################################

                            if newCloseStart != None:
                                potentialStarts.add(newCloseStart)
                            if newFarStart != None:
                                potentialStarts.add(newFarStart)
                        # right now we do nothing...
                        else:
                            print "Neither one starts at 1..."


                    else:
                        print "Gene", gene['id'], "has no blastp hit for golden gene", gold_id, gene['hits']

            if potentialStarts: # if set is not empty
                bestStart = findBestStart(db, gene, potentialStarts, ideal_move_distance)
                updateStart(db, gene['id'], bestStart, gene['rev_comp']) Uncomment when ready

###############################################################################################################


parser = argparse.ArgumentParser()
parser.add_argument("--clustalo_cutoff", help="Minimum percent identity when clustering", default=32.5)
parser.add_argument("--blastp_cutoff", help="Minimum e-value when clustering", default=1e-35)
args = parser.parse_args()

# Hard coded for dev
start_codons = ['ATG','GTG','TTG'] #CCG
revcomp_start_codons = ['CAT', 'CAC', 'CAA']
stop_codons = ['TAG', 'TAA', 'TGA']
revcomp_stop_codons = ['CTA', 'TTA', 'TCA']
cluster_id = 17
db = connect_db("geneysis.db")
cluster = get_cluster(db, cluster_id, args)
golden_phage_id = 5

adjust_cluster(db,cluster,golden_phage_id,start_codons,stop_codons,revcomp_start_codons,revcomp_stop_codons)
