from functions import *
import argparse


#makes the best possible adjustments for a given cluster, aligning any genes that do not belong to 
def adjust_cluster(db,cluster,golden_phage_id,start_codons):
    #first we need to make a list of all the golden phage proteins that are in this create_cluster
    golden_genes = get_golden_genes(cluster,golden_phage_id)

    #print golden_genes
    
    if len(golden_genes) == 0: #make sure there is at least one gene from the golden phage in the cluster
        return
    
    for gene in cluster:
        print
        print "New Gene"
        if not gene['adjusted'] and gene['phage_id'] != golden_phage_id: # only adjust non-golden genes that have not been adjusted
            golden_hit = None
            gold_id = None
            # find the closest golden gene blastp his
            for hit in gene['hits']:
                print hit
                # make sure we're dealing with a blastp hit is to a golden gene
                if hit['type'] == "blastp" and hit['subject_id'] in golden_genes:
                    # first golden gene blastp hit
                    if golden_hit is None:
                        gold_id = hit['subject_id']
                        golden_hit = hit
                    # if the evalue of this hit is smaller than the best one so far, then this is out new best one
                    elif hit['e_value'] < golden_hit['e_value']:
                        gold_id = hit['subject_id']
                        golden_hit = hit


            print "Gold ID:", gold_id

            if gold_id == None:
                print "There was no matches in blastp to the golden genes"
            else:
                # IN A BLASTP CLUSTER, THERE IS A SUBJECT_START AND A QUERY_START
                golden_start = golden_hit['subject_start']
                gene_start = golden_hit['query_start']

                print golden_start
                print gene_start

                 #our gene is too short and we need to move the start upstream
                if gene_start == 1 and golden_start == 1:
                    print "They are already perfectly aligned!"
                elif gene_start == 1:
                    print "Option 1"
                #our gene is too long and we need to trim it down
                elif golden_start == 1:
                    print "Option 2"
                # right now we do nothing...
                else:
                    print "Neither one starts at 1..."

###############################################################################################################



parser = argparse.ArgumentParser()
parser.add_argument("--clustalo_cutoff", help="Minimum percent identity when clustering", default=32.5)
parser.add_argument("--blastp_cutoff", help="Minimum e-value when clustering", default=0)
args = parser.parse_args()

# Hard coded for dev
start_codons = ['P'] #CCG
cluster_id = 7
db = connect_db("geneysis.db")
cluster = get_cluster(db, cluster_id, args)
golden_phage_id = 5
print cluster[4]['phage_id']
#print get_gene(db, 1579)

adjust_cluster(db,cluster,golden_phage_id,start_codons)