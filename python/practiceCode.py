from functions import *
import argparse


#makes the best possible adjustments for a given cluster, aligning any genes that do not belong to 
def adjust_cluster(db,cluster,golden_phage_id,start_codons):
    #first we need to make a list of all the golden phage proteins that are in this create_cluster
    golden_genes = get_golden_genes(cluster,golden_phage_id)

    print golden_genes
    
    if len(golden_genes) == 0: #make sure there is at least one gene from the golden phage in the cluster
        return
    
    for gene in cluster:
        print
        print "New Gene"
        if not gene['adjusted'] and gene['phage_id'] != golden_phage_id: # only adjust non-golden genes that have not been adjusted
            golden_hit = None
            gold_id = None
            gold_in_thing = False
            # find the closest golden gene blastp his
            for hit in gene['hits']: 
                print hit
                print
                # WHY IS THIS ONLY LOOKING FOR BLASTP? WHAT INFO IS IMPORTANT FOR CLUSTALO?
                # MY HARD CODED EXAMPLE ONLY HAS CLUSTALO STUFF
                # make sure we're dealing with a blastp hit is to a golden gene
                if hit['type'] == "blastp" and hit['subject_id'] in golden_genes:
                    # first golden gene blastp hit
                    if golden_hit is None:
                        gold_id = hit['subject_id']
                        golden_hit = hit
                    # if the evalue of this hit is smaller than the best one so far, then this is out new best one
                    elif hit['e_value'] < golden_hit['e_value']:
                        golden_hit = hit
                if hit['type'] == "clustalo" and hit['subject_id'] in golden_genes: # make sure we're dealing with a blastp hit is to a golden gene
                    # first golden gene blastp hit
                    if golden_hit is None:
                        golden_hit = hit
                    # if the ident of this hit is smaller than the best one so far, then this is out new best one
                    elif hit['ident'] > golden_hit['ident']: 
                        golden_hit = hit
                if hit['subject_id'] in golden_genes:
                    gold_in_thing = True
            print "There is a golden hit:", gold_in_thing
            print "---------------------------------------------------------"

            # golden_hit DOESN'T HAVE 'subject_start' OR 'query_start'...
            #golden_start = golden_hit['subject_start']
            #gene_start = golden_hit['query_start']

            # THIS THROWS AN ERROR. HOW COULD GOLD_G NOT EXIST?
            gold_g = get_gene(db, gold_id)
            golden_start = gold_g['start']
            gene_start = gene['start']

            print golden_start
            print gene_start

             #our gene is too short and we need to move the start upstream
            if gene_start == 1:
                None
            #our gene is too long and we need to trim it down
            elif golden_start == 1: 
                None
            # right now we do nothing...
            else:
                None

parser = argparse.ArgumentParser()
parser.add_argument("--clustalo_cutoff", help="Minimum percent identity when clustering", default=32.5)
parser.add_argument("--blastp_cutoff", help="Minimum e-value when clustering", default=1e-50)
args = parser.parse_args()

# Hard coded for dev
start_codons = ['P'] #CCG
cluster_id = 7
db = connect_db("geneysis.db")
cluster = get_cluster(db, cluster_id, args)
golden_phage_id = 7

print get_gene(db, 1579)

adjust_cluster(db,cluster,golden_phage_id,start_codons)