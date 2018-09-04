from Bio import SeqIO
from Bio.Seq import Seq, translate
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA, IUPACProtein
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
import sys, glob, os, sqlite3, subprocess
import json
import numpy as np

def loadProject(wd):
    project = {}
    with open("{d}project.json".format(d=wd),'r') as project_json:
        project = json.load(project_json)
    return project

def writeProject(wd,project):
    with open("{d}project.json".format(d=wd),'w') as project_json:
        json.dump(project,project_json)

def connect_db(db_file):
    return sqlite3.connect(db_file)

def create_db(db_file):
    if os.path.isfile(db_file):
        os.system("rm -rf " + db_file)

    con = connect_db(db_file)

    con.execute('''
        CREATE TABLE `phage` (
            `id`	INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
            `name`	INTEGER NOT NULL,
            `orig_name`	TEXT NOT NULL,
            `organism`	TEXT,
            `definition`	TEXT,
            `golden` INTEGER NOT NULL,
            `seq`	TEXT NOT NULL
        );
        ''')
    con.execute('''
        CREATE TABLE `gene` (
            `id`	INTEGER PRIMARY KEY AUTOINCREMENT,
            `phage_id`	INTEGER,
            `start`	INTEGER NOT NULL,
            `orig_start` INTEGER NOT NULL,
            `end`	INTEGER NOT NULL,
            `orig_end` INTEGER NOT NULL,
            `product`	TEXT NOT NULL,
            `note`	TEXT,
            `locus_tag`	TEXT NOT NULL,
            `old_locus_tag`	TEXT,
            `translation`	TEXT NOT NULL,
            `cluster`	INTEGER,
            `adjusted` INTEGER,
            `rev_comp` INTEGER,
            FOREIGN KEY(`phage_id`) REFERENCES `phage`(`id`),
            FOREIGN KEY(`cluster`) REFERENCES cluster(id)
        );
    ''')
    con.execute('''
        CREATE TABLE `clustalo` (
            `id`	INTEGER PRIMARY KEY AUTOINCREMENT,
            `query_id`	INTEGER NOT NULL,
            `subject_id`	INTEGER NOT NULL,
            `percent_ident`	REAL NOT NULL,
            FOREIGN KEY(`query_id`) REFERENCES gene(id),
            FOREIGN KEY(`subject_id`) REFERENCES gene(id)
        );
    ''')
    con.execute('''
        CREATE TABLE `blastp` (
            `id`	INTEGER PRIMARY KEY AUTOINCREMENT,
            `query_id`	INTEGER NOT NULL,
            `subject_id`	INTEGER NOT NULL,
            `percent_ident`	REAL NOT NULL,
            `e_value`	REAL NOT NULL,
            `query_start`	INTEGER NOT NULL,
            `subject_start`	INTEGER NOT NULL,
            FOREIGN KEY(`query_id`) REFERENCES phage(id),
            FOREIGN KEY(`subject_id`) REFERENCES phage(id)
        );
    ''')
    con.execute('''
        CREATE TABLE `cluster` (
            `id`	INTEGER PRIMARY KEY AUTOINCREMENT,
            `name`	TEXT UNIQUE,
            `adjusted` INTEGER
        );
    ''')
    return con

def insert_phage(db, phage):
    cur = db.cursor()
    cur.execute("INSERT INTO `phage` (name,orig_name,organism,definition,golden,seq) VALUES(?,?,?,?,0,?)",
                (phage.name, phage.name, phage.annotations['organism'], phage.description, str(phage.seq)))
    db.commit()
    return cur.lastrowid

def insert_gene(db, gene, phage_id):
    cur = db.cursor()
    if "old_locus_tag" in gene.qualifiers.keys():
        old_locus = gene.qualifiers["old_locus_tag"][0]
    else:
        old_locus = ""

    if "note" in gene.qualifiers.keys():
        note = gene.qualifiers["note"][0]
    else:
        note = ""

    if gene.strand == 1:
        cur.execute("INSERT INTO `gene` (phage_id,start,orig_start,end,orig_end,product,note,locus_tag,old_locus_tag,translation,adjusted,rev_comp) "
                "VALUES(?,?,?,?,?,?,?,?,?,?,0,0)",
                (phage_id, gene.location.start,gene.location.start, gene.location.end,gene.location.end, gene.qualifiers["translation"][0],
                 note, gene.qualifiers["locus_tag"][0], old_locus, str(gene.qualifiers["translation"][0])))
    else:
        cur.execute("INSERT INTO `gene` (phage_id,start,orig_start,end,orig_end,product,note,locus_tag,old_locus_tag,translation,adjusted,rev_comp) "
                "VALUES(?,?,?,?,?,?,?,?,?,?,0,1)",
                (phage_id, gene.location.start,gene.location.start, gene.location.end,gene.location.end, gene.qualifiers["translation"][0],
                 note, gene.qualifiers["locus_tag"][0], old_locus, str(gene.qualifiers["translation"][0])))
    db.commit()
    return cur.lastrowid

def insert_blastp_hits(db, hits):
    cur = db.cursor()
    cur.executemany("INSERT INTO `blastp` (query_id,subject_id,e_value,query_start,subject_start,percent_ident) "
                    "VALUES(?,?,?,?,?,?)", hits)
    db.commit()

def insert_clustalo_percents(db, query, subject, identity):
    cur = db.cursor()
    cur.execute("INSERT INTO `clustalo` (query_id,subject_id,percent_ident) "
                "VALUES(?,?,?)", (query, subject, identity))
    db.commit()

def get_phage(db, id):
    result = db.execute("SELECT * from `phage` where id = " + str(id))
    if result.arraysize != 1:
        raise Exception("Attempted to find ID " + str(id) + " but found improper number of results.")
    phage = {}
    for row in result:
        phage['id'] = int(row[0])
        phage['name'] = row[1]
        phage['orig_name'] = row[2]
        phage['organism'] = row[3]
        phage['definition'] = row[4]
        phage['golden'] = int(row[5])
        phage['seq'] = row[6]
    return phage

def get_gene(db, id):
    result = db.execute("SELECT * from `gene` where id = " + str(id))
    if result.arraysize != 1:
        raise Exception("Attempted to find ID " + str(id) + " but found improper number of results.")
    gene = {}
    for row in result:
        gene['id'] = int(row[0])
        gene['phage_id'] = int(row[1])
        gene['start'] = int(row[2])
        gene['orig_start'] = int(row[3])
        gene['end'] = int(row[4])
        gene['orig_end'] = int(row[5])
        gene['product'] = row[6]
        gene['note'] = row[6]
        gene['locus_tag'] = row[8]
        gene['old_locus_tag'] = row[9]
        gene['translation'] = row[10]
        if row[11] is None:
            gene['cluster'] = row[11]
        else:
            gene['cluster'] = int(row[11])
        if row[12] == 0:
            gene['adjusted'] = False
        else:
            gene['adjusted'] = True
        if row[13] == 0:
            gene['rev_comp'] = False
        else:
            gene['rev_comp'] = True
    return gene

def get_all_genes(db):
    result = db.execute("SELECT * from `gene`")
    genes = []
    for row in result:
        gene = dict()
        gene['id'] = int(row[0])
        gene['phage_id'] = int(row[1])
        gene['start'] = int(row[2])
        gene['orig_start'] = int(row[3])
        gene['end'] = int(row[4])
        gene['orig_end'] = int(row[5])
        gene['product'] = row[6]
        gene['note'] = row[6]
        gene['locus_tag'] = row[8]
        gene['old_locus_tag'] = row[9]
        gene['translation'] = row[10]
        if row[11] is None:
            gene['cluster'] = row[11]
        else:
            gene['cluster'] = int(row[11])
        if row[12] == 0:
            gene['adjusted'] = False
        else:
            gene['adjusted'] = True
        if row[13] == 0:
            gene['rev_comp'] = False
        else:
            gene['rev_comp'] = True
        genes.append(gene)
    return genes

def get_gene_ids(db):
    result = db.execute("SELECT id from `gene` ORDER BY id")
    ids = [row[0] for row in result]
    # Depreciated
    '''
    ids = []
    for row in result:
        ids.append(row[0])
    '''
    return ids

def get_blastp_hits(db, id, args):

    result = db.execute("SELECT * from `blastp` where query_id = %d and e_value <= %s" %
                        (id, str(args.blastp_cutoff)))
    hits = []
    for row in result:
        hit = {'type': "blastp", 'query_id': int(row[1]), 'subject_id': int(row[2]), 'ident': float(row[3]), 'e_value': float(row[4]),
               'query_start': int(row[5]), 'subject_start': int(row[6])}
        hits.append(hit)
    return hits

def get_all_blastp_hits(db, id):

    result = db.execute("SELECT * from `blastp` where query_id = %d" %
                        (id))
    hits = []
    for row in result:
        hit = {'type': "blastp", 'query_id': int(row[1]), 'subject_id': int(row[2]), 'ident': float(row[3]), 'e_value': float(row[4]),
               'query_start': int(row[5]), 'subject_start': int(row[6])}
        hits.append(hit)
    return hits

def get_clustalo_hits(db, id, args):
    result = db.execute(
        "SELECT * from `clustalo` where query_id = %d or subject_id = %d and percent_ident >= %f" % (id, id, args.clustalo_cutoff))
    hits = []
    for row in result:
        hit = {'type': "clustalo", 'query_id': int(row[1]), 'subject_id': int(row[2]), 'ident': float(row[3])}
        hits.append(hit)
    return hits

def get_all_hits(db, id, args):
    hits = get_clustalo_hits(db, id, args)
    for hit in get_blastp_hits(db, id, args):
        hits.append(hit)
    return hits

#creates a new cluster with a given name
def create_cluster(db, name):
    cur = db.cursor()
    cur.execute("INSERT INTO `cluster` (name,adjusted) VALUES('%s',0)" % name)
    db.commit()
    return cur.lastrowid

#sets the cluster for a gene
def update_gene_cluster(db, gene_id, cluster_id):
    cur = db.cursor()
    cur.execute("UPDATE `gene` set cluster = ? where id = ?",
                (cluster_id, gene_id))
    db.commit()

# gets all gene id's in a given cluster
def get_cluster_genes(db, cluster_id):
    result = db.execute("SELECT id from `gene` where cluster = %d" % cluster_id)
    cluster = [row[0] for row in result]
    return cluster

#gets all the information about each gene in a cluster
def get_cluster(db,cluster_id):
    result = db.execute("SELECT id from `gene` where cluster = %d" % cluster_id)
    cluster = []
    for row in result:
        gene = get_gene(db,row[0])
        #gene['hits'] = get_all_hits(db,gene['id'],args)
        # USE ONLY THE BLASTP HITS
        gene['hits'] = get_all_blastp_hits(db,gene['id'])
        cluster.append(gene)
    return cluster

# returns id of closest cluster, or id of newly created cluster
def get_closest_cluster(db, gene_id, args, i):
    close_clusters = {}
    hits = get_all_hits(db, gene_id, args)
    checked = []
    # go through all hits and if any of them are already in a cluster, it adds the cluster to the "close_clusters" dictionary
    for hit in hits:
        if hit['subject_id'] != gene_id:
            hit_id = hit['subject_id']
        else:
            hit_id = hit['query_id']

        if hit_id not in checked:
            checked.append(hit_id)
            hit_gene = get_gene(db, hit_id)
            if hit_gene['cluster'] is not None:
                if hit_gene['cluster'] not in close_clusters.keys():
                    close_clusters[hit_gene['cluster']] = [hit['ident']]
                else:
                    close_clusters[hit_gene['cluster']].append(hit['ident'])

    if len(close_clusters.keys()) > 0:
        closest_cluster = -1
        closest_identity = 0
        for cluster_id in close_clusters.keys():
            avg = np.mean(close_clusters[cluster_id])
            if avg > closest_identity:
                closest_cluster = cluster_id
                closest_identity = avg
        if closest_cluster == -1:
            raise Exception("Closest cluster not identified.")
        return closest_cluster
    else:
        return create_cluster(db, "cluster_%d" % i)

#returns the golden phages in order, from highest priority to lowest
#for golden numbers, just add 1 to the index of the phage_id in the list
def get_golden_phages(db):
    result = db.execute("select id from phage where golden != 0 order by golden asc;")
    golden_phages = [row[0] for row in result]
    return golden_phages

def get_golden_genes(golden_phages,cluster):
    golden_ids = {}
    for gene in cluster:
        if gene['phage_id'] in golden_phages:
            golden_number = golden_phages.index(gene['phage_id'])
            if golden_number not in  golden_ids.keys():
                golden_ids[golden_number] = []
            golden_ids[golden_number].append(gene['id'])
    return golden_ids

#makes the best possible adjustments for a given cluster, aligning any genes that do not belong to
def adjust_cluster(db,cluster,start_codons,stop_codons):
    #first we need to make a list of all the golden phage proteins that are in this create_cluster
    golden_phages = get_golden_phages(db)
    golden_genes = get_golden_genes(golden_phages, cluster)

    revcomp_start_codons = ['CAT', 'CAC', 'CAA']
    revcomp_stop_codons = ['CTA', 'TTA', 'TCA']

    if len(golden_genes) == 0: #make sure there is at least one gene from the golden phage in the cluster
        return

    for gene in cluster:
        if gene['phage_id'] not in golden_phages and gene['adjusted'] == 0:
            potentialStarts = [] #List to iterate instead of set type# Because if we have muliple golds we will have many different starts to try?
            codonShift = []
            farCodonShift = None
            closeCodonShift = None
            print
            print "New Gene ", gene['id']
            for index, gold_ids in enumerate(golden_genes.values()):
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
                            print(gene["id"])
                        elif gene_start == 1 and blastp_hit['ident'] > 50 and blastp_hit['e_value'] < 1e-9:
                            print "Too Short"
                            ideal_move_distance = np.abs(golden_start - gene_start)
                            newCloseStart, newFarStart, farCodonShift, closeCodonShift = tooShort(db, gene, ideal_move_distance, start_codons, stop_codons, revcomp_start_codons, revcomp_stop_codons)
                            if newCloseStart != None:
                                potentialStarts.append(newCloseStart)
                                codonShift.append(closeCodonShift)
                            if newFarStart != None:
                                potentialStarts.append(newFarStart)
                                codonShift.append(farCodonShift)
                        # our gene is too long and we need to trim it down
                        elif golden_start == 1 and blastp_hit['ident'] > 50 and blastp_hit['e_value'] < 1e-9:
                            print "Too Long"
                            ideal_move_distance = np.abs(gene_start - golden_start)
                            newCloseStart, newFarStart, farCodonShift, closeCodonShift = tooLong(db, gene, ideal_move_distance, start_codons, stop_codons, revcomp_start_codons,revcomp_stop_codons)
                            if newCloseStart != None:
                                potentialStarts.append(newCloseStart)
                                codonShift.append(closeCodonShift)
                            if newFarStart != None:
                                potentialStarts.append(newFarStart)
                                codonShift.append(farCodonShift)
                        # right now we do nothing...
                        else:
                            print "Neither one starts at 1..."
                    else:
                        print "Gene", gene['id'], "has no blastp hit for golden gene", gold_id, ##gene['hits']

            if potentialStarts: # if set is not empty
                print(potentialStarts)
                bestStart = findBestStart(db, gene, potentialStarts, ideal_move_distance, codonShift)
                updateStart(db, gene['id'], bestStart, gene['rev_comp']) #Uncomment when ready

def tooShort(db, gene, ideal_move_distance, start_codons, stop_codons, revcomp_start_codons, revcomp_stop_codons):
    phage = get_phage(db, gene['phage_id'])
    phageGenome = phage['seq']
    currentStart = gene['start']
    if not gene['rev_comp']:
        print "Forward"
        return tooShortForward(db, gene, ideal_move_distance, start_codons, stop_codons)
    elif gene['rev_comp']:
        print "Reverse Compliment"
        return tooShortRevComp(db, gene, ideal_move_distance, revcomp_start_codons, revcomp_stop_codons)

def tooShortRevComp(db, gene, ideal_move_distance, start_codons, stop_codons):
    # Init bestGeneStart
    farBestGeneStart = None
    closeBestGeneStart = None
    farCodonShift = None
    closeCodonShift = None
    # Run through all the potential starts
    for i in xrange(1,ideal_move_distance*2): # doubled to we have equal search space on both sides
        currentStart = gene['end'] + (3 * i) # increase our start 3 at a time
        phage = get_phage(db, gene['phage_id'])
        phageGenome = phage['seq']
        codon = phageGenome[currentStart-3:currentStart]
        ##codon = codon[::-1] # reverse the codon

        if codon in stop_codons:
            print "Found stop codon at {}".format(currentStart)
            break
        if codon in start_codons and i > ideal_move_distance:
            print "far"
            farBestGeneStart = currentStart
            farCodonShift = np.abs(gene['end'] - currentStart)
            break
        elif codon in start_codons and i <= ideal_move_distance:
            print "on or before"
            closeBestGeneStart = currentStart
            closeCodonShift = np.abs(gene['end'] - currentStart)
    return closeBestGeneStart, farBestGeneStart, farCodonShift, closeCodonShift

def tooShortForward(db, gene, ideal_move_distance, start_codons, stop_codons):
    # Init bestGeneStart
    farBestGeneStart = None
    closeBestGeneStart = None
    farCodonShift = None
    closeCodonShift = None
    # Run through all the potential starts
    for i in xrange(1,ideal_move_distance*2): # doubled to we have equal search space on both sides
        currentStart = gene['start'] - (3 * i) # decrease our start 3 at a time
        phage = get_phage(db, gene['phage_id'])
        phageGenome = phage['seq']
        codon = phageGenome[currentStart:currentStart+3]

        if codon in stop_codons:
            print "Found stop codon at {}".format(currentStart)
            break
        if codon in start_codons and i > ideal_move_distance:
            print "far"
            farBestGeneStart = currentStart
            farCodonShift = np.abs(gene['start'] - currentStart)
            break
        elif codon in start_codons and i <= ideal_move_distance:
            print "on or before"
            closeBestGeneStart = currentStart
            closeCodonShift = np.abs(gene['start'] - currentStart)
    return closeBestGeneStart, farBestGeneStart, farCodonShift, closeCodonShift

def tooLong(db, gene, ideal_move_distance, start_codons, stop_codons,revcomp_start_codons,revcomp_stop_codons):
    phage = get_phage(db, gene['phage_id'])
    phageGenome = phage['seq']
    currentStart = gene['start']
    if not gene['rev_comp']:
            print "Forward"
            return tooLongForward(db, gene, ideal_move_distance, start_codons, stop_codons)
    elif gene['rev_comp']:
            print "Reverse Compliment"
            return tooLongRevComp(db, gene, ideal_move_distance, revcomp_start_codons, revcomp_stop_codons)

def tooLongRevComp(db, gene, ideal_move_distance, start_codons, stop_codons):
    # Init bestGeneStart
    farBestGeneStart = None
    closeBestGeneStart = None
    farCodonShift = None
    closeCodonShift = None
    # Run through all the potential starts
    for i in xrange(1,ideal_move_distance*2): # doubled to we have equal search space on both sides
        currentStart = gene['end'] - (3 * i) # decrease our start 3 at a time
        phage = get_phage(db, gene['phage_id'])
        phageGenome = phage['seq']
        codon = phageGenome[currentStart-3:currentStart] # codon is going backward
        ##codon = codon[::-1] # reverse the codon

        if codon in stop_codons:
            print "Found stop codon at {}".format(currentStart)
            break
        if codon in start_codons and i > ideal_move_distance:
            print "far"
            farBestGeneStart = currentStart
            farCodonShift = np.abs(gene['end'] - currentStart)
            break
        elif codon in start_codons and i <= ideal_move_distance:
            print "on or before"
            closeBestGeneStart = currentStart
            closeCodonShift = np.abs(gene['end'] - currentStart)
    return closeBestGeneStart, farBestGeneStart, farCodonShift, closeCodonShift

def tooLongForward(db, gene, ideal_move_distance, start_codons, stop_codons):
    # Init bestGeneStart
    farBestGeneStart = None
    closeBestGeneStart = None
    farCodonShift = None
    closeCodonShift = None
    # Run through all the potential starts
    for i in xrange(1,ideal_move_distance*2): # doubled to we have equal search space on both sides
        currentStart = gene['start'] + (3 * i) # increase our start 3 at a time
        phage = get_phage(db, gene['phage_id'])
        phageGenome = phage['seq']
        codon = phageGenome[currentStart:currentStart+3] # codon is going forward

        if codon in stop_codons:
            print "Found stop codon at {}".format(currentStart)
            break
        if codon in start_codons and i > ideal_move_distance:
            print "far"
            farBestGeneStart = currentStart
            farCodonShift = np.abs(gene['start'] - currentStart)
            break
        elif codon in start_codons and i <= ideal_move_distance:
            print "on or before"
            closeBestGeneStart = currentStart
            closeCodonShift = np.abs(gene['start'] - currentStart)
    return closeBestGeneStart, farBestGeneStart, farCodonShift, closeCodonShift

def updateStart(db, gene_id, newStart, rev_comp):
    cur = db.cursor()
    if rev_comp == False:
        cur.execute("UPDATE gene SET start = " + str(newStart) + ", adjusted = 1 WHERE id = " + str(gene_id))
    else:
        cur.execute("UPDATE gene SET end = " + str(newStart) + ", adjusted = 1 WHERE id = " + str(gene_id))
    db.commit()

def findBestStart(db, gene, potentialStarts, ideal_move_distance, codonShift):
    if gene['rev_comp']:
        check = gene['end']
    else:
        check = gene['start']
    imd = ideal_move_distance

    for item in potentialStarts:
        if not isinstance(item, int):
            potentialStarts.pop(potentialStarts.index(item))
    diffs = [np.abs(s - imd) for s in codonShift]
    return potentialStarts[np.argmin(diffs)]

def fillGap(db, cluster, phage):
        golden_phages = get_golden_phages(db)
        golden_genes = get_golden_genes(golden_phages, cluster)

        tempGB = open("gapGenes.gb", "w")
        new_phage = get_phage(db, phage)
        new_phage_seq = new_phage["seq"]
        new_phage_seq = Seq(new_phage_seq, IUPACUnambiguousDNA())
        record = SeqRecord(new_phage_seq, id = "0", name = "Temp_Genome", description = "Temp_file")

        if len(golden_genes) == 0: #make sure there is at least one gene from the golden phage in the cluster
            return

        if len(cluster) == len(golden_genes):
            for gene in cluster:
                querySeq = gene["product"]
                query = open("query.FASTA", "w")
                query.write(">query\n{}".format(querySeq))
                query.close()
                frame1, frame2, frame3 = getReadingFrames(gene, db, phage)
                print
                print(gene["locus_tag"])

                counter = 0
                resultProtein = None
                while counter < 3:
                    if counter == 0:
                        subjectProtein = str(frame1)
                    elif counter == 1:
                        subjectProtein = str(frame2)
                    elif counter == 2:
                        subjectProtein = str(frame3)
                    subject = open("subject.FASTA", "w")
                    subject.write(">subject\n{}".format(subjectProtein))
                    subject.close()
                    resultProtein = blastGap(gene, "query.FASTA", "subject.FASTA")
                    if resultProtein != None:
                        break
                    else:
                        counter += 1

                numAminoAcids = len(resultProtein)
                resultProtein = Seq(resultProtein, IUPACProtein())
                newFeature = makeNewGene(numAminoAcids, resultProtein, frame1, frame2, frame3, gene, new_phage_seq)
                print(newFeature)
                record.features.append(newFeature)

        SeqIO.write(record, tempGB, "genbank")
        tempGB.close()

        genome = SeqIO.read("gapGenes.gb", "genbank")
        features = genome.features
        for feature in features:
            if feature.type == "CDS":
                if feature.qualifiers['translation'][0] == '-':
                    feature.qualifiers['translation'][0] = feature.qualifiers['translation'][0][1:]
                gene_id = insert_gene(db, feature, phage)

def blastGap(gene, query, subject):
    e_value = 1e-5
    gapBlast = NcbiblastpCommandline(query = query, subject = subject, outfmt=5, out='gap.xml')
    stdout, stderr = gapBlast()
    result_handle = open('gap.xml')
    blast_record = NCBIXML.read(result_handle)
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < e_value:
                print("****ALIGNMENT****")
                print("sequence:", alignment.title)
                print("length", alignment.length)
                print("e-value:", hsp.expect)
                print(hsp.query[0:75] + "...")
                print(hsp.match[0:75] + "...")
                print(hsp.sbjct[0:75] + "...")
                return hsp.sbjct

def getReadingFrames(gene, db, phage):
    if gene['rev_comp'] != 1:
        new_phage = get_phage(db, phage)
        new_phage_seq = new_phage['seq']
        seq = Seq(new_phage_seq, IUPACUnambiguousDNA())
        frame1 = seq.translate()
        frame2 = seq[1:(len(seq) - 1)]
        frame2 = frame2.translate()
        frame3 = seq[2:(len(seq) - 1)]
        frame3 = frame3.translate()
    else:
        new_phage = get_phage(db, phage)
        new_phage_seq = new_phage['seq']
        seq = Seq(new_phage_seq, IUPACUnambiguousDNA())
        seq = seq.reverse_complement()
        frame1 = seq.translate()
        frame2 = seq[1:(len(seq) - 1)]
        frame2 = frame2.translate()
        frame3 = seq[2:(len(seq) - 1)]
        frame3 = frame3.translate()
    return frame1, frame2, frame3

def makeNewGene(length, protein, frame1, frame2, frame3, gene, seq):
    locusTag = "LG_{}".format(gene["locus_tag"])
    counter = 0
    while counter < 3:
        if counter == 0:
            start = frame1.find(protein)
            if start != -1:
                break
            else:
                counter += 1
        elif counter == 1:
            start = frame2.find(protein)
            if start != -1:
                break
            else:
                counter += 1
        elif counter == 2:
            start = frame3.find(protein)
            if start != -1:
                break
            else:
                counter += 1
    if start != -1:
        nucStart = counter + 3 * start
        nucEnd = nucStart + 3 * length
        if gene["rev_comp"] == 1:
            geneSeq = seq[nucStart:nucEnd]
            geneSeq = geneSeq.reverse_complement()
            seq = seq.reverse_complement()
            nucStart = seq.find(geneSeq)
            nucEnd = nucStart + 3 * length
            newFeature = SeqFeature(FeatureLocation(nucStart, nucEnd, strand = -1), type = "CDS")
            newFeature.qualifiers["locus_tag"] = locusTag
            newFeature.qualifiers["translation"] = protein
        else:
            newFeature = SeqFeature(FeatureLocation(nucStart, nucEnd, strand = 1), type = "CDS")
            newFeature.qualifiers["locus_tag"] = locusTag
            newFeature.qualifiers["translation"] = protein
    return newFeature
