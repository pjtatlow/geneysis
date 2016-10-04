#!/usr/bin/env python
from functions import *
import argparse
import datetime
import time

settings = {
    "wd": "/opt/geneysis/",
    "db": "db/geneysis.db",
    "clustalo_cutoff": 32.5,
    "blastp_cutoff": 1e-50
}


parser = argparse.ArgumentParser()
parser.add_argument("task", help="Specify which task you want run")
parser.add_argument("--wd", help="Directory containing all required files")
parser.add_argument("--file", help="Genbank file to work with")
parser.add_argument("--clustalo_cutoff", help="Minimum percent identity when clustering", default=32.5)
parser.add_argument("--blastp_cutoff", help="Minimum e-value when clustering", default=1e-50)

args = parser.parse_args()

if args.wd is None:
    print "NO WD GIVEN"
    sys.exit(1);
elif args.wd[len(args.wd)-1] != "/":
    args.wd += "/"

#print args

if args.task == "create":
    for directory in ["blastp", "clustalo", "db", "fasta", "genbank"]:
        os.system("mkdir -p {wd}{dir}".format(wd=args.wd, dir=directory))
    create_db(args.wd + "db/geneysis.db")
    dirs = args.wd.split('/')
    name = dirs[len(dirs) - 2]
    project = {"name":name,"path":args.wd,"created":time.time(),"updated":time.time(),"phages":[],"events":[]}

    writeProject(args.wd, project)

    print json.dumps(project)
    
elif args.task == "import":

    if args.file is None:
        print >> sys.stderr, "No file given"
        sys.exit(1)


    db = connect_db(args.wd + "db/geneysis.db")

    abs_path = os.path.abspath(args.file)

    for phage in SeqIO.parse(abs_path, "genbank"):
        os.system("cp " + abs_path.replace(" ", "\ ") + " " + args.wd + "genbank/" + phage.name + ".gbk")

        phage_id = insert_phage(db, phage)

        for gene in phage.features:
            if gene.type == "CDS":
                if gene.qualifiers['translation'][0][0] == '-':
                    gene.qualifiers['translation'][0] = gene.qualifiers['translation'][0][1:]
                gene_id = insert_gene(db, gene, phage_id)
                
        project = loadProject(args.wd)
 
        project['phages'].append(phage.name)
        
        writeProject(args.wd, project)
        
        print json.dumps(project)
        
    db.close()

elif args.task == "create_fasta":
    with open(args.wd + "fasta/geneysis.fasta", 'w') as fasta:
        db = connect_db(args.wd + "db/geneysis.db")

        genes = get_all_genes(db)
        print len(genes), "genes"
        for gene in genes:
            fasta.write(">{id}\n{translation}\n".format(id=gene['id'], translation=gene['translation']))

        db.close()

elif args.task == "blast":

    db = connect_db(args.wd + "db/geneysis.db")

    os.system("makeblastdb -in " + args.wd + "fasta/geneysis.fasta -out " + args.wd +
              "blastp/proteinsdb -dbtype prot")
    print "Blasting all proteins"
    blastp = subprocess.Popen(
        "blastp -db /opt/geneysis/blastp/proteinsdb -query /opt/geneysis/fasta/geneysis.fasta -outfmt '6 qseqid sseqid evalue qstart sstart pident' ",
        shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    blast_results, err = blastp.communicate()
    blast_results = blast_results.split('\n')
    hits = []
    i = 0
    for hit in blast_results:
        hit = hit.split('\t')
        if len(hit) == 6 and hit[0] != hit[1]:
            hit[0] = int(hit[0])
            hit[1] = int(hit[1])
            hit[2] = float(hit[2])
            hit[3] = int(hit[3])
            hit[4] = int(hit[4])
            hit[5] = float(hit[5])
            hits.append(hit)
            i += 1
    insert_blastp_hits(db, hits)

    db.close()


elif args.task == "clustalo":
    os.system("clustalo -i " + args.wd + "fasta/geneysis.fasta --outfmt=clu --distmat-out=" + args.wd +
              "clustalo/percent.id --full --percent-id --force -o " + args.wd + "clustalo/align.clu")

    db = connect_db(args.wd + "db/geneysis.db")

    with open(args.wd + "clustalo/percent.id", 'r') as clustalo:
        clustalo.readline()  # take out header
        for line in clustalo:
            line = line.strip().split()
            gene = int(line[0])
            for x in xrange(gene + 1, len(line)):
                if float(line[x]) >= args.clustalo_cutoff:
                    insert_clustalo_percents(db, gene, x, line[x])

    db.close()

elif args.task == "cluster":
    db = sqlite3.connect(args.wd + "db/geneysis.db")

    ids = get_gene_ids(db)
    i = 1
    for id in ids:
        cluster = get_closest_cluster(db, id, args, i)
        if cluster == i:
            i += 1
        update_gene_cluster(db, id, cluster)

    db.close()

else:
    print args.task, "not a valid task!"
    sys.exit(1)

# if "-in" in sys.argv:
#     index = sys.argv.index("-in")
#     file_path = sys.argv[index + 1]
#     if file_path[len(file_path) - 1] != "/":
#         file_path += "/"
#     input_files = glob.glob(file_path + "*")
#
# else:
#     print "Missing Genbank file path.\nUsage: geneysis.py --in path-to-genbank-files -out"
#     sys.exit(1)
#
# if "-out" in sys.argv:
#     index = sys.argv.index("-out")
#     out_path = sys.argv[index + 1]
#
# fasta = open(settings['wd'] + "fasta/geneysis.fasta", 'w')
#
# db = create_db(settings)
#
# for name in input_files:
#     for phage in SeqIO.parse(name, "genbank"):
#         os.system("cp " + name.replace(" ", "\ ") + " " + settings["wd"] + "genbank/" + phage.name + ".gbk")
#
#         phage_id = insert_phage(db, phage)
#
#         for gene in phage.features:
#             if gene.type == "CDS":
#                 if gene.qualifiers['translation'][0][0] == '-':
#                     gene.qualifiers['translation'][0] = gene.qualifiers['translation'][0][1:]
#                 gene_id = insert_gene(db, gene, phage_id)
#                 fasta.write(">" + str(gene_id) + "\n" + gene.qualifiers['translation'][0] + "\n")
# fasta.flush()
# fasta.close()
#
# os.system("makeblastdb -in " + settings["wd"] + "fasta/geneysis.fasta -out " + settings["wd"] +
#           "blastp/proteinsdb -dbtype prot")
# print "Blasting all proteins"
# blastp = subprocess.Popen(
#     "blastp -db /opt/geneysis/blastp/proteinsdb -query /opt/geneysis/fasta/geneysis.fasta -outfmt '6 qseqid sseqid evalue qstart sstart pident' ",
#     shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
# blast_results, err = blastp.communicate()
# blast_results = blast_results.split('\n')
# hits = []
# i = 0
# for hit in blast_results:
#     hit = hit.split('\t')
#     if len(hit) == 6 and hit[0] != hit[1]:
#         hit[0] = int(hit[0])
#         hit[1] = int(hit[1])
#         hit[2] = float(hit[2])
#         hit[3] = int(hit[3])
#         hit[4] = int(hit[4])
#         hit[5] = float(hit[5])
#         hits.append(hit)
#         i += 1
# print "Inserting blast hits"
# insert_blastp_hits(db, hits)
#
#
# # os.system("clustalo -i " + settings["wd"] + "fasta/geneysis.fasta --outfmt=clu --distmat-out=" + settings["wd"] +
# #           "clustalo/percent.id --full --percent-id --force -o " + settings["wd"] + "clustalo/align.clu")
# print "Parsing Clustalo output"
# with open(settings["wd"] + "clustalo/percent.id", 'r') as clustalo:
#     clustalo.readline()  # take out header
#     for line in clustalo:
#         line = line.strip().split()
#         gene = int(line[0])
#         for x in xrange(gene + 1, len(line)):
#             if float(line[x]) >= settings['clustalo_cutoff']:
#                 insert_clustalo_percents(db,gene,x,line[x])

# print "Finished loading database"
#
# db = sqlite3.connect(settings['wd'] + settings['db'])
# print "Clustering proteins"
# ids = get_gene_ids(db)
# i = 1
# for id in ids:
#     cluster = get_closest_cluster(db,id,settings,i)
#     if cluster == i:
#         i += 1
#     update_gene_cluster(db, id, cluster)
# db.close()
