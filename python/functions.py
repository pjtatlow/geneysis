from Bio import SeqIO
import sys, glob, os, sqlite3, subprocess
import json


def mean(x):
    total = 0.0
    for y in x:
        total += float(y)
    return total / float(len(x))

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

    con = sqlite3.connect(db_file)

    con.execute('''
        CREATE TABLE `phage` (
            `id`	INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
            `name`	INTEGER NOT NULL,
            `orig_name`	TEXT NOT NULL,
            `organism`	TEXT,
            `definition`	TEXT,
            `seq`	TEXT NOT NULL
        );
        ''')
    con.execute('''
        CREATE TABLE `gene` (
            `id`	INTEGER PRIMARY KEY AUTOINCREMENT,
            `phage_id`	INTEGER,
            `start`	INTEGER NOT NULL,
            `end`	INTEGER NOT NULL,
            `product`	TEXT NOT NULL,
            `note`	TEXT,
            `locus_tag`	TEXT NOT NULL UNIQUE,
            `old_locus_tag`	TEXT,
            `translation`	TEXT NOT NULL,
            `cluster`	INTEGER,
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
            `name`	TEXT UNIQUE
        );
    ''')
    return con


def insert_phage(db, phage):
    cur = db.cursor()
    cur.execute("INSERT INTO `phage` (name,orig_name,organism,definition,seq) VALUES(?,?,?,?,?)",
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

    cur.execute("INSERT INTO `gene` (phage_id,start,end,product,note,locus_tag,old_locus_tag,translation) "
                "VALUES(?,?,?,?,?,?,?,?)",
                (phage_id, gene.location.start, gene.location.end, gene.qualifiers["product"][0],
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


def get_gene(db, id):
    result = db.execute("SELECT * from `gene` where id = " + str(id))
    if result.arraysize != 1:
        raise Exception("Attempted to find ID " + str(id) + " but found improper number of results.")
    gene = {}
    for row in result:
        gene['id'] = row[0]
        gene['phage'] = row[1]
        gene['start'] = row[2]
        gene['end'] = row[3]
        gene['product'] = row[4]
        gene['note'] = row[5]
        gene['locus_tag'] = row[6]
        gene['old_locus_tag'] = row[7]
        gene['translation'] = row[8]
        gene['cluster'] = row[9]
    return gene


def get_all_genes(db):
    result = db.execute("SELECT * from `gene`")
    genes = []
    for row in result:
        gene = dict()
        gene['id'] = row[0]
        gene['phage'] = row[1]
        gene['start'] = row[2]
        gene['end'] = row[3]
        gene['product'] = row[4]
        gene['note'] = row[5]
        gene['locus_tag'] = row[6]
        gene['old_locus_tag'] = row[7]
        gene['translation'] = row[8]
        gene['cluster'] = row[9]
        genes.append(gene)
    return genes

def get_gene_ids(db):
    result = db.execute("SELECT id from `gene` ORDER BY id")
    ids = []
    for row in result:
        ids.append(row[0])
    return ids


def get_blastp_hits(db, id, args):
    result = db.execute("SELECT * from `blastp` where query_id = %d and e_value <= %f" %
                        (id, args.blastp_cutoff))
    hits = []
    for row in result:
        hit = {'type': "blastp", 'query_id': row[1], 'subject_id': row[2], 'ident': row[3], 'e_value': row[4],
               'query_start': row[5], 'subject_start': row[6]}
        hits.append(hit)
    return hits


def get_clustalo_hits(db, id, args):
    result = db.execute(
        "SELECT * from `clustalo` where query_id = %d or subject_id = %d and percent_ident >= %f" % (id, id, args.clustalo_cutoff))
    hits = []
    for row in result:
        hit = {'type': "clustalo", 'query_id': row[1], 'subject_id': row[2], 'ident': row[3]}
        hits.append(hit)
    return hits


def get_all_hits(db, id, args):
    hits = get_clustalo_hits(db, id, args)
    for hit in get_blastp_hits(db, id, args):
        hits.append(hit)
    return hits


def create_cluster(db, name):
    cur = db.cursor()
    cur.execute("INSERT INTO `cluster` (name) VALUES('%s')" % name)
    db.commit()
    return cur.lastrowid


def update_gene_cluster(db, gene_id, cluster_id):
    cur = db.cursor()
    cur.execute("UPDATE `gene` set cluster = ? where id = ?",
                (cluster_id, gene_id))
    db.commit()


def get_cluster_genes(db, cluster_id):
    result = db.execute("SELECT id from `gene` where cluster = %d" % cluster_id)
    cluster = []
    for row in result:
        cluster.append(row[0])
    return cluster


# returns id of closest cluster, or id of newly created cluster
def get_closest_cluster(db, gene_id, args, i):
    close_clusters = {}
    hits = get_all_hits(db, gene_id, args)
    checked = []
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
            avg = mean(close_clusters[cluster_id])
            if avg > closest_identity:
                closest_cluster = cluster_id
                closest_identity = avg
        if closest_cluster == -1:
            raise Exception("Closest cluster not identified.")
        return closest_cluster
    else:
        return create_cluster(db, "cluster_%d" % i)

