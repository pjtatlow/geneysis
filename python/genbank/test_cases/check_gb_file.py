from Bio import SeqIO
from sys import argv

gb_file = argv[1]
# Note: only one record in our gb files, so the loop runs only once
for record in SeqIO.parse(gb_file, 'genbank'):
    print record.name
    print record.annotations['organism']
    print record.description
    print len(record.seq)
    for feature in record.features:
        if feature.type == 'CDS':
            gene = feature
            print gene.qualifiers['translation']
