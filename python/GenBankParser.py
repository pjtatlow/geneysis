import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import FeatureLocation, SeqFeature
import random

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="input GenBank file")
parser.add_argument("-g", "--genBankOut", help="output GenBank file")
parser.add_argument("-c", "--csvOutput", help="output text file")

args = parser.parse_args()

startCodons = ["ATG", "GTG", "TTG"]
revStartCodons = ["CAT", "CAC", "CAA"]

def randShiftForward(seqLoc, genomeSeq):
    upOrDown = random.randint(1,2)
    if upOrDown == 2:
        start = int(seqLoc.start)
        upperRange = start - 150
        searchSeq = str(genomeSeq[upperRange:start])
        newStartCodon = -1
        for codon in startCodons:
            if newStartCodon == -1:
                newStartCodon = searchSeq.find(codon)
        newStart = upperRange + newStartCodon
    else:
        start = int(seqLoc.start) + 2
        lowerRange = start + 150
        searchSeq = str(genomeSeq[start:lowerRange])
        newStartCodon = -1
        for codon in startCodons:
            if newStartCodon == -1:
                newStartCodon = searchSeq.find(codon)
        newStart = start + newStartCodon
    newLoc = FeatureLocation(newStart, seqLoc.end, strand=seqLoc.strand)
    return newLoc

def randShiftReverse(seqLoc, genomeSeq):
    upOrDown = random.randint(1,2)
    if upOrDown == 2:
        start = int(seqLoc.end) - 2
        upperRange = start - 150
        searchSeq = str(genomeSeq[upperRange:start])
        newStartCodon = -1
        for codon in revStartCodons:
            if newStartCodon == -1:
                newStartCodon = searchSeq.find(codon)
        newStart = upperRange + newStartCodon + 3
    else:
        start = int(seqLoc.end)
        lowerRange = start + 150
        searchSeq = str(genomeSeq[start:lowerRange])
        newStartCodon = -1
        for codon in revStartCodons:
            if newStartCodon == -1:
                newStartCodon = searchSeq.find(codon)
        newStart = start + newStartCodon + 3
    newLoc = FeatureLocation(seqLoc.start, newStart, strand=seqLoc.strand)
    return newLoc

outputFile1 = open(args.csvOutput, "w")
outputFile2 = open(args.genBankOut, "w")
genBankFile = args.input
genome = SeqIO.read(genBankFile, "genbank")
seq = genome.seq
features = genome.features
record = SeqRecord(seq, id = genome.id, name = genome.name, description = genome.description)
outputFile1.write("Gene ID\tOld Start\tNew Start\n")
newGene = None
for feature in features:
    quals = feature.qualifiers
    if feature.type == "gene":
        local = feature.location
        geneSeq = seq[local.start:local.end]
        if local.strand == 1:
            newGene = randShiftForward(local, seq)
        elif local.strand == -1:
            newGene = randShiftReverse(local, seq)
        name = str(feature.qualifiers["locus_tag"])
        featureWO = SeqFeature(newGene, type = "gene", qualifiers = quals)
        record.features.append(featureWO)
        name = name.strip("['").strip("']")
        outputFile1.write("{}\t{}\t{}\n".format(name, local.start, newGene.start))
    elif feature.type == "CDS":
        featureWO = SeqFeature(newGene, type = "CDS", qualifiers = quals)
        record.features.append(featureWO)
    else:
        record.features.append(feature)
SeqIO.write(record, outputFile2, "genbank")
outputFile1.close()
outputFile2.close()
