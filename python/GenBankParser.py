import argparse
from Bio import SeqIO
from Bio.Seq import Seq, translate
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import FeatureLocation, SeqFeature
import random
import warnings
from Bio import BiopythonWarning
warnings.simplefilter("ignore", BiopythonWarning)

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="input GenBank file")
parser.add_argument("-g", "--genBankOut", help="output GenBank file")
parser.add_argument("-c", "--csvOutput", help="output text file")

args = parser.parse_args()

startCodons = ["ATG", "GTG", "TTG"]
revStartCodons = ["CAT", "CAC", "CAA"]
stopCodons = ["TAA", "TAG", "TGA"]
revStopCodons = ["TTA", "CTA", "TCA"]

def lookUpstream(genomeSeq, start, startCodonList, stopCodonList, codonDirection):
    print("Looking for start codon upstream...")
    for i in xrange(1,50):
        currentStart = start - (3 * i)
        if codonDirection == 3:
            codon = genomeSeq[currentStart:currentStart+codonDirection]
        elif codonDirection == -3:
            codon = genomeSeq[currentStart+codonDirection:currentStart]
        if codon in startCodonList:
            print("Found new start codon!")
            return currentStart
            break
        elif codon in stopCodonList:
            print("Found stop codon")
            return -1
            break
    print("No new start codon found")
    return -1

def lookDownstream(genomeSeq, start, startCodonList, stopCodonList, codonDirection):
    print("Looking for start codon downstream...")
    for i in xrange(1,50):
        currentStart = start + (3 * i)
        if codonDirection == 3:
            codon = genomeSeq[currentStart:currentStart+codonDirection]
        elif codonDirection == -3:
            codon = genomeSeq[currentStart+codonDirection:currentStart]
        if codon in startCodonList:
            print("Found new start codon")
            return currentStart
            break
        elif codon in stopCodonList:
            print("Found stop codon")
            return -1
            break
    print("No new start codon found")
    return -1

def randShift(seqLoc, genomeSeq):
    upOrDown = random.randint(1,2)
    newStartCodon = -1;
    counter = 0;
    while newStartCodon == -1:
        if upOrDown == 2:
            start = int(seqLoc.start)
            end = int(seqLoc.end)
            searchSeq = str(genomeSeq)
            if genomeSeq[start:start+3] in startCodons:
                newStartCodon = lookUpstream(searchSeq, start, startCodons, stopCodons, 3)
            else:
                newStartCodon = lookUpstream(searchSeq, end, revStartCodons, revStopCodons, -3)
        elif upOrDown == 1:
            start = int(seqLoc.start)
            end = int(seqLoc.end)
            searchSeq = str(genomeSeq)
            if genomeSeq[start:start+3] in startCodons:
                newStartCodon = lookDownstream(searchSeq, start, startCodons, stopCodons, 3)
            else:
                newStartCodon = lookDownstream(searchSeq, end, revStartCodons, revStopCodons, -3)
        if newStartCodon == -1:
            if upOrDown == 1:
                upOrDown = 2
            else:
                upOrDown = 1
        counter += 1
        if counter > 2:
            print("***infinite loop***")
            break
    if counter <= 2:
        if seqLoc.strand == 1:
            newLoc = FeatureLocation(newStartCodon, seqLoc.end, strand=seqLoc.strand)
        else:
            newLoc = FeatureLocation(seqLoc.start, newStartCodon, strand=seqLoc.strand)
    else:
        newLoc = FeatureLocation(seqLoc.start, seqLoc.end, strand=seqLoc.strand)
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
        print(feature.qualifiers["locus_tag"])
        local = feature.location
        geneSeq = seq[local.start:local.end]
        newGene = randShift(local, seq)
        name = str(feature.qualifiers["locus_tag"])
        newSeq = seq[newGene.start:newGene.end]
        strSeq = str(seq)
        if strSeq[newGene.start:newGene.start+3] in startCodons:
            newProduct = newSeq.translate()
            name = name.strip("['").strip("']")
            outputFile1.write("{}\t{}\t{}\n".format(name, local.start, newGene.start))
        elif strSeq[newGene.end-3:newGene.end] in revStartCodons:
            tempSeq = newSeq.reverse_complement()
            newProduct = tempSeq.translate()
            name = name.strip("['").strip("']")
            outputFile1.write("{}\t{}\t{}\n".format(name, local.end, newGene.end))
        quals["translation"] = newProduct
        featureWO = SeqFeature(newGene, type = "gene", qualifiers = quals)
        record.features.append(featureWO)
    elif feature.type == "CDS":
        newSeq = seq[newGene.start:newGene.end]
        strSeq = str(seq)
        if strSeq[newGene.start:newGene.start+3] in startCodons:
            newProduct = newSeq.translate()
        elif strSeq[newGene.end-3:newGene.end] in revStartCodons:
            tempSeq = newSeq.reverse_complement()
            newProduct = tempSeq.translate()
        quals["translation"] = newProduct
        featureWO = SeqFeature(newGene, type = "CDS", qualifiers = quals)
        record.features.append(featureWO)
    else:
        record.features.append(feature)
SeqIO.write(record, outputFile2, "genbank")
outputFile1.close()
outputFile2.close()

