from sys import argv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import FeatureLocation, SeqFeature
import random

startCodons = ["ATG", "GTG", "TTG"]
revStartCodons = ["CAT", "CAC", "CAA"]

def randShiftForward(seqLoc, genomeSeq):
    upOrDown = random.randint(1,2)
    streamShift = ""
    if upOrDown == 2:
        start = int(seqLoc.start)
        upperRange = start - 150
        searchSeq = str(genomeSeq[upperRange:start])
        newStartCodon = -1
        for codon in startCodons:
            if newStartCodon == -1:
                newStartCodon = searchSeq.find(codon)
        newStart = upperRange + newStartCodon
        streamShift = "up"
    else:
        start = int(seqLoc.start) + 2
        lowerRange = start + 150
        searchSeq = str(genomeSeq[start:lowerRange])
        newStartCodon = -1
        for codon in startCodons:
            if newStartCodon == -1:
                newStartCodon = searchSeq.find(codon)
        newStart = start + newStartCodon
        streamShift = "down"
    newLoc = FeatureLocation(newStart, seqLoc.end, strand=seqLoc.strand)
    return newLoc, streamShift

def randShiftReverse(seqLoc, genomeSeq):
    upOrDown = random.randint(1,2)
    streamShift = ""
    if upOrDown == 2:
        start = int(seqLoc.end) - 2
        upperRange = start - 150
        searchSeq = str(genomeSeq[upperRange:start])
        newStartCodon = -1
        for codon in revStartCodons:
            if newStartCodon == -1:
                newStartCodon = searchSeq.find(codon)
        newStart = upperRange + newStartCodon + 3
        streamShift = "down"
    else:
        start = int(seqLoc.end)
        lowerRange = start + 150
        searchSeq = str(genomeSeq[start:lowerRange])
        newStartCodon = -1
        for codon in revStartCodons:
            if newStartCodon == -1:
                newStartCodon = searchSeq.find(codon)
        newStart = start + newStartCodon + 3
        streamShift = "up"
    newLoc = FeatureLocation(seqLoc.start, newStart, strand=seqLoc.strand)
    return newLoc, streamShift

def shiftDescription(oldLoc, newLoc, streamShift):
    if oldLoc.strand == 1:
        nucShift = abs(oldLoc.start - newLoc.start)
        description = "{} nuclotide shift {}stream ({} => {})".format(nucShift, streamShift, oldLoc.start, newLoc.start)
    elif oldLoc.strand == -1:
        nucShift = abs(oldLoc.end - newLoc.end)
        description = "{} nuclotide shift {}stream ({} => {})".format(nucShift, streamShift, oldLoc.end, newLoc.end + 1)
    return description

outputFile1 = open(argv[2], "w")
outputFile2 = open(argv[3], "w")
genBankFile = argv[1]
genome = SeqIO.read(genBankFile, "genbank")
seq = genome.seq
features = genome.features
record = SeqRecord(seq, id = genome.id, name = genome.name)
for feature in features:
    if feature.type == "CDS":
        local = feature.location
        geneSeq = seq[local.start:local.end]
        if local.strand == 1:
            newGene, streamShift = randShiftForward(local, seq)
        elif local.strand == -1:
            newGene, streamShift = randShiftReverse(local, seq)
        name = str(feature.qualifiers["locus_tag"])
        featureWO = SeqFeature(newGene, type = name)
        record.features.append(featureWO)
        name = name.strip("['").strip("']")
        shiftDescript = shiftDescription(local, newGene, streamShift)
        outputFile1.write("New {}: {}\n".format(name, shiftDescript))
SeqIO.write(record, outputFile2, "genbank")
outputFile1.close()
outputFile2.close()
