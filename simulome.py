###############################################################################################################################################
# Copyright (c) 2016 - Adam Price
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this 
# software and associated documentation files (the "Software"), to deal in the Software without restriction, 
# including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, 
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##############################################################################################################################################

import StringIO
import sys
import os
import csv
import copy
import numpy
import random
from optparse import OptionParser
from optparse import OptionGroup
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

#Defaults for optional parameters.
desired_gene_count = 100 - 1
spacer_length_default = 500
random_spacer_len = 0
feature_type = "gene"
feature_type_len = len(feature_type)
strict_duplicates = 1
isVerbose = 0
run_snp = 0
run_dup = 0
run_del = 0
run_ins = 0
dup_percent = 0.0
operon_percent = 0

#Other globals.
anno_header_offset = 0
footer = []
genome_outfile = ""
gff_outfile = ""
mut_genome_outfile = ""
mut_gff_outfile = ""

##############
# writeGenome:  This function will write a provided annotation file and genome file. 
#               The isMut parameter is a flag to modify the filename for mutated variants of the simulated genome.
##############	
def writeGenome(targetGff, targetGenome, isMut):
    #Determine if this is the mutated genome or the original simulation.
    if isMut == 1:
        gff_outfile_location = mut_gff_outfile 
        genome_outfile_location = mut_genome_outfile 
    else:
        gff_outfile_location = gff_outfile 
        genome_outfile_location = genome_outfile 

    print "\tWriting annotation file: " + gff_outfile_location
    print "\tWriting FASTA file: " + genome_outfile_location

    #Write annotation files.
    try:
        outfile = open(gff_outfile_location,"w")
        for line in targetGff:
            for item in line:
                outfile.write(str(item) + "\t")
            outfile.write("\n")
        outfile.write(footer[0][0] + "\n")
    except Exception,e:
        print "Error writing mutated annotation file."
        print str(e)
        sys.exit()
    finally:
        outfile.close()

    #Write genome files.
    try:
        curGenomeString = ""
        outfile = open(genome_outfile_location,"w")
        genome = curGenomeString.join(targetGenome)
        outfile.write(">" + targetGff[anno_header_offset+2][0] + "\n")
        outfile.write(genome)
    except Exception,e:
        print "Error writing simulated genome FASTA file."
        print str(e)
        sys.exit()
    finally:
        outfile.close()


###############
# getInsertSeq:  This function generates a random sequence that is used to simulate insertion events.
###############	
def getInsertSeq(insLen):
    insertSeq = ""
    for i in range(insLen):
        randomSelection = random.randint(1,4)
        if randomSelection == 1:
            nucleotide = "A"
        if randomSelection == 2:
            nucleotide = "C"
        if randomSelection == 3:
            nucleotide = "G"
        if randomSelection == 4:
            nucleotide = "T"
        
        insertSeq = insertSeq + nucleotide
    return insertSeq


#############
# mutateBase:  Given a nucleotide base, this function will return a randomly selected mutation base.
#############	
def mutateBase(base):
    randomSelection = random.randint(1,4)
    if randomSelection == 1:
        nucleotide = "A"
    if randomSelection == 2:
        nucleotide = "C"
    if randomSelection == 3:
        nucleotide = "G"
    if randomSelection == 4:
        nucleotide = "T"
    
    #Don't change a base to itself.
    if nucleotide == base:
        return mutateBase(base)  #Fancy recursive loop.
    else:
        return nucleotide


###############
# randomSpacer:  This function simulates random sequences to be used as intergenic regions. 
###############	
def randomSpacer(spacerLen, curGenome):
    spacer = ""
    curGenomeString = ""
    for i in range(spacerLen):
        randomSelection = random.randint(1,4)
        if randomSelection == 1:
            nucleotide = "A"
        if randomSelection == 2:
            nucleotide = "C"
        if randomSelection == 3:
            nucleotide = "G"
        if randomSelection == 4:
            nucleotide = "T"
        
        spacer = spacer + nucleotide

    if len(simulated_genome) == 0:
        if isVerbose == 2:
            print "Appending opening intergenic region..." 
        curGenome.append(spacer)
        return curGenome
    
    #BLAST the sequence we just created so we don't accidently create unwanted duplicate regions.
    seq1 = SeqRecord(Seq(curGenomeString.join(curGenome)), id = "simulatedGenome")
    seq2 = SeqRecord(Seq(spacer), id = "seq2") 
    if isVerbose == 2:
        print "Simulating intergenic region.\n\tCurrent position at: " + str(len(curGenomeString.join(curGenome)))
    
    seq1_file = output_dir + "/" + "curGenomeSim.fasta"
    seq2_file = output_dir + "/" + "queryGeneSim.fasta"
    SeqIO.write(seq1, seq1_file, "fasta")
    SeqIO.write(seq2, seq2_file, "fasta")

    blastNucOutput = NcbiblastnCommandline(query=seq1_file, subject=seq2_file, outfmt=5)()[0]
    blast_nuc_record = NCBIXML.read(StringIO.StringIO(blastNucOutput))

    alignments = len(blast_nuc_record.alignments)
    if alignments == 0:
        if isVerbose == 2:
            print "Appending intergenic region."
        curGenome.append(str(spacer))
    else:
        if isVerbose == 2:
            print "Intergenic region generated aligns to position in simulated genome.  Re-attempting."
            for alignment in blast_nuc_record.alignments:
                for hsp in alignment.hsps:
                    print "\tIdentity: " + str(hsp.identities)
                    print "\tScore: " + str(hsp.score)
                    print "\tAlignment Length: " + str(hsp.align_length)
                    print "\tMismatches: " + str(hsp.align_length - hsp.identities)
                    print "\n"
        return randomSpacer(spacerLen, curGenome)

    return curGenome


#################
# simulateDelete:  This function simulates deletion events.  It will remove sequence from each gene 
#                  of the specified length until the number of desired deletion events is reached.  
#################	
def simulateDelete(curGff, simGenome, delLen, numDels):
    delGenomeString = ""
    delGff = copy.copy(curGff)
    delGffOrig = copy.copy(curGff)
    delGenome = copy.copy(simGenome)   
    delSet = []
    newDelGenome = []                  
    curStart = 0
    curEnd   = 0
    curRemovedLength = 0
    totalRemovedLength = 0
    expectedDelLength = delLen * numDels
	
    #Iterate the annotation file, the extra anno_header_offset is to account for header lines we want to ignore.
    for i in range((len(delGff)-anno_header_offset)):
        #Intergenic Region
        if isVerbose == 2:
            print "=====INTERGENIC REGION====="
        if i == 0:
            curStart = 0
            curEnd = delGffOrig[i+anno_header_offset][3]
        else:
            curStart = delGffOrig[i+anno_header_offset-1][4]  #anno_header_offset is the header cancellation, -1 more is the previously examined gene.	
            curEnd = delGffOrig[i+anno_header_offset][3]
        newDelGenome.append(delGenomeString.join(delGenome)[curStart:curEnd])
        if isVerbose == 2:
            print "Adding " + str(len(delGenomeString.join(delGenome)[curStart:curEnd])) + " bases to genome.  Total length is now " + str(len(delGenomeString.join(newDelGenome)))
            print "--------------"
            print "=====GENE====="
	
        curStart = delGff[i+anno_header_offset][3]
        curEnd = delGff[i+anno_header_offset][4]
        newStart = 0
        newEnd = 0
		
        origGene = delGenomeString.join(delGenome)[curStart:curEnd]
        newSeq = copy.copy(origGene)
        origGeneLen = len(origGene)
        newGeneLen = origGeneLen
        origGeneAnno = delGff[i+anno_header_offset]
        newGeneAnno = copy.copy(origGeneAnno)
		
        if origGeneLen <= expectedDelLength:
            print "Warning: " + origGeneAnno[0] + " length is smaller than expected deletion.  Omitting this gene."
            continue
			
        #Create deletions.
        for j in range(numDels):
            delPos = random.randint(0, origGeneLen-curRemovedLength-(delLen + 1))
            newSeq = newSeq[:delPos] + newSeq[delPos+delLen:]
            curRemovedLength = curRemovedLength + delLen

        #Update annotation information to reflect changes in the sequence.
        newStart = curStart - totalRemovedLength
        newEnd = curEnd - totalRemovedLength - curRemovedLength
        totalRemovedLength = totalRemovedLength + curRemovedLength
        newGeneAnno[3] = newStart
        newGeneAnno[4] = newEnd
        curRemovedLength = 0
        delGff[i+anno_header_offset] = newGeneAnno
        newDelGenome.append(newSeq)

        if isVerbose == 2:
            print "--------------"
		
        #Handle the final spacer at the end.
        if i == (len(delGff)-(anno_header_offset + 1)):
            if isVerbose == 2:
                print "=====INTERGENIC REGION====="
            curStart = delGffOrig[i+anno_header_offset][4]
            curEnd = len(delGenomeString.join(delGenome))  #Genome ends with spacer, so genome length is appropriate for the final endpoint.
            newDelGenome.append(delGenomeString.join(delGenome)[curStart:curEnd])
            if isVerbose == 2:
                print "Adding " + str(len(delGenomeString.join(delGenome)[curStart:curEnd])) + " bases to genome.  Total length is now " + str(len(delGenomeString.join(newDelGenome)))
                print "--------------"

    #Sanity check that genome shrunk the expected size.
    if (len(delGenomeString.join(delGenome)) - totalRemovedLength) != (len(delGenomeString.join(newDelGenome))):
        print "Error: Unexpected genome length after deletion events. Terminating."
        print str((len(delGenomeString.join(delGenome)) - totalRemovedLength))
        print str((len(delGenomeString.join(newDelGenome))))
        sys.exit()
    else:
        if isVerbose >= 1:
            print "Deletions simulated."

    return [delGff, newDelGenome]


#################
# simulateInsert:  This function simulates insertion events.  It will insert random sequences into genes 
#                  of the specified length until the number of desired insertion events is reached.  
#################	
def simulateInsert(curGff, simGenome, insertLen, numInsertions, isCopyEvent):
    insGenomeString = ""
    insGff = copy.copy(curGff)
    insGffOrig = copy.copy(curGff)
    insGenome = copy.copy(simGenome)   
    insSet = []
    newInsGenome = []                  
    curStart = 0
    curEnd   = 0
    curAddedLength = 0
    totalAddedLength = 0

    #Iterate the annotation file, the extra 6 is to account for header lines we want to ignore.
    for i in range((len(insGff)-anno_header_offset)):
        #Intergenic Region
        if isVerbose == 2:
            print "=====INTERGENIC REGION====="

        if i == 0:
            curStart = 0
            curEnd = insGffOrig[i+anno_header_offset][3]
        else:
            curStart = insGffOrig[i+anno_header_offset-1][4]  #i-6 is the header cancellation, -1 more is the previously examined gene.	
            curEnd = insGffOrig[i+anno_header_offset][3]
        newInsGenome.append(insGenomeString.join(insGenome)[curStart:curEnd])
        
        if isVerbose == 2:
            print "Adding " + str(len(insGenomeString.join(insGenome)[curStart:curEnd])) + " bases to genome.  Total length is now " + str(len(insGenomeString.join(newInsGenome)))
            print "--------------"
            print "=====GENE====="

        curStart = insGff[i+anno_header_offset][3]
        curEnd = insGff[i+anno_header_offset][4]
        newStart = 0
        newEnd = 0
		
        origGene = insGenomeString.join(insGenome)[curStart:curEnd]
        newSeq = copy.copy(origGene)
        origGeneLen = len(origGene)
        newGeneLen = origGeneLen
        origGeneAnno = insGff[i+anno_header_offset]
        newGeneAnno = copy.copy(origGeneAnno)
		
        #Create insertions.
        for j in range(numInsertions):
            insertPos = random.randint(0, origGeneLen-1)
            if isCopyEvent == 1:
                if isVerbose == 2:
                    print "Copying sequence of length " + str(insertLen) + " randomly from simulated genome."
					
                copyStart = random.randint(1, len(insGenomeString.join(simGenome))-(insertLen+1))
                insertionSeq = insGenomeString.join(insGenome)[copyStart:copyStart+insertLen]
            else:
                insertionSeq = getInsertSeq(insertLen)
            newSeq = newSeq[:insertPos] + insertionSeq + newSeq[insertPos:]
            curAddedLength = curAddedLength + insertLen
		
        #Update annotation information to reflect changes in the sequence.
        newStart = curStart + totalAddedLength
        newEnd = curEnd + totalAddedLength + curAddedLength
        totalAddedLength = totalAddedLength + curAddedLength
        newGeneAnno[3] = newStart
        newGeneAnno[4] = newEnd
        curAddedLength = 0
        insGff[i+anno_header_offset] = newGeneAnno
        newInsGenome.append(newSeq)

        if isVerbose == 2:
            print "--------------"
		
        #Handle the final spacer at the end.
        if i == (len(insGff)-(anno_header_offset+1)):
            if isVerbose == 2:
                print "=====INTERGENIC REGION====="
            curStart = insGffOrig[i+anno_header_offset][4]
            curEnd = len(insGenomeString.join(insGenome))  #Genome ends with spacer, so genome length is appropriate for the final endpoint.
            newInsGenome.append(insGenomeString.join(insGenome)[curStart:curEnd])
            if isVerbose == 2:
                print "Adding " + str(len(insGenomeString.join(insGenome)[curStart:curEnd])) + " bases to genome.  Total length is now " + str(len(insGenomeString.join(newInsGenome)))
                print "--------------"

    #Sanity check that genome grew the expected size.
    if (len(insGenomeString.join(insGenome)) + totalAddedLength) != (len(insGenomeString.join(newInsGenome))):
        print "Error: Unexpected genome length after insertion. Terminating."
        print str((len(insGenomeString.join(insGenome)) + totalAddedLength))
        print str((len(insGenomeString.join(newInsGenome))))
        print str(totalAddedLength)
        sys.exit()
    else:
        if isVerbose >= 1:
            print "Insertion simulated."

    return [insGff, newInsGenome]

	
##############
# simulateSNP:  This function simulates SNPs.  A window is selected and random positions in the window are mutated until the 
#               desired number of SNPs is properly simulated.  If the window size exceeds the gene size, window size is set equal
#               to the size of the gene.
##############	
def simulateSNP(curGff, simGenome, windowLen, numSNP):
    snpGenomeString = ""
    snpGff = copy.copy(curGff)
    snpGenome = copy.copy(simGenome)   #This is a copy of the original genome that is safe to work with. Because pointers.
    snpSet = []
    newSNPGenome = []                  #This is a modified version of the original genome containing SNPs.
    curStart = 0
    curEnd   = 0

    #Iterate the annotation file, the extra 6 is to account for header lines we want to ignore.
    for i in range((len(snpGff)-anno_header_offset)):	
		#Intergenic Region
        if isVerbose == 2:
            print "=====INTERGENIC REGION====="
        if i == 0:
            curStart = 0
            curEnd = snpGff[i+anno_header_offset][3]
        else:
            curStart = snpGff[i+anno_header_offset-1][4] #i-6 is the header cancellation, -1 more is the previously examined gene.	
            curEnd = snpGff[i+anno_header_offset][3]
        newSNPGenome.append(snpGenomeString.join(snpGenome)[curStart:curEnd])

        if isVerbose == 2:
            print "Adding " + str(len(snpGenomeString.join(snpGenome)[curStart:curEnd])) + " bases to genome.  Total length is now " + str(len(snpGenomeString.join(newSNPGenome)))
            print "--------------"
            print "=====GENE====="

        #Gene
        curStart = snpGff[i+anno_header_offset][3]
        curEnd = snpGff[i+anno_header_offset][4]
		
        origGene = snpGenomeString.join(snpGenome)[curStart:curEnd]
        origGeneLen = len(origGene)
		
        #If the desired window size is larger than the gene, shrink the window to be equal to the size of the gene.
        resized = 0
        origWindowLen = windowLen
        if isVerbose == 2:
            print "WindowLen: " + str(windowLen)
            print "GeneLen: " + str(origGeneLen)
        if origGeneLen <= windowLen:
            if isVerbose >= 1:
                print "Warning: Window length is larger than target gene.  Adjusting window size to gene length."
            resized = 1
            windowLen = origGeneLen

        #Identify a window, establish bit mask to keep track of mutated positions, and extract sequence from that window.
        windowStart = random.randint(0, origGeneLen-windowLen)
        windowEnd = windowStart + windowLen
        if isVerbose == 2:
            print "Window selected: " + str(windowStart) + " : " + str(windowEnd)
        windowBitMask = [0] * windowLen
        windowSeq = origGene[windowStart:windowEnd]
        numMutated = 0
        mutatedSeq = ""

        #Randomly mutate positions in the target window until the desired number of SNPs is achieved.
        while numMutated < numSNP:
            for j in range(len(windowSeq)):
                mutateChance = random.randint(0,windowLen)
                if mutateChance <= numSNP:
                    #Don't mutate the same position twice.
                    if windowBitMask[j] == 1:
                        continue
                    if isVerbose == 2:
                        print "MUTATE " + str(j) + " : " + windowSeq[j]
                    newBase = mutateBase(windowSeq[j])
                    windowSeq = windowSeq[:j] + newBase + windowSeq[j+1:]
                    windowBitMask[j] = 1
                    numMutated = numMutated + 1
                    if numMutated >= numSNP:
                        break
        
        #Reassemble the mutated gene and put it back in the genome.
        mutatedGene = origGene[0:windowStart] + windowSeq + origGene[windowEnd:origGeneLen]
        if origGeneLen != len(mutatedGene):
            print "Error: Sanity check failure in simulateSNP.  Mutated gene length doesn't match original gene length."
            sys.exit()
        else:
            newSNPGenome.append(mutatedGene)
            if isVerbose == 2:
                print "Adding " + str(len(mutatedGene)) + " bases to genome.  Total length is now " + str(len(snpGenomeString.join(newSNPGenome)))

        #If the window size was adjusted previously, set it back to it's original value.
        if resized == 1:
            windowLen = origWindowLen
            resized = 0

        if isVerbose == 2:       
            print "--------------"
		
        #Handle the final spacer at the end.
        if i == (len(snpGff)-(anno_header_offset+1)):
            if isVerbose == 2:
                print "=====INTERGENIC REGION====="
            curStart = snpGff[i+anno_header_offset][4]
            curEnd = len(snpGenomeString.join(snpGenome)) #Genome ends with spacer, so genome length is appropriate for the final endpoint.
            newSNPGenome.append(snpGenomeString.join(snpGenome)[curStart:curEnd])
            if isVerbose == 2:
                print "Adding " + str(len(snpGenomeString.join(snpGenome)[curStart:curEnd])) + " bases to genome.  Total length is now " + str(len(snpGenomeString.join(newSNPGenome)))
                print "--------------"

    #Sanity check on genome size.
    if len(snpGenomeString.join(snpGenome)) != len(snpGenomeString.join(newSNPGenome)):
        print "Error: SNP genome length doesn't match original genome length. Terminating."
        sys.exit()
    else:
        if isVerbose >= 1:
            print "SNP genome simulated."

    return [snpGff, newSNPGenome]


####################
# simulateDuplicate:  This function creates duplicate regions of a size specified by the user.  
####################		
def simulateDuplicate(curGff, simGenome, percentDuplicate):
    dupGenomeString = ""
    genomeLen = len(dupGenomeString.join(simGenome))
    dupGff = copy.copy(curGff)
    dupGenome = copy.copy(simGenome)
    #Be sure 'percentDuplicate' is in decimal format.
    dupLen = percentDuplicate * genomeLen
    dupSet = []
    curDupLen = 0

    if isVerbose == 2:
        print "Duplicating " + str(dupLen) + " positions for a " + str(percentDuplicate*100) + " percent replicate genome." 
    
    #Randomly copy genes from the simulated genome until the right amount of duplication is reached.
    while curDupLen <= dupLen:
        randomGene = copy.copy(random.choice(dupGff))
        if len(randomGene) == 1 or randomGene[2] == "region" or randomGene[2] == "chromosome":
            continue

        randomGeneLen = randomGene[4] - randomGene[3]
        curDupLen = curDupLen + randomGeneLen
        dupSet.append(randomGene)

    #Update the GFF data to reflect which genes are duplicates.
    for i in range(len(dupSet)):
        idString = dupSet[i][8]
        newIdString = idString[:idString.index(";")] + "_copy" + str(int(percentDuplicate*100)) + idString[idString.index(";"):] + "_copy" + str(int(percentDuplicate*100))
        firstSlice = newIdString[:newIdString.index(feature_type)] + feature_type + "=copy" + str(int(percentDuplicate*100)) + "_"
        secondSlice = newIdString[newIdString.index(feature_type)+feature_type_len:]
        newIdString = firstSlice + secondSlice
        dupSet[i][8] = newIdString

    if isVerbose == 2:
        print "Total duplicate length is " + str(curDupLen)

    #Add sequences for the selected duplicate regions to a copy of the simulated genome.
    for i in range(len(dupSet)):
        #Append the sequence to the growing genome for the selected duplicate.
        dupSeq = dupGenomeString.join(dupGenome)[dupSet[i][3]:dupSet[i][4]]
        dupGenome.append(dupSeq)

        #Update the GFF file with the new gene.
        curDup = dupSet[i]
        curDup[3] = genomeLen
        curDup[4] = genomeLen + len(dupSeq)
        dupGff.append(curDup)
        genomeLen = genomeLen + len(dupSeq)

        #Add a spacer.
        if random_spacer_len == 1:
            spacer_length = random.randint(0, 2000)
        else:
            spacer_length = spacer_length_default

        randomSpacer(spacer_length,dupGenome)
        genomeLen = genomeLen + spacer_length

    return [dupGff, dupGenome]


################
# Main function:  This function handles menu controls, flow control, and sampling and creation of the initial simulated genome.
################
if __name__ == '__main__':
    usage = "%prog -f <genome.fasta> -a <genome.gff> -o <destination> <RUN MODE ARGUMENTS> <OPTIONAL ARGUMENTS>"
    parser = OptionParser(usage)

    required = OptionGroup(parser, "Required Arguments")
    required.add_option("-f", "--genome", dest="genome_file", help="File representing genome. FASTA nucleotide format. ")
    required.add_option("-a", "--annotation", dest="anno_file", help="File containing genome annotation information in GTF/GFF3 format. ")
    required.add_option("-o", "--output_prefix", dest="prefix", help="Creates a folder named with the supplied prefix containing output files. ")

    snpmode = OptionGroup(parser, "SNP Run Mode")
    snpmode.add_option("-1", "--snp", dest="snp_mode", help="Mutated genome will contain SNPs. TRUE/FALSE.")
    snpmode.add_option("-s", "--num_snp", dest="num_snp", help="The number of SNPs to simulate in each gene. Required for SNP mode.")
    snpmode.add_option("-w", "--window", dest="snp_window", help="Window size for SNP mutations. Allows for adjustment of SNP density. (OPTIONAL). ")

    indelmode = OptionGroup(parser, "Insertion/Deletion Run Mode")
    indelmode.add_option("-2", "--indel", dest="indel_mode", help="Mutated genome will contain insertion/deletion mutations.  [1 = Insertions only; 2 = Deletions only; 3 = Both insertions and deletions.]")
    indelmode.add_option("-n", "--insert_length", dest="ins_len", help="Length of insertion mutations. Required for insertion mode.")
    indelmode.add_option("-i", "--num_insert", dest="num_ins", help="The number of insertions to simulate in each gene. (DEFAULT = 1)")
    indelmode.add_option("-k", "--isCopyEvent", dest="is_copy_event", help="Insertion sequences will be randomly copied from the other parts of the genome. TRUE/FALSE. (DEFAULT = FALSE)")
    indelmode.add_option("-m", "--delete_length", dest="del_len", help="Length of deletion mutations. Required for deletion mode.")
    indelmode.add_option("-d", "--num_delete", dest="num_del", help="The number of deletions to simulate in each gene. (DEFAULT = 1)")

    dupmode = OptionGroup(parser, "Duplication Run Mode")
    dupmode.add_option("-3", "--duplicate", dest="dup_mode", help="Mutated genome will contain duplicate regions. TRUE/FALSE. ")
    dupmode.add_option("-c", "--percent", dest="percent_dup", help="Percent of duplicate regions in the genome. Required for duplication mode.")

    optionals = OptionGroup(parser, "Optional Arguments")
    optionals.add_option("-g", "--num_genes", dest="num_genes", help="Number of genes to simulate. (DEFAULT=100).")
    optionals.add_option("-l", "--intergenic_len", dest="spacer_len", help="Length of intergenic regions. For random length intergenic regions, specify 0 for this option. (DEFAULT=500).")
    optionals.add_option("-p", "--operon_level", dest="operon_level", help="Simulate operons.  Input should be approximate percentage of desired operon content. (DEFAULT=0)")
    optionals.add_option("-x", "--seed", dest="seed", help="Specifies a seed for random number generator.  (DEFAULT=RANDOM).")
    optionals.add_option("-t", "--type", dest="type", help="Feature type to simulate from annotation file. I.E: gene, exon, CDS. Case sensitive.(DEFAULT=gene).")
    optionals.add_option("-y", "--strict_dup", dest="strict_duplicates", help="Allow duplicate sequence regions to exist in the initial genome simulation.  TRUE/FALSE.  (DEFAULT=TRUE).")
    optionals.add_option("-v", "--verbose", dest="isVerbose", help="Verbose level. [0 = Quiet, 1 = Verbose, 2 = Very Verbose] (DEFAULT=1).")

    parser.add_option_group(required)
    parser.add_option_group(snpmode)
    parser.add_option_group(indelmode)
    parser.add_option_group(dupmode)
    parser.add_option_group(optionals)
    (options, args) = parser.parse_args()
	
    if len(args) != 0 or not options.genome_file or not options.anno_file or not options.prefix:
        parser.error("Required arguments have not been supplied.  \n\t\tUse -h to get more information.")
        sys.exit()

    #SNP mode user controls.
    if options.snp_mode:
        if options.snp_mode.upper() == "TRUE" or options.snp_mode.upper() == "FALSE":
            if options.snp_mode.upper() == "TRUE":
                run_snp = 1
            else:
                run_snp = 0
        else:
            print "SNP mode must be TRUE or FALSE."
            sys.exit()
        if not options.num_snp:
            print "Number of SNPs per gene is required for SNP run mode. (-s)"
            sys.exit()
        else:
            if not options.num_snp.isdigit() or int(options.num_snp) <= 0:
                print "Number of SNPs must be an integer greater than 0."
                sys.exit()
        if options.snp_window:
            if not options.snp_window.isdigit() or int(options.snp_window) <= 1:
                print "Invalid selection for SNP window size.  Please enter an integer greater than 1."
                sys.exit()
            if int(options.num_snp) >= int(options.snp_window):
                print "Number of SNPs must be smaller than SNP window size."
                sys.exit()
    
    #Indel mode user controls.	
    if options.indel_mode:
        if int(options.indel_mode) < 1 or int(options.indel_mode) > 3:
            print int(options.indel_mode)
            print "Invalid selection for Insertion/Deletion run mode. \n[1 = Insertions only; 2 = Deletions only; 3 = Both insertions and deletions.]"
            sys.exit()
        else:
            run_indel = int(options.indel_mode)
            if run_indel == 1 or run_indel == 3:
                if not options.ins_len:
                    print "Insert length is required for Insert run mode. (-n)"
                    sys.exit()
                else:
                    if not options.ins_len.isdigit() or int(options.ins_len) <= 0:
                        print "Insert length must be a positive integer."
                        sys.exit()
                if options.num_ins:
                    if not options.num_ins.isdigit() or int(options.num_ins) <= 0:
                        print "Invalid number of insertions.  Please enter a positive integer."
                        sys.exit()
                else:
                    options.num_ins = 1 
                if options.is_copy_event:
                    if options.is_copy_event.upper() == "TRUE" or options.is_copy_event.upper() == "FALSE":
                        if options.is_copy_event.upper() == "TRUE":
                            options.is_copy_event = 1
                        else:
                            options.is_copy_event = 0
                else:
                    options.is_copy_event = 0
                run_ins = 1
            if run_indel == 2 or run_indel == 3:
                if not options.del_len:
                    print "Deletion length is required for Delete run mode. (-m)"
                    sys.exit()
                else:
                    if not options.del_len.isdigit() or int(options.del_len) <= 0:
                        print "Deletion length must be a positive integer."
                        sys.exit()
                if options.num_del:
                    if not options.num_del.isdigit() or int(options.num_del) <= 0:
                        print "Invalid number of deletions.  Please enter a positive integer."
                        sys.exit()
                else:
                    options.num_del = 1
                run_del = 1
            if run_indel == 1 and (options.del_len or options.num_del):
                print "Deletion parameters invalid for insertion only mode."
                sys.exit()
            if run_indel == 2 and (options.ins_len or options.num_ins):
                print "Insertion parameters invalid for deletion only mode."
                sys.exit()            

    #Duplicate mode user controls.
    if options.dup_mode:
        if options.dup_mode.upper() == "TRUE" or options.dup_mode.upper() == "FALSE":
            if options.dup_mode.upper() == "TRUE":
                run_dup = 1
            else:
                run_dup = 0
        else:
            print "Duplication mode must be TRUE or FALSE."
            sys.exit()
        if not options.percent_dup:
            print "Duplication percentage required for duplication run mode. (-c)"
            sys.exit()
        if not options.percent_dup.isdigit():
            print "Duplication percent must be an integer between 1-100"
            sys.exit()
        if int(options.percent_dup) < 1 or int(options.percent_dup) > 100:
            print "Duplication percent must be an integer between 1-100"
            sys.exit()
        dup_percent = float(options.percent_dup) * 100.0

    #Operon option controls.
    if options.operon_level:
        if not options.operon_level.isdigit():
            print "Operon percent must be an integer between 1-100"
            sys.exit()
        if int(options.operon_level) < 1 or int(options.operon_level) > 100:
            print "Operon percent must be an integer between 1-100"
            sys.exit()        
        operon_percent = int(options.operon_level)

    #Desired gene count user controls.
    if options.num_genes:
        if options.num_genes.isdigit() and int(options.num_genes) > 1:
            desired_gene_count = int(options.num_genes) + 1
        else:
            print "Error: Number of genes must be a positive number greater than 1."
            sys.exit()

    #Intergenic region size user controls.
    if options.spacer_len:
        if options.spacer_len.isdigit() and int(options.spacer_len) >= 0:
            spacer_length_default = int(options.spacer_len)
            if options.spacer_len == 0:
                random_spacer_len = 1
        else:
            print "Error: Intergenic length must be a positive number."
            sys.exit()

    #Seed controls.
    if options.seed:
        random.seed(options.seed)

    #Annotation type controls.
    if options.type:
        feature_type = options.type
        feature_type_len = len(feature_type)

    #Strict duplicate controls.
    if options.strict_duplicates:
        if options.strict_duplicates.upper() == "TRUE" or options.strict_duplicates.upper() == "FALSE":
            if options.strict_duplicates.upper() == "TRUE":
                strict_duplicates = 1
            else:
                strict_duplicates = 0
        else:
            print "Error: Option, --strict_duplicates must be either TRUE or FALSE."
            sys.exit()

    #Verbose controls.
    if options.isVerbose:
        if options.isVerbose.isdigit():
            if int(options.isVerbose) >= 0 and int(options.isVerbose) <= 2:
                isVerbose = int(options.isVerbose)
            else:
                print "Error: Option, --isVerbose must be one of the following: [0 = Quiet, 1 = Verbose, 2 = Very Verbose]. "
                sys.exit()
        else:
            print "Error: Option, --isVerbose must be one of the following: [0 = Quiet, 1 = Verbose, 2 = Very Verbose]. "
            sys.exit()
			
    #Output file locations.
    output_dir = options.prefix
    genome_outfile = output_dir + "/" + options.prefix + "_simulated.fasta"
    gff_outfile = output_dir + "/" + options.prefix + "_simulated.gff"
    mut_genome_outfile = output_dir + "/" + options.prefix + "_mutated_simulation.fasta"
    mut_gff_outfile = output_dir + "/" + options.prefix + "_mutated_simulation.gff"

    #Create a directory for our output files if it doesn't already exist.
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    #Variables
    annotations = []
    anno_subset = []
    anno_header = []
    gene_len = []
    gene_len_data = {}
    simulated_genome = []
    simulated_gff = []
    selected = []

    #Prepare input files.
    try:
        if isVerbose >= 1:
            print "Reading FASTA files..."
        features = list(SeqIO.parse(options.genome_file, "fasta"))
    except Exception,e:
        print "Error reading genome FASTA file. Terminating execution."
        print str(e)
        sys.exit()

    try:
        gotHeader = 0
        if isVerbose >= 1:
            print "Reading annotation files..."
        with open(options.anno_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            annotations = list(reader)
        for i in range(len(annotations)):
            if (len(annotations[i]) == 1) and gotHeader == 0:
                anno_header.append(annotations[i])
                anno_header_offset = anno_header_offset + 1
                continue
            if (len(annotations[i]) == 1) and gotHeader == 1:
                footer.append(annotations[i])
                break
            if (annotations[i][2] == "region" or annotations[i][2] == "chromosome") and gotHeader == 0:
                anno_header.append(annotations[i])
                anno_header_offset = anno_header_offset + 1
                gotHeader = 1
                continue
            if annotations[i][2] == feature_type:
                anno_subset.append(annotations[i])
                gene_len_data[int(annotations[i][3])] = (int(annotations[i][4]) - int(annotations[i][3]))
                gene_len.append((int(annotations[i][4]) - int(annotations[i][3])))
    except Exception,e:
        print "Error reading genome annotation file. Terminating execution."
        print str(e)
        sys.exit()

    #The maximum simulation genome size can only be the size of the input genome.
    if len(anno_subset) <= desired_gene_count:
        print "Error: Desired gene count for simulated genome exceeds is actual gene count." 
        print "\t" + str(len(anno_subset)) + " features of type '" + feature_type + "' exist in input genome. \nRequested " + str(desired_gene_count) + " genes. "
        sys.exit()

    if len(gene_len) == 0:
        print "Error: No sequences found for feature type: " + feature_type
        sys.exit()

    #Select the proper number of genes in a normal distribution around their lengths. 
    deviation1_start = numpy.mean(gene_len, axis=0) - (numpy.std(gene_len, axis=0)/2)
    deviation1_end = numpy.mean(gene_len, axis=0) + (numpy.std(gene_len, axis=0)/2)
    deviation2_left_start = (numpy.mean(gene_len, axis=0) - (numpy.std(gene_len, axis=0)/2)) - (numpy.std(gene_len, axis=0)/2)
    deviation2_left_end = (numpy.mean(gene_len, axis=0) - (numpy.std(gene_len, axis=0)/2)) - 1
    deviation2_right_start = (numpy.mean(gene_len, axis=0) + (numpy.std(gene_len, axis=0)/2)) + 1
    deviation2_right_end = (numpy.mean(gene_len, axis=0) + (numpy.std(gene_len, axis=0)/2)) + (numpy.std(gene_len, axis=0)/2)
    deviation3_left_start = 1
    deviation3_left_end = (numpy.mean(gene_len, axis=0) - (numpy.std(gene_len, axis=0)/2)) - (numpy.std(gene_len, axis=0)/2) - 1
    deviation3_right_start = (numpy.mean(gene_len, axis=0) + (numpy.std(gene_len, axis=0)/2)) + (numpy.std(gene_len, axis=0)/2) - 1
    deviation3_right_end = max(gene_len) - 1

    #Make a list of all the genes falling within the first, second, and third standard deviation.  
    deviation1_potentials = []
    deviation2_potentials = []
    deviation3_potentials = []
    for startPos, geneLen in gene_len_data.iteritems():
        if geneLen > deviation1_start and geneLen < deviation1_end:
            deviation1_potentials.append(startPos)
            continue
        if (geneLen > deviation2_left_start and geneLen < deviation2_left_end) or (geneLen > deviation2_right_start and geneLen < deviation2_right_end):
            deviation2_potentials.append(startPos)
            continue
        if (geneLen > deviation3_left_start and geneLen < deviation3_left_end) or (geneLen > deviation3_right_start and geneLen < deviation3_right_end):
            deviation3_potentials.append(startPos)

    #Determine how many genes we want to get from each deviation to create a normal distribution of lengths.  
    deviation1_genes_count = int(.683 * desired_gene_count)
    deviation2_genes_count = int(.271 * desired_gene_count)
    deviation3_genes_count = int(.046 * desired_gene_count)
    remainder = desired_gene_count - (deviation1_genes_count + deviation2_genes_count + deviation3_genes_count)
    deviation1_genes_count = deviation1_genes_count + remainder

    #Report about genome structure and selection of genes from a normal distribution.
    if isVerbose >= 1:
        print "==========================================================="
        print "Composition: " + str(len(gene_len)) + " total genes."
        print "\tMean: " + str(numpy.mean(gene_len, axis=0))
        print "\tStandard Dev: " + str(numpy.std(gene_len, axis=0))
        print "\t\tFirst standard deviation range: " + str(deviation1_start) + " - " + str(deviation1_end)
        print "\t\t\t" + str(len(deviation1_potentials)) + " genes identified in this range.  " + str((len(deviation1_potentials)*100)/len(gene_len)) + "% of genes."
        print "\t\t\tSelecting " + str(deviation1_genes_count) + " genes from this range." 
        print "\t\tSecond standard deviation range (left): " + str(deviation2_left_start) + " - " + str(deviation2_left_end)
        print "\t\tSecond standard deviation range (right): " + str(deviation2_right_start) + " - " + str(deviation2_right_end)
        print "\t\t\t" + str(len(deviation2_potentials)) + " genes identified in this range.  " + str((len(deviation2_potentials)*100)/len(gene_len)) + "% of genes."
        print "\t\t\tSelecting " + str(deviation2_genes_count) + " genes from this range." 
        print "\t\tThird standard deviation range (left): " + str(deviation3_left_start) + " - " + str(deviation3_left_end)
        print "\t\tThird standard deviation range (right): " + str(deviation3_right_start) + " - " + str(deviation3_right_end)
        print "\t\t\t" + str(len(deviation3_potentials)) + " genes identified in this range.  " + str((len(deviation3_potentials)*100)/len(gene_len)) + "% of genes."
        print "\t\t\tSelecting " + str(deviation3_genes_count) + " genes from this range." 
        print "\t\t" + str(deviation1_genes_count + deviation2_genes_count + deviation3_genes_count) + " genes total to be selected."
        print "==========================================================="

    #Pick which genes we want to use randomly, but in a normal distribution based on gene length.
    if isVerbose >= 1:
        print "Picking random genes..."
    selectStarts = []
    for i in range(deviation1_genes_count):
        selectStarts.append(random.choice(deviation1_potentials))
    for i in range(deviation2_genes_count):
        selectStarts.append(random.choice(deviation2_potentials))
    for i in range(deviation3_genes_count):
        selectStarts.append(random.choice(deviation3_potentials))
    if isVerbose >= 1:
        print str(len(selectStarts)) + " genes to be selected." 

    for i in range(len(anno_subset)):
        if int(anno_subset[i][3]) in selectStarts:
            selected.append(anno_subset[i])
    if isVerbose >= 1:
        print "Randomly selected " + str(len(selected)) + " genes."      


    #Get the sequences for the selected genes from the FASTA file.
    for i in range(len(selected)):
        startPos = int(selected[i][3])
        endPos = int(selected[i][4])
        selected[i] = [selected[i],features[0].seq[startPos:endPos]]
    if isVerbose >= 1:
        print "Sequences for selected genes acquired..."


    #Create the simulated genome. 
    curPos = 0
    gene_total_len = 0
    simulated_gff = anno_header
    for i in range(len(selected)):
        curGenomeString = ""
        if isVerbose == 2:
            print "========================Processing " + str(i) + "============================"
            print "Current simulated genome length is: " + str(len(curGenomeString.join(simulated_genome)))

        #Add a spacer at the start of the simulated genome.
        if len(simulated_genome) == 0:
            if random_spacer_len == 1:
                spacer_length = random.randint(0, 2000)
            else:
                spacer_length = spacer_length_default
            simulated_genome = randomSpacer(spacer_length,simulated_genome)

            if isVerbose == 2:
                print "Adding beginning spacer."
                print "Current position: " + str(curPos)

            curPos = curPos + spacer_length
            continue

        seq1 = SeqRecord(Seq(curGenomeString.join(simulated_genome)), id = "simulatedGenome")
        seq2 = SeqRecord(selected[i][1], id = "seq2") 

        if isVerbose == 2:
            print "Simulation at " + str(len(seq1.seq)) + " nucleotides and growing..."
            print "GFF at position: " + str(curPos)

        seq1_file = output_dir + "/" + "curGenomeSim.fasta"
        seq2_file = output_dir + "/" + "queryGeneSim.fasta"
        SeqIO.write(seq1, seq1_file, "fasta")
        SeqIO.write(seq2, seq2_file, "fasta")

        blastNucOutput = NcbiblastnCommandline(query=seq1_file, subject=seq2_file, outfmt=5)()[0]
        blast_nuc_record = NCBIXML.read(StringIO.StringIO(blastNucOutput))
		
        #Check if the user is okay with duplicated regions, and grow the simulation accordingly.
        alignments = len(blast_nuc_record.alignments)
        if alignments == 0:
            simulated_genome.append(str(selected[i][1]))
            curGeneLen = len(selected[i][1])
            gene_total_len = gene_total_len + curGeneLen
            curGff = selected[i][0]
            curGff[3] = curPos
            curGff[4] = curPos + curGeneLen
            simulated_gff.append(curGff)
            curPos = curPos + curGeneLen

            if isVerbose == 2:
                print "Current gene length: " + str(len(selected[i][1]))
                print "Gene total length, spacers omitted: " + str(gene_total_len)
                print "Current position: " + str(curPos)
                print selected[i][1]
				
        else:
            if isVerbose == 2:
                print "Natural duplicated region identified in initial genome: "
                for alignment in blast_nuc_record.alignments:
                    for hsp in alignment.hsps:
                        print "\tIdentity: " + str(hsp.identities)
                        print "\tScore: " + str(hsp.score)
                        print "\tAlignment Length: " + str(hsp.align_length)
                        print "\tMismatches: " + str(hsp.align_length - hsp.identities)
                        print "\n"
            if strict_duplicates == 1:  
                simulated_genome.append(str(selected[i][1]))
                curGeneLen = len(selected[i][1])
                gene_total_len = gene_total_len + curGeneLen
                curGff = selected[i][0]
                curGff[3] = curPos
                curGff[4] = curPos + curGeneLen
                simulated_gff.append(curGff)
                curPos = curPos + curGeneLen

        #If we want operons, simply don't include intergenic regions between some genes.
        if operon_percent > 0:
            operon_chance = random.randint(0,100)
            if operon_chance <= operon_percent:
                if isVerbose == 2:
                    print "Operon created at position " + str(curPos)
                continue

        #If random intergenic regions is turned on, create an appropriate length.  2000 is the max for now.
        if random_spacer_len == 1:
            spacer_length = random.randint(0, 2000)
        else:
            spacer_length = spacer_length_default

        randomSpacer(spacer_length,simulated_genome)
        curPos = curPos + spacer_length
        if isVerbose == 2:
            print "Added spacer. Current length: " + str(curPos)
            print "Simulation at " + str(len(curGenomeString.join(simulated_genome))) + " nucleotides and growing..."
            print "GFF at position: " + str(curPos)
            print "========================END GENE " + str(i) + "=============================="

        #Clean up our temporary files.
        os.remove(seq1_file)
        os.remove(seq2_file)

    if isVerbose >= 1:
        print "Simulation complete at " + str(len(curGenomeString.join(simulated_genome))) + " nucleotides."
        print "Gene total length: " + str(gene_total_len)

    #Write the initial simulation data.
    writeGenome(simulated_gff, simulated_genome, 0)

    #Perform any mutations requested and write that data as well.	
    doMutation = 0
    dupData = [simulated_gff, simulated_genome]
    if run_del == 1:
        dupData = simulateDelete(dupData[0], dupData[1], int(options.del_len), int(options.num_del))
        doMutation = 1
    if run_ins == 1:
        dupData = simulateInsert(dupData[0], dupData[1], int(options.ins_len), int(options.num_ins), int(options.is_copy_event))
        doMutation = 1
    if run_snp == 1:
        dupData = simulateSNP(dupData[0], dupData[1], int(options.snp_window), int(options.num_snp))
        doMutation = 1
    if run_dup == 1:
        dupData = simulateDuplicate(dupData[0], dupData[1], dup_percent)
        doMutation = 1

    if doMutation == 1:
        writeGenome(dupData[0], dupData[1], 1)

    print "Process Complete."
	
