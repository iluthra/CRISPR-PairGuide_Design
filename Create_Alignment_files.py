#This program just creates alignment files for guides in the input region
import os
import sys
import argparse
import tables
import string
import numpy as np
import cProfile
import fileinput
dna_comp = None

def comp(seq_str):
    """complements the provided DNA sequence and returns it"""
    global dna_comp

    if dna_comp is None:
        dna_comp = string.maketrans("ATCGMRWSYKNatcgmrwsykn",
                                    "TAGCKYWSRMNtagckywsrmn")
    return seq_str.translate(dna_comp)


def revcomp(seq_str):
    """returns reverse complement of provided DNA sequence"""
    return comp(seq_str)[::-1]
    

def str_from_nparray(vals):
    """converts a numpy array into a sequence string"""
    return "".join(chr(x) for x in vals)


def get_seq(chrom,start,end,strand,genome_file):
     

    if chrom not in genome_file.root:
        known_chrom = [node.name for node in genome_file.root]
        raise ValueError("unknown chromosome %s, possible chromosomes are: "
                         + ",".join(known_chrom))

    
    # numpy array representation of sequence and convert to string
    chrom_node = genome_file.getNode("/%s" %chrom)
    np_seq = chrom_node[start-1:end]
    seq_str = str_from_nparray(np_seq)

    if strand == "-":
        # reverse-complement sequence
        seq_str = revcomp(seq_str)


    return seq_str

#searches input region for gRNAs
def finding_guides(seq, inStart, inEnd,gRNA_length):
    guideRna = []
    #has list with RNAs with proper direction
    guideRnaOutput = []
    location = []
    direction = []
    i = 0
    while i < (len(seq)-1):
        
        if seq[i] == "G":
            
            if seq[i+1] == "G":
               
                #only if we are at least 23 bases in can we store a guide RNA, so check and store in list
                if i >= gRNA_length+1:
                    #gRNA must start with a G
                    if seq[i-(gRNA_length+1)] == "G":
                        guideRna.append(seq[i-(gRNA_length+1):i-1])
                        guideRnaOutput.append(seq[i-(gRNA_length+1):i-1])
                        #guideRnaPam.append(seq[i-21:i-1] + " " + seq[i-1:i+2])
                        start = (i - (gRNA_length+1)) + inStart
                        end = (i - 1) + inStart
                        location.append(str(start) + '-' + str(end))
                        direction.append('fwd')
        #increment counter
        i = i + 1


    #check reverse strand for guides
    w = 0
    forward = ""
    while w < (len(seq)-1):
        
        if seq[w] == "C":
            
            if seq[w+1] == "C":
                #only if we are at least 23 bases in can we store a guide RNA, so check and store in list
                if (w+(gRNA_length+3)) <= (len(seq)):
                   
                    forward = seq[w+3: w+(gRNA_length+3)]
                    CODE={'A':'T','T':'A','C':'G','G':'C','N':'N'} 
                    minus_seq=''
                    for c in forward:
                        minus_seq=minus_seq+CODE[c]
                        reverse_seq=minus_seq[::-1]
                        #make sure the gRNA starts with a G
                    if reverse_seq[:1] == "G":
                        guideRnaOutput.append(reverse_seq)
                        guideRna.append(forward)
                        #guideRnaPam.append(seq[w:w+3] + " " + seq[w+3:w+23])
                        start = w+3+inStart
                        end = w+(gRNA_length+3)+inStart
                        location.append(str(start) + '-' + str(end))
                        direction.append('rev')
        #increment counter
        w = w + 1


    #print guideRna
    #print guideRnaPam
    #print guideRnaOutput
    #print location
    #print direction
    b = 0
    #remove any guides that have N's in them
    while b < len(guideRna):
        if "N" in guideRna[b]:
            del guideRna[b]
            del guideRnaOutput[b]
            del location[b]
            del direction[b]
            b = b + 1
        else:
            b = b +1 
    
    #print len(guideRna)
    return guideRna, guideRnaOutput, location, direction

def main():

    #read in input region from command line  
    if len(sys.argv) < 2:
      sys.stderr.write("usage: %s <chrom> [<start> <end>]\n" % sys.argv[0])
      exit(2)

    inChrom = sys.argv[1]

    if len(sys.argv) > 2:
      inStart = int(sys.argv[2])
      inEnd = int(sys.argv[3])
    else:
      start = ""
      end = ""

    sys.stderr.write("%s %d %d\n" % (inChrom, inStart, inEnd))
    inStrand = "+"
    #rename depending on input gene
    fastainput = open("Input_GATA3.fasta", 'w')
    #change path to hdf5 file you want to use 
    seq_h5 = tables.openFile("/iblm/netapp/data1/external/GRC37/GRC37.h5", "r")
    
    inputRegion = get_seq(inChrom, inStart, inEnd, inStrand, seq_h5)
    
    fastainput.write(">" + inChrom + '\n')
    fastainput.write(inputRegion)
    fastainput.close()
    #reopen the fasta input file from above
    f = open("Input_GATA3.fasta", "r")
    
    
    file = f.readlines()
    #create output file and rename depending on gene of interest
    outfile = open('GuideRNAs_GATA3','w')
    outfile2 = open('FQGuides_GATA3.fq' , 'w')
  
    #declare empty lists
    sequence = []
    header = []
    
    seq = ""
    header = ""
    GC = []

    #store all the headers in a list
    for f in file:

        if f.startswith('>'):
            header = header + f
            header = header.splitlines()[0]
        #get ride of new line charaters and spaces
        else:
            f = f.replace(" ", "")
            f = f.replace("\n", "")
            seq = seq + f

    #make it all upper case, easier to parse
    seq = seq.upper()
    
    
    #call function to find all 20, 19 and 21 bp guide RNAs that start with G
    gRNA_length = 20
    guideRna1, guideRnaOutput1, location1, direction1 = finding_guides(seq,inStart, inEnd,gRNA_length)
    gRNA_length = 19
    guideRna2,guideRnaOutput2, location2, direction2 = finding_guides(seq,inStart, inEnd,gRNA_length)
    gRNA_length = 21
    guideRna3, guideRnaOutput3, location3, direction3 = finding_guides(seq,inStart, inEnd,gRNA_length)
    #store them all together
    guideRna = guideRna1 + guideRna2 + guideRna3
    guideRnaOutput = guideRnaOutput1 + guideRnaOutput2 + guideRnaOutput3
    location = location1 + location2 + location3
    direction = direction1 + direction2 + direction3
    #print guideRna
    #print guideRnaOutput

    #print len(guideRna)
    #create .fq file used for alignment
    j = 0
    outfile.write(header)
    while j<len(guideRna):

        outfile.write(guideRna[j] + '\n')
        outfile2.write('@' + str(j) + '\n')
        outfile2.write(guideRna[j] + '\n')
        outfile2.write('+' + '\n')
        score = "I"*len(guideRna[j])
        outfile2.write(score + '\n')

    
        j = j+1


    #if we don't close and reopen FQGuides.fq we cannot excute the commands since FQGuides.fq was only a write file
    outfile2.close()
    #reopen the .fq file created above
    open("FQGuides_GATA3.fq", "r")
    print "We created a .fq file with all the guides"

    #call bwa aligner from this code and get information from the sam file
    #change input reference genome and output file name here
    cmd = 'bwa aln -o 0 -n 3 -N -t 5 hg37.fa FQGuides_GATA3.fq > aln_GATA3.sai'
    os.system(cmd)
    print "Alignment Step has finished"

    #change input reference genome and output file name here
    cmd2 = 'bwa samse -n 1000 hg37.fa aln_GATA3.sai FQGuides_GATA3.fq > aln_GATA3.sam'
    os.system(cmd2)
    
    print "Sam file has been created"
    import re
 
    seq_h5.close()


        
    
main()
