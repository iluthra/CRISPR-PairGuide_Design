#Create a Database for all single guides, specificity is calculated
import os
import time
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

#calculates specificity for all guides
def specificity_calc(listIndex, infoList_current,guide_current,seq_h5):
    
    Shitlist = []
    Shitlist = 0
    #infoList1 now has everything from .sam file per guide
    current_guide_length = len(guide_current)
    while listIndex < len(infoList_current):
            #get sequence for alignment of guide
            refSEQ = get_seq(infoList_current[listIndex][0], int(infoList_current[listIndex][2]), int(infoList_current[listIndex][2]) + (current_guide_length-1) , infoList_current[listIndex][1], seq_h5)            
            
            mismatchPos = []
            m = 0

            while m < len(refSEQ):
                if refSEQ[m] != guide_current[m]:
                    mismatchPos.append(m)
                    m = m + 1
                else:
                    m = m +1


            #set W matrix depending on guide length
            #Zhang Lab
            if current_guide_length == 20:        
                W = [0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445,
                    0.508, 0.613, 0.851, 0.731, 0.828, 0.615, 0.804, 0.685, 0.583]
            elif current_guide_length == 21:
                #added a 0
                W = [0, 0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445,
                    0.508, 0.613, 0.851, 0.731, 0.828, 0.615, 0.804, 0.685, 0.583]
            elif current_guide_length == 19:
                #removed a zero
                W = [0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445,
                    0.508, 0.613, 0.851, 0.731, 0.828, 0.615, 0.804, 0.685, 0.583]    

            d = 0
            l = (current_guide_length-1)


            k = 0
            if len(mismatchPos) == 0:

                listIndex = listIndex + 1
                
            if len(mismatchPos) == 3:
                
                d = ((mismatchPos[2] - mismatchPos[1]) + (mismatchPos[1] - mismatchPos[0])) / (len(mismatchPos) - 1 )
                a = (1 - W[mismatchPos[2]])*(1 - W[mismatchPos[1]])*(1 - W[mismatchPos[0]])
                b = (1.0/(((1-d*1.0/l)*4.0)+1))
                c = 1.0/9.0
                SHit = a*b*c
                Shitlist = SHit + Shitlist 
                listIndex = listIndex + 1
                
            if len(mismatchPos) == 2:

                d = ((mismatchPos[1] - mismatchPos[0])) / (len(mismatchPos) - 1 )
                a =  (1 - W[mismatchPos[1]])*(1 - W[mismatchPos[0]])
                b = (1/(((1-d*1.0/l)*4.0)+1))
                c = 1.0/4.0
                SHit = a*b*c
                Shitlist = SHit + Shitlist
                listIndex = listIndex + 1
                
            if len(mismatchPos) == 1:

                d = 0
                SHit = (1 - W[mismatchPos[0]])*(1/((1-d*1.0/l)*4)+1)
                Shitlist = SHit + Shitlist
                listIndex = listIndex + 1
            #this shouldn't ever happen, but will catch it if it does    
            if len(mismatchPos) > 3:
                
                listIndex = listIndex + 1
    #score for current guide          
    Sguide = 100/  (1 + Shitlist)
    return Sguide

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
                    #gRNA must start with G
                    if seq[i-(gRNA_length+1)] == "G":
                        guideRna.append(seq[i-(gRNA_length+1):i-1])
                        
                        guideRnaOutput.append(seq[i-(gRNA_length+1):i-1])
                        #guideRnaPam.append(seq[i-21:i-1] + " " + seq[i-1:i+2])
                        start = (i - (gRNA_length+1)) + inStart
                        end = (i - 1) + inStart -1
                        location.append(str(start) + '-' + str(end))
                        direction.append('fwd')
        #increment counter
        i = i + 1




    #find guides on the reverse DNA strand
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
                    
                    if reverse_seq[:1] == "G":
                        guideRnaOutput.append(reverse_seq)
                        guideRna.append(forward)
                        #guideRnaPam.append(seq[w:w+3] + " " + seq[w+3:w+23])
                        start = w+3+inStart
                        end = w+(gRNA_length+3)+inStart -1
                        location.append(str(start) + '-' + str(end))
                        direction.append('rev')
        #increment counter
        w = w + 1

    #uncomment to see the actual sequences    
    #print guideRna
    #print guideRnaPam
    #print guideRnaOutput
    #print location
    #print direction
    b = 0
    #remove any guidse with N
    while b < len(guideRna):
        if "N" in guideRna[b]:
            del guideRna[b]
            del guideRnaOutput[b]
            del location[b]
            del direction[b]
            b = b + 1
        else:
            b = b +1 
    
    return guideRna, guideRnaOutput, location, direction

def main():
    
    #user should input the same region they used in Create_alignment_files.py
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
    
    file_fasta = f.readlines()
    #file for Database of single guides
    Database_file = open('Database_GATA3.txt', 'w')

    #write header
    Database_file.write('Location' + '\t' + '\t' + '\t' + '\t' + 'gRNA' + '\t' + '\t' + '\t' + 'Direction' + '\t' + '#other perfect aligns' + '\t' +  'Specificity Score' + '\t' + '%GC Content' + '\t' + '#MM1'+ '\t' + '#MM2'+ '\t' + '#MM3' +  '\n')

    #declare empty lists
    sequence = []
    header = []
 
    seq = ""
    header = ""
    GC = []

    #store all the headers in a list
    for f in file_fasta:

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
    j = 0

    #calcuate %GC content for all the guides and store in a list
    while j<len(guideRna):

        totalgc = guideRnaOutput[j].count("G") + guideRnaOutput[j].count("C")

        gccontent = (totalgc / float(len(guideRna[j]))) * 100.00
 
        GC.append(gccontent)
        j = j+1


    import subprocess
    #open the sam file created in Create_alignment_files.py
    file_sam = open("aln_GATA3.sam", "r")
    print "Sam file has been opened"
    import re
    
    infoListFull = []
    alnList = []
    q = 0
   
    i = 0
   
   
    #parse sam file line by line
    #each line contains all alignments for each guide
    for line in file_sam:
        #skip header lines
        if line.startswith('@'):
            continue
    
        sys.stderr.write("We are on guide number %d\n" % q)
        #extract what we need for calculations
        alnList = re.findall('(chr[\dMXYIVLR]+),([+-])(\d+),\d+M,(\d)', line)

        #this give list for the qth GUIDE (guide we are currently on)
        infoList_current = alnList
        #uncomment to see number of alignments for this guide
        #print len(infoList1)

        
        if len(infoList_current) != 0:
            #removes all reference sequences that match perfectly somewhere else == 0 mismatches

            numPerfAlns = 0
            mm1 = 0
            mm2 = 0
            mm3 = 0

            for i in range(len(infoList_current)):
                if infoList_current[i][3] == '0':
                    numPerfAlns = numPerfAlns + 1
                elif infoList_current[i][3] == '1':
                    mm1 = mm1 + 1
                elif infoList_current[i][3] == '2':
                    mm2 = mm2 + 1
                elif infoList_current[i][3] == '3':
                    mm3 = mm3 + 1
                else:
                    # should not see these, but you never know with BWA
                    mm3 = mm3 + 1

            #skip calculation if more tan one perfect alignment
            if numPerfAlns > 0:
                Database_file.write(header + ':' +  location[q] + '\t' + '\t' + guideRnaOutput[q] + '\t'  + direction[q] + '\t' + '\t' +  str(numPerfAlns) + '\t' + '\t' + '\t' + '1' + '\t' + '\t' + '\t' + str(GC[q]) + '\t' +  '\t' + str(mm1) + '\t' + str(mm2) + '\t' + str(mm3) + '\n')
                q = q +1
                #move onto next guide
                continue


            
            #if any of these conditions are true set specificity score to a low value approx 1.
            #Don't do any calculation and move onto next guide
            if mm1 > 5 or mm2 > 40 or mm3 > 500:
                Database_file.write(header + ':' +  location[q] + '\t' + '\t' + guideRnaOutput[q] + '\t'  + direction[q] + '\t' + '\t' +  str(numPerfAlns) + '\t' + '\t' + '\t' + '1' + '\t' + '\t' + '\t' + str(GC[q]) + '\t' +  '\t' + str(mm1) + '\t' + str(mm2) + '\t' + str(mm3) + '\n')
                q = q +1
                continue

            listIndex = 0
            #current guide from list
            guide_current = guideRna[q]
            #call function to complete calculation
            SGUIDE = specificity_calc(listIndex,infoList_current,guide_current,seq_h5)

            Database_file.write(header + ':' +  location[q] + '\t' + '\t' + guideRnaOutput[q] + '\t'  + direction[q] + '\t' + '\t' +  str(numPerfAlns) + '\t' + '\t' + '\t' + str(SGUIDE) + '\t' + '\t' + str(GC[q]) + '\t' +  '\t' + str(mm1) + '\t' + str(mm2) + '\t' + str(mm3) + '\n')
            q = q +1
    
        #this condition is met when BWA alinments excedes 1000 (we set this as the max to output)
        elif len(infoList_current) == 0:
            
            Database_file.write(header + ':' +  location[q] + '\t' + '\t' + guideRnaOutput[q] + '\t'  + direction[q] + '\t' + '\t' +  '-' + '\t' + '\t' + '\t' + '1' + '\t' + '\t'+ '\t' + str(GC[q]) + '\t' +  '\t' + 'too many' + '\t' + 'too many' + '\t' + 'too many' + '\n')
            q = q +1
            
           
      
    print "End of program"
    seq_h5.close()


main()        
