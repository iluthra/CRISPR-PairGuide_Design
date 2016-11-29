#Negative Control Generator
import random
import string
import os
import sys

#add same checks here for BSMBI fwd and rev sites
#rename depending on gene
outfile = open('Negative_control_oligos.txt', 'w')

def random_generator_20mer(num_random,length_seq):
    random_list_20mer = []

    z = 0
    while z < num_random:
        random_current = (''.join(random.choice('ATCG') for _ in range(length_seq)))
        
        random_list_20mer.append(random_current)
        z = z +1

    return random_list_20mer

def restriction_site_check(negative_guides, BsmBI, BsmBI_rev,spacer1,spacer2,u6):
    z =0
    random_20_noBsmbl = []
    while z < len(negative_guides):
        
        if negative_guides[z][0] == 'G':
            z = z +1
        elif negative_guides[z][0] != 'G':
            negative_guides[z] = 'G' +negative_guides[z]
            z = z +1
    h = 0
    while h < len(negative_guides):
        left_string = u6+negative_guides[h]+spacer1+BsmBI_rev[:-1]
        right_string = BsmBI[1:] + spacer2 + negative_guides[h] + scaffold
        if BsmBI in left_string:
            h = h + 1
        
        elif BsmBI_rev in left_string:
            h = h + 1
        elif BsmBI in right_string:
            h = h + 1
        elif BsmBI_rev in right_string:
            h = h + 1
        else:
            random_20_noBsmbl.append(negative_guides[h])
            h = h + 1
    print random_20_noBsmbl
    return random_20_noBsmbl

    

def random_generator_8mer(num_random,length_seq,BsmBI,BsmBI_rev):
    #need to check (BSMBI_rev-1bpfrom left + random8 + BSMBI-1bpfromright) for BSMBI fwd and rev
    random_list_8mer = []
   

    z = 0
    while z < num_random:
 

        random_current = (''.join(random.choice('ATCG') for _ in range(length_seq)))
        center_string = BsmBI_rev[1:] + random_current + BsmBI[:-1]
        totalgc = 0
        totalgc = random_current.count("G") + random_current.count("C")
        totalgc = int(totalgc)
        if totalgc < 4 or totalgc > 6:
            continue
        elif "AAAA" in random_current:
            continue
        elif "TTTT" in random_current:
            continue
        elif "GGGG" in random_current:
            continue
        elif "CCCC" in random_current:
            continue
        elif BsmBI in center_string:
            continue
        elif BsmBI_rev in center_string:
            continue
        else:
            random_list_8mer.append(random_current)
            z = z + 1

    print random_list_8mer
    return random_list_8mer





BsmBI_file = open("BsmBI.txt", "r")
input_info = []

for l in BsmBI_file:
    l = l.strip()
    input_info.append(l.split(':')[1])
BsmBI = input_info[0]

spacer1 = input_info[1]
spacer2 = input_info[2]
u6 = input_info[3]
scaffold = input_info[4]
#find reverse complement of BsmBI
CODE={'A':'T','T':'A','C':'G','G':'C','N':'N'} 
minus_seq=''
for c in BsmBI:
    minus_seq=minus_seq+CODE[c]
    BsmBI_rev=minus_seq[::-1]

Pair_guides = []
#add code here that creates negative control pairs and appends to Pair_guides list
#this number is arbitrary
num_random_oligos = 100000
random_20 = []
random_20_len = 20
z = 0
#create random guides for negative controls
random_20 = random_generator_20mer(num_random_oligos,random_20_len)

print random_20
j = 0
#create .fq file for alignment of negative controls
negative_control = open("Negative_control.fq",'w')
while j<len(random_20):

   
    negative_control.write('@' + str(j) + '\n')
    negative_control.write(random_20[j] + '\n')
    negative_control.write('+' + '\n')
    negative_control.write('IIIIIIIIIIIIIIIIIIII' + '\n')
    
    j = j +1


negative_control.close()

open("Negative_control.fq", "r")

cmd = 'bwa aln -o 0 -n 4 -N -t 5 hg37.fa Negative_control.fq > aln_GATA3_negative_controls.sai'
os.system(cmd)

print "Alignment Step has finished"

cmd2 = 'bwa samse -n 8000 hg37.fa aln_GATA3_negative_controls.sai Negative_control.fq > aln_GATA3_negative_controls.sam'
os.system(cmd2)

f2 = open("aln_GATA3_negative_controls.sam", "r")
print "Sam file has been opened"
import re
current_aln = 0
aln_negative = []
unmapped_num = 0
for line in f2:
        
        #print line
        if line.startswith('@'):
            continue


        alnList = re.findall('(chr[\dMXYIVLR]+),([+-])(\d+),\d+M,(\d)', line)
      
        i = 0
        for i in range(len(alnList)):
            if alnList[i][3] == '0' or alnList[i][3] == '1' or alnList[i][3] == '2' or  alnList[i][3] == '3':
                #count = count + 1
                break
            
            elif alnList[i][3] == '4':
                i = i + 1

        if i == len(alnList) & len(alnList) != 0:
            aln_negative.append(line)
           
            unmapped_num = unmapped_num + 1
       
    
            
print unmapped_num
negative_guides = []
f = 0
while f < len(aln_negative):
    negative_guides.append([x for x in aln_negative[f].split('\t')][9])
    f = f + 1
    
random_20_noBsmbl = restriction_site_check(negative_guides, BsmBI, BsmBI_rev,spacer1,spacer2,u6)
unmapped_num = len(random_20_noBsmbl)
print random_20_noBsmbl


random_8_len = 8

random_8 = random_generator_8mer(unmapped_num,random_8_len,BsmBI,BsmBI_rev)


print len(random_8)

spacer1_and_BsmBI_rev = spacer1 + BsmBI_rev
BsmBI_and_spacer2 = BsmBI + spacer2


b = 0

if (len(random_20_noBsmbl)%2) != 0:
    del random_20_noBsmbl[0]
    del random_8[0]


while b < len(random_20_noBsmbl):
    outfile.write(u6+ random_20_noBsmbl[b] +spacer1_and_BsmBI_rev + random_8[b/2] +BsmBI_and_spacer2 + random_20_noBsmbl[b+1] + scaffold+ '\n')
    b = b +2
print b/2
    

