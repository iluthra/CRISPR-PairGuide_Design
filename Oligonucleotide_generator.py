#final oligo generator
import random
import string
import os
import sys

#open file created in Guide_pairing.py
pairs = open("Guide_Pairs_S=65_D=1000s_GLs_varying_GATA3.txt","r")
#create file that has all info per pair on each line
outfile_formated = open("Formated_GuidePairs_GLs_varying_GATA3", 'w')
#rename depending on gene
f = open('Final_oligos_GATA3_GLS_varying.txt', 'w')

def random_generator(num_random,length_seq,BsmBI,BsmBI_rev):
    #need to check (BSMBI_rev-1bpfrom left + random8 + BSMBI-1bpfromright) for BSMBI fwd and rev
    random_list = []


    z = 0
    while z < num_random:
        #add check GC content and 4bp repeats
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
            random_list.append(random_current)
            z = z + 1

 
    return random_list


def main():
    Pair_guides = []
    #open file with BsmBI, spacers, u6 and scaffold
    BsmBI = open("BsmBI.txt", "r")
    input_info = []
    
    for l in BsmBI:
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

    for l in pairs:

        Pair_guides.append(l.split('\t'))




    k = 1
    pair_guides_clean = []
    #filter out empty lists or list with just '\n'
    while k < len(Pair_guides):
        if len(Pair_guides[k]) > 16:
            pair_guides_clean.append(Pair_guides[k])

        k = k +1

    #delete header line
    if pair_guides_clean[0][2] == 'Location':

        del pair_guides_clean[0]
    
    print len(pair_guides_clean)/2
    random_8_list = []
    num_random_oligos = len(pair_guides_clean)/2
    count = 0
    random_8_len = 8
    #generative 8bp randomers
    random_8 = random_generator(num_random_oligos,random_8_len,BsmBI,BsmBI_rev)

    spacer1_and_BsmBI_rev = spacer1 + BsmBI_rev
    BsmBI_and_spacer2 = BsmBI + spacer2



    b = 0
    oligos = []
    #write oligos to file
    while b < len(pair_guides_clean):
        oligos.append(u6+ pair_guides_clean[b][3]+spacer1_and_BsmBI_rev + random_8[b/2] +BsmBI_and_spacer2 + pair_guides_clean[b+1][3] + scaffold)
        f.write(u6+ pair_guides_clean[b][3]+spacer1_and_BsmBI_rev + random_8[b/2] +BsmBI_and_spacer2 + pair_guides_clean[b+1][3] + scaffold+ '\n')
        b = b +2
    print b/2

    print len(oligos)
    
    outfile_formated.write('Location1' + '\t'+'Location2' + '\t' + 'gRNA1' + '\t'+ 'gRNA2' + '\t'  + 'Direction1' + '\t'+ 'Direction2' + '\t' + 'other.perfect.aligns1' + '\t'+ 'other.perfect.aligns2' + '\t' +  'Specificity.Score1' + '\t' +  'Specificity.Score2' + '\t' + '%GC.Content1' + '\t' + '%GC.Content2' + '\t' + 'g1#MM1'+ '\t' + 'g2#MM1'+ '\t' + 'g1#MM2'+  '\t' + 'g2#MM2' +'\t' + 'g1#MM3' +'\t' + 'g2#MM3' + '\t' + 'C.Start1' + '\t' + 'C.Start2' + '\t' + 'C.End1'+ '\t' + 'C.End2'+ '\t' + 'ChromFree1' + '\t' + 'ChromFree2'+  '\n')

    #print pairs[0]
    d=0
    print len(oligos)
    
    #write to formated file, easy to paste into excel and check pairs created
    while d < len(pair_guides_clean):
        p = 0
        pair_guides_clean[d] = filter(None, pair_guides_clean[d])
        pair_guides_clean[d+1] = filter(None, pair_guides_clean[d+1])
        while p < len(pair_guides_clean[d]):
            pair_guides_clean[d][p] = pair_guides_clean[d][p].strip()
            pair_guides_clean[d+1][p] = pair_guides_clean[d+1][p].strip()
            if pair_guides_clean[d][p] != '\n':
                outfile_formated.write(str(pair_guides_clean[d][p]) + '\t '+ str(pair_guides_clean[d+1][p]) + '\t')
                #p = p +1
            p = p +1
        outfile_formated.write(oligos[d/2])
        outfile_formated.write('\n')
        d = d +2
    

main()
