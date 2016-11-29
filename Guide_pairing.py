#this code will pair the single guides
import bisect
import numpy as np
#open either filtered or chromatin database


def create_graphing_files(pairs,region_start,region_end):
    cut_column = 15
    #name files depending on S and D used
    outfile2 = open("Cut_length_S=65_D=1000s_GLs_varying_GATA3.txt" , 'w')
    outfile = open("Step_size_S=65_D=1000s_GLs_varying_GATA3.txt" , 'w')
    outfile_coverage = open("Coverage_of_guides_S=60_D=1000s_GLs_varying_GATA3.txt", 'w')
    stepsize = []
    m = 0
    #write file with all step sizes
    while m < (len(pairs)-2):
        temp = 0
        final = int(pairs[m+2][cut_column])
        initial = int(pairs[m][cut_column])
        temp = final - initial

        stepsize.append(temp)
        outfile.write(str(temp) + '\n')
        m = m + 2

    k = 0

    cut_length = []

    #create file for cut lenth to graph in R
    #cut length = PairGuide2cutStart - PairGuide1CutStart
    while k < len(pairs)-1:
        temp2 = int(pairs[k+1][cut_column]) - int(pairs[k][cut_column])
        cut_length.append(temp2)
        outfile2.write(str(temp2) + '\n')
        k = k + 2


    #create index list for input region

    region_len = region_end - region_start + 1 
    sum_overlap = np.zeros(region_len)
    index = range(region_start,region_end)


    #parse through single guides list and find how many guides overlap each base in input region
    f = 1
    overlap = []
    long_list = []

    f = 0

    z = 0
    # loop over guide-pairs to determine deletion coverage
    while z < len(pairs):

        start_index = int(pairs[z][cut_column])- region_start
        end_index = int(pairs[z+1][cut_column]) - region_start + 1
        sum_overlap[start_index:end_index] += 1
        z = z +2


    b = 0

    while b < len(index):

        outfile_coverage.write(str(index[b]) + '\t' + str(sum_overlap[b]) + '\n')
        b = b +1

    return None







def main():


    Database = open("GATA3_Database_ChromatinRegions_GLs_varying.txt","r")
    #declare empty lists
    SingleGuides_total = []
    SingleGuides = []

    f = open("Guide_Pairs_S=65_D=1000s_GLs_varying_GATA3.txt" , 'w')

    #change these depeding on input region
    region_start   = 7096667
    region_end = 9117164

    #create a list of lists that have all info. about single guides
    for l in Database:

        SingleGuides_total.append(l.split('\t'))

    pairs = []
    #these are sequnces that cannot be present in any of our guides
    i = 1
    #read this in from a text file
    BsmBI_file = open("BsmBI.txt", "r")
    input_info = []


    for l in BsmBI_file:
        l = l.strip()
        input_info.append(l.split(':')[1])

    #store strings from input file
    BsmBI = input_info[0]
    spacer1 = input_info[1]
    spacer2 = input_info[2]
    u6 = input_info[3]
    scaffold = input_info[4]

    #calculate reverse complement for BsmBI
    CODE={'A':'T','T':'A','C':'G','G':'C','N':'N'} 
    minus_seq=''
    for c in BsmBI:
        minus_seq=minus_seq+CODE[c]
        BsmBI_rev=minus_seq[::-1]

    print "#Guides before removing restrictions"
    print len(SingleGuides_total)


    count = 0
    #create new list that doesn't have guides with "bad sequences"
    while count < len(SingleGuides_total):
            left_string = u6+SingleGuides_total[count][2]+spacer1+BsmBI_rev[:-1]
            right_string = BsmBI[1:] + spacer2 + SingleGuides_total[count][2] + scaffold
            if BsmBI in left_string:
                    count = count +1
                    continue
            elif BsmBI_rev in left_string:
                    count= count+1
                    continue
            elif BsmBI in right_string:
                    count= count+1
                    continue
            elif BsmBI_rev in right_string:
                    count= count+1
                    continue
            else:
                    SingleGuides.append(SingleGuides_total[count])
                    count = count +1


    print "#Guides after removing restrictions"
    print len(SingleGuides)



    #Set your step and deletion sizes here
    S = 65
    D = 1000
    cut_column = 15
    while i < len(SingleGuides)-1:
            #we always start with the first guide, so we append that first always

            if i == 1:
                    #print "He"
                    last_cut = 0
                    pairs.append(SingleGuides[i])

                    cut_start = int(SingleGuides[i][cut_column])
                    last_cut = cut_start
                    #minimum cut distance
                    next_min = cut_start + D
                    j = i+1
                    while j < len(SingleGuides):

                            next_cut = int(SingleGuides[j][cut_column])
                            if  next_cut >= next_min:
                                    pairs.append(SingleGuides[j])
                                    break
                            j = j +1
                    if len(pairs) != (2*i):
                            del pairs[i-1]
                    i = i +1

            if i > 1:

                    #only add the next pair after the first if it the cut start site is at least step size away from the previous one
                    if (int(SingleGuides[i][15]) - last_cut) >= S:

                            cut_start = int(SingleGuides[i][cut_column])

                            next_min = cut_start + D
                            j = i+1
                            while j < len(SingleGuides):

                                    next_cut = int(SingleGuides[j][cut_column])
                                    if  next_cut >= next_min:
                                            pairs.append(SingleGuides[i])
                                            pairs.append(SingleGuides[j])
                                            last_cut = int(SingleGuides[i][cut_column])
                                            break
                                    j = j +1


                    i = i +1


    print "The number of pairs are"
    print len(pairs)/2

    #write file that has all paired guides


    k = 0
    f.write('Guide Pairs' + '\n')
    f.write('\t')
    f.write('\t'+'Location' + '\t' + '\t' + '\t' + '\t' + 'gRNA' + '\t'  + 'Direction' + '\t' + 'other.perfect.aligns' + '\t' +  'Specificity.Score' + '\t' + '%GC.Content' + '\t' + '#MM1'+ '\t' + '#MM2'+ '\t' + '#MM3' + '\t' + 'C.Start' + '\t' + 'C.End'+ '\t' + 'ChromFree'  '\n')
    f.write('\t')
    line_counter = 0
    for element in pairs:

        line_counter = line_counter + 1
        for item in element:
            item.strip('\n')
            f.write(str(item)+'\t')
        #add spacer line between each pair of guides
        if line_counter % 2 == 0:

            f.write('\n')
            f.write('\t')

    #this function will create the files needed for graphing data
    #set file names in the function
    files = create_graphing_files(pairs,region_start,region_end)
    

main()
