#added ref as col 10 in python numbers

#base qual, 0-94, ! = 0:
#possible quals: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
possible_quals = '!\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'
cutoff = 14 #qualcutoff = 15

hist_of_ref_quals = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] #output what the quals are
hist_of_snp_quals = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] #output what the quals are
outf2 = open("/home/mshahandeh/BCsynA/sim_bams/sim_pileups/SAMPLE.3L_quals.txt","w")    

inf = open("/home/mshahandeh/BCsynA/sim_bams/sim_pileups/SAMPLE.3L_sim-sorted.pileup","r")


outf = open("/home/mshahandeh/BCsynA/sim_bams/sim_pileups/SAMPLE.3L_SNP.txt","w")    
   
 
old_chrom = ''

i = inf.readline()  	
if i != '': 				
    i = i.replace('\n','')
    i = i.split('\t')

done = 1 #this will kick in when there are no more line of data in the file

while done > 0:
    if i == '':
        done = 0
    else:
        options = [0,0,0,0,0,0] #a,g,c,t,d,i FOR PRINTING 
        chrom = i[0] #which reftig is this?
        if chrom != old_chrom:  #if we are on a new reftig, initialize as the first base of the reftig
            base = 1
            old_chrom = chrom
        i_base = int(i[1])
        ref = i[2]
        if i_base == base:   #here you are in the game with data. Otherwise, you will just print zeros for this base. 
            cov = i[3]
            ref = i[2]
            alleles = i[4]
            quals = i[5]
            which_qual_are_you_on = 0 #this will be used to pair up the qual score with the base call
            allele_list = []  #pre options list
            indel_state = 0
            possible_large_indel = 0
#            outf_ttt.write(chrom+'\t'+str(i[4])+'\n')
            for a in i[4]:
#                print indel_state
#                print a
                if possible_large_indel == 1:
                    possible_large_indel = 0
                    if a != 'A':   #this is in case the indel is >9 bases
                        if a != 'G':
                            if a != 'C':
                                if a != 'T':
                                    if a != 'c':
                                        if a != 'a':
                                            if a != 'g':
                                                if a != 't':
                                                    if a != 'n':
                                                        if a!= 'N':
                                                            sdf = indel_state * 10
                                                            indel_state = int(a) + sdf + 1
                if indel_state == -1:  #the last char was +-
                    indel_state = int(a) +1 #number of char to ignore
                    possible_large_indel = 1
                if a == '^':
                    indel_state = -20 #means ignore this next char, this char was carrot
                if indel_state == 0:
                    if a == '.':
                        this_qual = quals[which_qual_are_you_on]
                        for i in range(len(possible_quals)):
                            if this_qual == possible_quals[i]: 
                                hist_of_ref_quals[i] += 1
                                which_qual_are_you_on += 1
                                if i >= cutoff: #this means the qual is above the threshold
                                    allele_list.append(ref)
                    if a == ',':
                        this_qual = quals[which_qual_are_you_on]
                        for i in range(len(possible_quals)):
                            if this_qual == possible_quals[i]: 
                                hist_of_ref_quals[i] += 1
                                which_qual_are_you_on += 1
                                if i >= cutoff: #this means the qual is above the threshold
                                    allele_list.append(ref)                
                    if a == '+':
                        indel_state = -1
                        allele_list.append('I')
                    if a == '-':
                        indel_state = -1
                        allele_list.append('D') #I added this Jan 2015
                    if a == '*': #is splat even an option? I think this might be archaic
                        allele_list.append('D')
                    if a == 'a':
                        this_qual = quals[which_qual_are_you_on]
                        for i in range(len(possible_quals)):
                            if this_qual == possible_quals[i]: 
                                hist_of_snp_quals[i] += 1
                                which_qual_are_you_on += 1
                                if i >= cutoff: #this means the qual is above the threshold
                                    allele_list.append('A')
                    if a == 'A':
                        this_qual = quals[which_qual_are_you_on]
                        for i in range(len(possible_quals)):
                            if this_qual == possible_quals[i]: 
                                hist_of_snp_quals[i] += 1
                                which_qual_are_you_on += 1
                                if i >= cutoff: #this means the qual is above the threshold
                                    allele_list.append('A')
                    if a == 'g':
                        this_qual = quals[which_qual_are_you_on]
                        for i in range(len(possible_quals)):
                            if this_qual == possible_quals[i]: 
                                hist_of_snp_quals[i] += 1
                                which_qual_are_you_on += 1
                                if i >= cutoff: #this means the qual is above the threshold
                                    allele_list.append('G')
                    if a == 'G':
                        this_qual = quals[which_qual_are_you_on]
                        for i in range(len(possible_quals)):
                            if this_qual == possible_quals[i]: 
                                hist_of_snp_quals[i] += 1
                                which_qual_are_you_on += 1
                                if i >= cutoff: #this means the qual is above the threshold
                                    allele_list.append('G')
                    if a == 'c':
                        this_qual = quals[which_qual_are_you_on]
                        for i in range(len(possible_quals)):
                            if this_qual == possible_quals[i]: 
                                hist_of_snp_quals[i] += 1
                                which_qual_are_you_on += 1
                                if i >= cutoff: #this means the qual is above the threshold
                                    allele_list.append('C')
                    if a == 'C':
                        this_qual = quals[which_qual_are_you_on]
                        for i in range(len(possible_quals)):
                            if this_qual == possible_quals[i]: 
                                hist_of_snp_quals[i] += 1
                                which_qual_are_you_on += 1
                                if i >= cutoff: #this means the qual is above the threshold
                                    allele_list.append('C')
                    if a == 't':
                        this_qual = quals[which_qual_are_you_on]
                        for i in range(len(possible_quals)):
                            if this_qual == possible_quals[i]: 
                                hist_of_snp_quals[i] += 1
                                which_qual_are_you_on += 1
                                if i >= cutoff: #this means the qual is above the threshold
                                    allele_list.append('T')
                    if a == 'T':
                        this_qual = quals[which_qual_are_you_on]
                        for i in range(len(possible_quals)):
                            if this_qual == possible_quals[i]: 
                                hist_of_snp_quals[i] += 1
                                which_qual_are_you_on += 1
                                if i >= cutoff: #this means the qual is above the threshold
                                    allele_list.append('T')
                if indel_state > 0:
                    indel_state -= 1
                if indel_state == -10:  #means ignore this character, last char was carrot
                    indel_state = -0            
                if indel_state == -20:  #means ignore this next char, this char was carrot
                    indel_state = -10
#                print indel_state
            for j in allele_list:
                if j == 'A':
                    options[0] += 1
                if j == 'G':
                    options[1] += 1
                if j == 'C':
                    options[2] += 1
                if j == 'T':
                    options[3] += 1
                if j == 'D':
                    options[4] += 1
                if j == 'I':
                    options[5] += 1
                if j == 'a':
                    options[0] += 1
                if j == 'g':
                    options[1] += 1
                if j == 'c':
                    options[2] += 1
                if j == 't':
                    options[3] += 1
                if j == 'd':
                    options[4] += 1
                if j == 'i':
                    options[5] += 1
            
            #read next line:
            i = inf.readline()  	
            if i != '':
                i = i.replace('\n','')
                i = i.split('\t')
        else:
            cov = 0
            ref = 'O'
        outf.write(chrom+'\t'+str(base)+'\t'+str(cov)+'\t'+str(options[0])+'\t'+str(options[1])+'\t'+str(options[2])+'\t'+str(options[3])+'\t'+str(options[4])+'\t'+str(options[5])+'\t'+str(ref)+'\n')
#chrom, base, cov, a,g,c,t,d,i, ref(but only the true ref for bases with non-zero coverage!!)
        base +=1


for i in range(len(hist_of_ref_quals)):
    ref = hist_of_ref_quals[i]
    snp = hist_of_snp_quals[i]
    outf2.write(str(i)+'\t'+str(ref)+'\t'+str(snp)+'\n')


outf.close()
outf2.close()
