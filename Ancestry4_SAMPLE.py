
windsize = 1000000
by = 100000

#these are the ends of the simC4_synA SNP files

len4 = 1100000

#inf1:
#chrom	 pos	ref	alt	simC4	synA
#4	1073	C	G	0	1
#4	1615	T	C	0	1
#4	1665	A	C	0	1
#4	1738	T	A	0	1
#4	1763	C	T	0	1
#4	2155	C	T	0	1
#4	2246	C	T	0	1

#inf2:
#chrom	pos		cov		a		g		c		t		i		d		ref
#4      28      6       0       0       5       0       0       0       C
#4      29      6       6       0       0       0       0       0       A
#4      30      6       6       0       0       0       0       0       A

#inf = open("/Users/michaelshahandeh/Documents/Sequencing/BCsynA_scripts/4_Analysis/simC4_synA.4_SNPs.txt", "r")
#inf2 = open("/Users/michaelshahandeh/Documents/Sequencing/BCsynA_scripts/4_Analysis/1.4_SNP.txt", "r")

inf = open("/home/mshahandeh/BCsynA/sim_bcfs/simC4_synA.4_SNPs.txt","r") 
inf2 = open("/home/mshahandeh/BCsynA/sim_bams/sim_pileups/SAMPLE.4_SNP.txt","r")    

#output will be the % synA ancestry, so 0 = simC4, 1 = synA, options in between

outf = open("/home/mshahandeh/BCsynA/Ancestry/4/ancestry_SAMPLE_4.txt","w")    
#outf = open("/Users/michaelshahandeh/Desktop/Ancestry_14.txt","w") 
#outf2 = open("/Users/michaelshahandeh/Desktop/rawdata_14.txt" , "w")
outf2 = open("/home/mshahandeh/BCsynA/Ancestry/4/raw/rawdata_SAMPLE_4.txt" , "w")

outf.write('SAMPLE' + '\n')
outf2.write('SAMPLE' + '\n')
 
# ---------------------------------------------------------------------------------------------------
 #first: read in a list of marker states and positions (snps different between simC4_synA in my case)
positions = []
simC4 = []
synA = []


done = 1 #this will kick in when there are no more line of data in the file
while done > 0:
    i = inf.readline()  	   
    if i == '':
        done = 0
    else:
        i = i.replace('\n','')
        i = i.split('\t')
        positions.append(i[1])
        if i[4] == '1':
        	simC4.append(i[3])
        if i[4] == '0':
        	simC4.append(i[2])
        if i[5] == '0':
        	synA.append(i[2])
        if i[5] == '1':
        	synA.append(i[3])
# ---------------------------------------------------------------------------------------------------

this_line_positions = []
this_line_ancestry = []

count = 0 #which marker are you on?
done = 1 
while done > 0:
    i = inf2.readline() #this file has one line for each base in the genom, see above
    if i == '': #if you have run out of chromosome, stop (this should be redundant)
        done = 0
    if count == len(positions): #if you have run out of markers, stop
        done = 0
    if done == 1:
        i = i.replace('\n','')
        i = i.split('\t')
        pos = int(i[1])      #what locus is this?
        if int(positions[count]) == pos:  #if this locus == the marker you are currently looking for, do the following
            #outf2.write(str(positions[count])+'\t'+str(simC4[count])+'\t'+str(synA[count])+'\t'+str(i)+'\n')  #this was in here as a debugging aid. Print and write statements are key.
            ancestry = 'NA'  #default ancestry = NA
            a =  float(i[3])  #the read counts for each base
            g =  float(i[4])
            c =  float(i[5])           
            t =  float(i[6])      
            if synA[count] == 'A':  #if the synA parent has an A at this position,
                if simC4[count] == 'G':  #and the simC4 parent has a G,
                    if (a+g) > 0:  #and you aren't going to divide by zero,
                        ancestry = a / (a+g)  #then this is the frequency of the synA snp. Usually 1 or 0 at low coverage, but sometimes .5 or .33 and that is useful info. Note that I use a+g as coverage, not a+g+c+t
                if simC4[count] == 'C':
                    if (a+c) > 0:
                        ancestry = a / (a+c)
                if simC4[count] == 'T':
                    if (a+t) > 0:
                        ancestry = a / (a+t)


            if synA[count] == 'G':
                if simC4[count] == 'A':
                    if (a+g) > 0:
                        ancestry = g / (a+g)
                if simC4[count] == 'C':
                    if (c+g) > 0:
                        ancestry = g / (g+c)
                if simC4[count] == 'T':
                    if (g+t) > 0:
                        ancestry = g / (g+t)

            if synA[count] == 'C':
                if simC4[count] == 'A':
                    if (a+c) > 0:
                        ancestry = c / (a+c)
                if simC4[count] == 'G':
                    if (g+c) > 0:
                        ancestry = c / (g+c)
                if simC4[count] == 'T':
                    if (t+c) > 0:
                        ancestry = c / (t+c)

            if synA[count] == 'T':
                if simC4[count] == 'A':
                    if (a+t) > 0:
                        ancestry = t / (a+t)
                if simC4[count] == 'C':
                    if (c+t) > 0:
                        ancestry = t / (c+t)
                if simC4[count] == 'G':
                    if (g+t) > 0:
                        ancestry = t / (g+t)

            if ancestry != 'NA':  #if you had some coverage at this location
                this_line_positions.append(pos)  #make a list of loci + freq of synA allele
                this_line_ancestry.append(ancestry)
            count+=1  #move on to the next marker

#debugging aids
#for i in range(len(this_line_positions)):
#    outf2.write(str(this_line_positions[i])+'\t'+str(this_line_ancestry[i])+'\n')
#outf2.close()



#------------------------------------------------------------------------
#windows:

wind_pos = [] #these will be output
wind_val = []

start = 0  #initializing variables
end = start + windsize
done = 0
temp2 = []
i = 0

while done < 1: # this kicks in when you are out of markers
    first = 'init' #you are starting a new window
    done2 = 0
    while done2 < 1: #for this window
        #print(i) #debugging
        if i >= (len(this_line_positions))-1: #if you are out of markers, stop
            done2 = 1
        else:
            if this_line_positions[i] > end: #if this locus is past the end of the window you are working on, stop
                done2 = 1
            else:
                if this_line_positions[i] > start: #now we are in the window
                    if first == 'init': #you just started a new window, so determine where it started
                        first = i  #to speed things up, I will start searching for the next window at the beginning of the previous window. This will work even if overlap is 99%; I am happy with this sol'n
                    if this_line_ancestry[i] != 'NA':  #you are in the window, so if != NA, append to temp2
                        temp2.append(this_line_ancestry[i])   
        i += 1 #on to the next locus
        
    #now you just finished a window -- you were kicked out of loop above b/c you are past the end value
    if len(temp2) > 9: #require data at 10 or more markers in a window
        mid = (start+end)/float(2)
        wind_pos.append(mid)
        ave = sum(temp2) / len(temp2)
        wind_val.append(ave)
        #print("A",start, end,mid,ave) #debugging
    else:
        mid = (start+end)/float(2)
        wind_pos.append(mid)
        wind_val.append('NA')
        #print("B",start, end)
    start += by #slide the window
    end += by
    temp2 = []
    if first != 'init': #this is incase we are at the end. There was a bug I had to deal with.
        i = first #reset the start position
    if end > len4:  ########################################################  #THIS IS JENKY: I have to manually input these values at the top
        done = 1

#print
for i in range(len(wind_pos)):
    if wind_val[i] == 'NA':
        outf.write('-'+'\n')
    else:
        if wind_val[i] > 0.75:
            a = 'A'
        if wind_val[i] <= 0.75:
            a = 'H'

        outf.write(a+'\n')
        outf2.write(str(wind_pos[i])+ '\t' + str(wind_val[i]) + '\n')
outf.close()
outf2.close()

