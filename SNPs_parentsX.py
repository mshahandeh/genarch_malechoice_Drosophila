in_file_simC4_synA = open("/home/mshahandeh/BCsynA/sim_bcfs/simC4_synA.varX_SNPs-final.vcf", "r")
out_file_simC4_synA = open("/home/mshahandeh/BCsynA/sim_bcfs/simC4_synA.X_SNPs.txt" , "w")

done = 1
while done > 0: 
    i = in_file_simC4_synA.readline()
    if i == '':
        done = 0
    elif i[0][0] == '#':
        continue
    else:
        i = i.split('\t')
        chrom = i[0]
        pos = i[1]
        ref = i[3]
        alt = i[4]
        sim_SNP = i[9]
        sech_SNP = i[10]
        sim_SNP = sim_SNP.replace(':' , '\t')
        sim_SNP = sim_SNP.replace(',' , '\t')
        sim_SNP = sim_SNP.replace('\n' , '')
        sim_SNP = sim_SNP.split('\t')
        sech_SNP = sech_SNP.replace(':' , '\t')
        sech_SNP = sech_SNP.replace(',' , '\t')
        sech_SNP = sech_SNP.replace('\n' , '')
        sech_SNP = sech_SNP.split('\t')
        if sech_SNP[0] == '0'  and sech_SNP[1] == '0' and int(sech_SNP[2]) >= 30 and sim_SNP[0] == '1' and int(sim_SNP[1]) >= 30 and sim_SNP[2] == '0':
            SNP = chrom + '\t' + pos + '\t' + ref + '\t' + alt + '\t' + '1' + '\t' + '0' + '\n'
        elif sim_SNP[0] == '0'  and sim_SNP[1] == '0' and int(sim_SNP[2]) >= 30 and sech_SNP[0] == '1' and int(sech_SNP[1]) >= 30 and sech_SNP[2] == '0':
            SNP = chrom + '\t' + pos + '\t' + ref + '\t' + alt + '\t' + '0' + '\t' + '1' + '\n'
        elif sech_SNP[0] == '0'  and sech_SNP[1] == '0' and sech_SNP[2] == '0' or sim_SNP[0] == '0'  and sim_SNP[1] == '0' and sim_SNP[2] == '0':
            continue
        else:
            continue
        out_file_simC4_synA.write(SNP)

in_file_simC4_synA.close()
out_file_simC4_synA.close()