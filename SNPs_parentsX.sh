# Make sure that we are in the same subdirectory as where the qsub command 
# is issued. 

cd $PBS_O_WORKDIR 	

#  Determine the nodes, num process, etc.

cat $PBS_NODEFILE > nodes

# Get number of nodes allocated
NO_OF_NODES=`cat $PBS_NODBEFILE | egrep -v '^#'\|'^$' | wc -l | awk '{print $1}'`
NODE_LIST=`cat $PBS_NODEFILE `#!/bin/sh

#PBS -l nodes=1

# Just for kicks, see which nodes we got.
echo $NODE_LIST

#samtools version of indexing the reference
/home/mshahandeh/samtools-1.1/samtools faidx /home/mshahandeh/D.sim_reference/dsimV2-Mar2012.fa

#generate a single pileup from multiple files for X:8,950,000-13,660,000 region; -Q = min base quality (30 = 1 in 1000 errors); -q = min read quality (15 = .05); -I = skip indel calling
/home/mshahandeh/samtools-1.1/samtools mpileup -C 50 -B -Q 30 -q 15 -I -vf /home/mshahandeh/D.sim_reference/dsimV2-Mar2012.fa -r X /home/mshahandeh/BCsynA/sim_bams/simC4_sim-sorted.bam /home/mshahandeh/synA.Full_sim-sorted.bam.bam > /home/mshahandeh/BCsynA/sim_bcfs/simC4_synA_SNPs-rawX.vcf

#Call SNPs -S followed by two column text file with sample names (FULL bam path) and ploidy (0-2)
/home/mshahandeh/bcftools-1.1/bcftools call -X -S /home/mshahandeh/BCsynA/BCsynAscripts/SNPs_simC4_synA_X.txt -vc /home/mshahandeh/BCsynA/sim_bcfs/simC4_synA_SNPs-rawX.vcf > /home/mshahandeh/BCsynA/sim_bcfs/simC4_synA_SNPs-varX.vcf

#Filter SNPs
/home/mshahandeh/bcftools-1.1/bcftools view /home/mshahandeh/BCsynA/sim_bcfs/simC4_synA_SNPs-varX.vcf | /home/mshahandeh/bcftools-1.1/vcfutils.pl varFilter > /home/mshahandeh/BCsynA/sim_bcfs/simC4_synA.varX_SNPs-final.vcf
