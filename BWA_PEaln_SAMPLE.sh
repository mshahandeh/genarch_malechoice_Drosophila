# Make sure that we are in the same subdirectory as where the qsub command 
# is issued. 

cd $PBS_O_WORKDIR 

#  Determine the nodes, num process, etc.
cat $PBS_NODEFILE > nodes

# Get number of nodes allocated
NO_OF_NODES=`cat $PBS_NODBEFILE | egrep -v '^#'\|'^$' | wc -l | awk '{print $1}'`
NODE_LIST=`cat $PBS_NODEFILE `

# Just for kicks, see which nodes we got.
echo $NODE_LIST

/home/mshahandeh/bwa-0.6.2/bwa aln /home/mshahandeh/D.sim_reference/dsimV2-Mar2012.fa /home/mshahandeh/Data/SAMPLE_*_R1_001.fastq > /home/mshahandeh/Data/SAMPLE.R1.sai

/home/mshahandeh/bwa-0.6.2/bwa aln /home/mshahandeh/D.sim_reference/dsimV2-Mar2012.fa /home/mshahandeh/Data/SAMPLE_*_R2_001.fastq > /home/mshahandeh/Data/SAMPLE.R2.sai

/home/mshahandeh/bwa-0.6.2/bwa sampe /home/mshahandeh/D.sim_reference/dsimV2-Mar2012.fa /home/mshahandeh/Data/SAMPLE.R1.sai /home/mshahandeh/Data/SAMPLE.R2.sai /home/mshahandeh/Data/SAMPLE_*_R1_001.fastq /home/mshahandeh/Data/SAMPLE_*_R2_001.fastq > /home/mshahandeh/BCsynA/sim_sams/SAMPLE_sim.sam

/home/mshahandeh/samtools-1.1/samtools view -Sb /home/mshahandeh/BCsynA/sim_sams/SAMPLE_sim.sam > /home/mshahandeh/BCsynA/sim_bams/SAMPLE_sim.bam

/home/mshahandeh/samtools-1.1/samtools sort /home/mshahandeh/BCsynA/sim_bams/SAMPLE_sim.bam  /home/mshahandeh/BCsynA/sim_bams/SAMPLE_sim-sorted

/home/mshahandeh/samtools-1.1/samtools index /home/mshahandeh/BCsynA/sim_bams/SAMPLE_sim-sorted.bam

