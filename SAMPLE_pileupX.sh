#!/bin/sh

#PBS -l nodes=1

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

/home/mshahandeh/samtools-1.1/samtools mpileup -f /home/mshahandeh/D.sim_reference/dsimV2-Mar2012.fa -r X /home/mshahandeh/BCsynA/sim_bams/SAMPLE_sim-sorted.bam  > /home/mshahandeh/BCsynA/sim_bams/sim_pileups/SAMPLE.X_sim-sorted.pileup
