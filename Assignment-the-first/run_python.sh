#!/bin/bash
#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=compute               #REQUIRED: which partition to use
#SBATCH --mail-user=sgolubev@uoregon.edu     #optional: if you'd like email
#SBATCH --mail-type=ALL                   #optional: must set email first, what type of email you want
#SBATCH --cpus-per-task=4                 #optional: number of cpus, default is 1
#SBATCH --mem=32GB                        #optional: amount of memory, default is 4GB
conda activate base
/usr/bin/time -v python initial_data_exploration.py -f '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz' -l 101 -his 'histR1'