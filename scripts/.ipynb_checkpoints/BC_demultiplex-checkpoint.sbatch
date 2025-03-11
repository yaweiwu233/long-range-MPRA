#!/usr/bin/bash

#SBATCH --job-name=BC_demultiplex
#SBATCH -N 1
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=4
#SBATCH -o slurm-%x_%A_%a.out
#SBATCH --mail-type=END
#SBATCH --mail-user=yawei.wu@wustl.edu

sample_ID=sample_ID.txt
eBC=total_library_final_with_BC.txt
dBC=dBCs.txt
input_folder=../fastq/GATA1_lib__HBG_2LP
output=../data/GATA1_lib_HBG_2LP

for sample in $(cat $sample_ID)
do
	python3 BC_demultiplex_supervised.py -f $input_folder/$sample -e $eBC -d $dBC -o $output --eBC-pattern "ACTAGT([ATCG]{10})CTCGAG[ATCGN]*GGATCC([ATCG]{16})CCTGCAGG" --dBC-pattern "CTTAAG([ATCG]{12})"
done
