# Head transcriptome differencies between males swarming and non-swarming mosquitoes (<i>Aedes aegypti</i>)
Whole head RNA-sequencing of swarming and non-swarming yellow fever mosquito (<i>Aedes aegypti</i>) males analysis.  

## Sequencing QC and read alignment to the reference transcriptome

The '01_snakemake_alignemnt_and_QC' folder contain the Snakemake code (Mölder et al., 2021) to run a FastQC (Andrews, 2015), align reads to the reference transcriptome using Kallisto () and sum up the results with MultiQC (). 
To run it, add the raw reads with correct names in the data/samples folder, and the reference transcriptome as data/ref_genome.fasta. In this analysis we used the following one: 

## Differential transcript expression analysis

The code '02_transcript_DEG_analysis.R' contain the analysis of RNA-seq at the trancript level. 

## Citations 

Andrews,S. (2010) FastQC: A Quality Control Tool for High Throughput Sequence Data.

Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., 2021. Sustainable data analysis with Snakemake. F1000Res 10, 33.






