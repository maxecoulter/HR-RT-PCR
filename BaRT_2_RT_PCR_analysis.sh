
#!/bin/bash
#$ -cwd
#$ -j yes
#$ -pe smp 28
#

# Script to:
# First generate a transcriptome DB ;
# Then Blastn the Primers against the DB;
# Then analyse results


###Files names and paths:
DataBaseDB="Transcriptome_DB2"

#BlastN_Output="/mnt/shared/scratch/mc42302/201903_RTD2/Benchmarking/RT_PCR/Primers_vs_AllTranscripts_BaRT_1p_12"
#BlastN_Output="/mnt/shared/scratch/mc42302/201903_RTD2/Benchmarking/RT_PCR/Primers_vs_AllTranscripts_BaRT_Illumina3_12"
#BlastN_Output="/mnt/shared/scratch/mc42302/201903_RTD2/Benchmarking/RT_PCR/Primers_vs_AllTranscripts_BaRT_Iso4_12"
BlastN_Output="/mnt/shared/scratch/mc42302/201903_RTD2/Benchmarking/RT_PCR/Primers_vs_AllTranscripts_BaRT_2_10_12"


#TranscriptomeFileName="/mnt/shared/scratch/mc42302/201903_RTD2/Benchmarking/transcriptomes/BaRT.1.0_Unsplit_exons_padded.fasta"
TranscriptomeFileName="/mnt/shared/scratch/mc42302/201903_RTD2/Benchmarking/transcriptomes/BaRT_2.0.fasta"



PrimersFileName="/mnt/shared/scratch/mc42302/201903_RTD2/Benchmarking/RT_PCR/primers_RT_PCR1_updated.fasta" # the primer names must end with "_R" (reverse) and "_F" (forward) and must form pairs with the same name (example: >Primer101_R ; >Primer101_F
#PathToSamplesFolder="/mnt/shared/scratch/mc42302/201903_RTD2/Pacbio_20_samples/salmon_quants/quantsBaRT_1p_12/"
#PathToSamplesFolder="/mnt/shared/scratch/mc42302/201903_RTD2/Pacbio_20_samples/salmon_quants/quantsBaRT_Illumina3_12/"
#PathToSamplesFolder="/mnt/shared/scratch/mc42302/201903_RTD2/Pacbio_20_samples/salmon_quants/quantsBaRT_Iso4_12/"
PathToSamplesFolder="/mnt/shared/scratch/mc42302/201903_RTD2/Pacbio_20_samples/salmon_quants/quantsBaRT_2_10_12/"

#rtPCRInput="/mnt/shared/scratch/mc42302/201903_RTD2/Benchmarking/RT_PCR/rtPCR_productsAndProportions_updated.txt"  ## ## Header: "Primer 	Size	INF2	LEAF	EMBRYO	INF1	NODE"
#rtPCRInput="/mnt/shared/scratch/mc42302/201903_RTD2/Benchmarking/RT_PCR/rtPCR_productsAndProportions_updated.txt"
#rtPCRInput="/mnt/shared/scratch/mc42302/201903_RTD2/Benchmarking/RT_PCR/rtPCR_productsAndProportions_updated.txt"
rtPCRInput="/mnt/shared/scratch/mc42302/201903_RTD2/Benchmarking/RT_PCR/rtPCR_productsAndProportions_updated.txt"
#rtPCRInput="/mnt/shared/scratch/mc42302/201903_RTD2/Benchmarking/RT_PCR/rtPCR_productsAndProportions_subset2.txt"
scatterplot_output="/mnt/shared/scratch/mc42302/201903_RTD2/Benchmarking/RT_PCR/BaRT_2_10_12_w6"


conda activate blast
echo "Making blast db file from transcriptome..."
makeblastdb -in $TranscriptomeFileName -out $DataBaseDB -dbtype nucl

#### Run Blastn  - Primers against transcriptome database
echo "Blasting..."
date	
###			
###			
blastn -task blastn-short -num_threads 28 -query $PrimersFileName \
-db $DataBaseDB -out $BlastN_Output -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen salltitles qcovs qcovhsp"  

conda deactivate
### Remove database
rm Transcriptome_DB.n*

		

#Now carry out RT-PCR analysis based on BLAST results
conda activate ggplot_python
echo "Now analysing data..."

python BaRT_RT_PCR_analysis.py -b $BlastN_Output -p $rtPCRInput -s $PathToSamplesFolder -w 6 -cw 0 -o $scatterplot_output

conda deactivate

echo "script done"




#################################################################################################################################################






