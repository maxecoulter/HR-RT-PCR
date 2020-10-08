
#!/bin/bash
#$ -cwd
#$ -j yes
#$ -pe smp 28


# Script to:
# First generate a transcriptome DB ;
# Then Blastn the Primers against the DB;
# Then analyse results


###Files names and paths:
DataBaseDB="Transcriptome_DB"

BlastN_Output="/path/to/Primers_vs_AllTranscripts"

TranscriptomeFileName="/path/to/Transcriptome.fasta"

PrimersFileName="/path/to/primers.fasta" # the primer names must end with "R" (reverse) and "F" (forward) and must form pairs with the same name (example: >101R ; >101F)

PathToSamplesFolder="/path/to/salmon_quants/"

rtPCRInput="/path/to/rtPCR_productsAndProportions.txt"  ## ## Header: "Primer 	Size	INF2	LEAF	EMBRYO	INF1	NODE"


manual="/path/to/RT_PCRcluster_results.txt" #Manually edit the product match file, optional

output="/path/to/output_prefix" #Output prefix



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

python BaRT_RT_PCR_analysis.py -b $BlastN_Output -p $rtPCRInput -s $PathToSamplesFolder -w 6 -o $output

conda deactivate

echo "script done"




#################################################################################################################################################






