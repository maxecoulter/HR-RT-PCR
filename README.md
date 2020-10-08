# HR-RT-PCR

Author: Max Coulter 2020
Contact: mecoulter@dundee.ac.uk

An updated and expanded version of the script described here (https://github.com/PauloFlores/RNA-Seq-validation) with the main matching and analysis in a python script.
This script uses a similar method to that described in the link, but with some changes:

* The python script carries out correlation analysis and produces figures as well as providing matches
* Matching RT-PCR products to transriptome products is a non-trivial problem in some cases. Therefore the algorithm for matching products together has been improved, by considering all possible combinations of matches with a specified window size and choosing the best match. Even then, other complex considerations from the raw RT-PCR data may need to be taken into account with matches added/removed. There is therefore an option to edit the product matches manually

* The script produces a correlation plot and spearman and pearson correlations as output, as well as RT-PCR:transcriptome product matches (if not using manual mode), transcriptome products identified for each primer set, and correlations for each individual primer set for troubleshooting

### Before you run the script

You need to run salmon for transcript quantification on your samples with your chosen reference transcriptome. You also need to have BLAST version 2.5.0+ installed. 

The bash script provided runs BLASTn to identify possible matches for the list of primers provided. Once BLAST is complete the python script will run and analyse the results from BLAST, Salmon and RT-PCR results. All input files with paths should be adjusted as required within the bash script.

### Input files:

1. **BlastN_Output** this is the location and name of the BLAST output you need to specify in line 17 of the bash script
1. **TranscriptomeFileName** this is the transcriptome you used for running Salmon. This is specified in line 19 of the bash script
1. 

