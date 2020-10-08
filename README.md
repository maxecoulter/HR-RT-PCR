# HR-RT-PCR

Author: Max Coulter 2020

Contact: mecoulter@dundee.ac.uk

An updated and expanded version of the script described here (https://github.com/PauloFlores/RNA-Seq-validation) with the main matching and analysis in a python script.
This script uses a similar method to that described in the link, but with some changes:

* The python script carries out correlation analysis and produces figures as well as providing matches
* Matching RT-PCR products to transriptome products is a non-trivial problem in some cases. Therefore the algorithm for matching products together has been improved, by considering all possible combinations of matches with a specified window size and choosing the best match. Even then, other complex considerations from the raw RT-PCR data may need to be taken into account with matches added/removed. There is therefore an option to edit the product matches manually

* The script produces a correlation plot and spearman and pearson correlations as output, as well as RT-PCR:transcriptome product matches (if not using manual mode), transcriptome products identified for each primer set, and correlations for each individual primer set for troubleshooting

### Before you run the script

The script is for comparing RT-PCR quantifications with salmon quantifications to determine accuracy of transcritome based quantifications. Therefore an input file with RT-PCR base product proportions is needed (details of format described below).

You need to run salmon for transcript quantification on your samples with your chosen reference transcriptome. You also need to have BLAST version 2.5.0+ installed. 

The bash script provided (**BaRT_2_RT_PCR_analysis.sh**) runs BLASTn to identify possible matches for the list of primers provided. Once BLAST is complete the python script (**BaRT_RT_PCR_analysis.py**) will run and analyse the results from BLAST, Salmon and RT-PCR results. All input files with paths should be adjusted as required within the bash script.

### Input variables to be set:

1. **BlastN_Output** this is the location and name of the BLAST output you need to specify in line 17 of the bash script
1. **TranscriptomeFileName** this is the transcriptome you used for running Salmon. This is specified in line 19 of the bash script
1. **PrimersFileName** the name and path to the primers used in RT-PCR analysis, specified in line 21 of the bash script. Primer names must end with F or R (for forward and reverse) e.g:

        >53F 
        GCAGTTAATCTCCTCTGG 
        >53R 
        GGCGTAGGTGGACGCGGTGAGG
    

1.**PathToSamplesFolder** this is the path to folder with salmon quantifications in it, specified in line 23 of bash script. Quantifications for each sample should be in folders in this directory, make sure the sample names for each folder match with the sample names given in RT-PCR input file.
1.**rtPCRInput** this is the path and input file with RT PCR product proportions, specified in line 25 of the bash script. The tsv file should be in the following format:

                Primer	Size	sample1	sample2	sample3	sample4	sample5	sample6	sample7	sample8	sample9	sample10	sample11	sample12
                53	130		0.245169601	0.149009806	0.162643258	0.224078159		0.172184007	0.179695698	0.188429752	0.184210526	0.142186166	0.12686687
                53	223		0.754830399	0.850990194	0.837356742	0.775921841		0.827815993	0.820304302	0.811570248	0.815789474	0.857813834	0.87313313
                54	110	0.241532767	0.41989802	0.452108593	0.440929428	0.501279834	0.23851939	0.464493683	0.493283391	0.301314199	0.503980049	0.390434631	0.553052
                54	332	0.758467233	0.58010198	0.547891407	0.559070572	0.498720166	0.76148061	0.535506317	0.506716609	0.698685801	0.496019951	0.609565369	0.446948

