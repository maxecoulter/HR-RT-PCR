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
2. **TranscriptomeFileName** this is the transcriptome you used for running Salmon. This is specified in line 19 of the bash script
3. **PrimersFileName** the name and path to the primers used in RT-PCR analysis, specified in line 21 of the bash script. Primer names must end with F or R (for forward and reverse) e.g:

        >53F 
        GCAGTTAATCTCCTCTGG 
        >53R 
        GGCGTAGGTGGACGCGGTGAGG
    

4. **PathToSamplesFolder** this is the path to folder with salmon quantifications in it, specified in line 23 of bash script. Quantifications for each sample should be in folders in this directory, make sure the sample names for each folder match with the sample names given in RT-PCR input file.
5. **rtPCRInput** this is the path and input file with RT PCR product proportions, specified in line 25 of the bash script. The tsv file should be in the following format:

                Primer	Size	sample1	sample2	sample3	sample4	sample5	sample6
                53	130	0.240	0.149	0.163	0.224	0.172	0.179	
                53	223	0.754	0.850	0.837	0.775	0.827	0.820
                54	110	0.241	0.419	0.452	0.440	0.501	0.238
                54	332	0.758	0.580	0.547	0.559	0.498	0.761
             
  a) Column 1, primer name without F or R termination;
  b) Column 2, product size;
  c) Columns 3,4... Individual per sample proportions for each product size

6. **manual** this is an optional file for specifying RT-PCR-transcriptome product matches, specified in line 28 of bash script. The file with the correct format is produced automatically when the script is first run output file name has the suffix (**RT_PCRcluster_results**) This should be checked after the first run to see if matches are correct. If adjustments are required, this file should be adjusted. The file is in this format:

                Primer	RT PCR product	Matched transcriptome products
                14      161	160
                14	502	502
                23	376	375
                23	447	449
                28	419	421
   a) Column 1, primer name without F or R termination;
   b) Column 2, the RT PCR product for matching to the transcriptome product(s);
   c) Column 3 (4,5...), the transcriptome product(s) for matching. This is usually just one product, but multiple products can be added, each in a seperate column.

The possible transcriptome products for each primer pair are specified in the output file with the suffix **cluster_results**. *Note that you can only use products identified for that primer pair. If you specify impossible product sizes the program will break*.

7. **output** this is the path and outfile prefix to determine name and location of outfiles, specified in line 30 of the bash script.


### Output files:

1. **primers_complete_match.txt** A list of all the primers that had product matches used in the correlation analysis
2. **\<output prefix>\_cluster_results** A file with information on each transcriptome product identified for each primer set. The file is in the following format:

                151	93	[[93, 'G30361;G30361.6(-)', ('perfect', 'perfect')], [93, 'G30361;G30361.9(-)', ('perfect', 'perfect')], [93, 'G30361;G30361.1(-)', ('perfect', 'perfect')], [93, 'G30361;G30361.5(-)', ('perfect', 'perfect')], [93, 'G30361;G30361.7(-)', ('perfect', 'perfect')]]
                151	96	[[96, 'G30361;G30361.3(-)', ('perfect', 'perfect')], [96, 'G30361;G30361.4(-)', ('perfect', 'perfect')]]
                151	122	[[122, 'G30361;G30361.2(-)', ('perfect', 'perfect')], [122, 'G30361;G30361.8(-)', ('perfect', 'perfect')]]
                75	215	[[215, 'G21499;G21499.4(-)', ('perfect', 'perfect')], [215, 'G21499;G21499.10(-)', ('perfect', 'perfect')], [215, 'G21499;G21499.6(-)', ('perfect', 'perfect')], [215, 'G21499;G21499.14(-)', ('perfect', 'perfect')], [215, 'G21499;G21499.1(-)', ('perfect', 'perfect')], [215, 'G21499;G21499.12(-)', ('perfect', 'perfect')], [215, 'G21499;G21499.7(-)', ('perfect', 'perfect')], [215, 'G21499;G21499.11(-)', ('perfect', 'perfect')], [215, 'G21499;G21499.5(-)', ('perfect', 'perfect')], [215, 'G21499;G21499.8(-)', ('perfect', 'perfect')]]
                75	220	[[220, 'G21499;G21499.9(-)', ('perfect', 'perfect')], [220, 'G21499;G21499.3(-)', ('perfect', 'perfect')], [220, 'G21499;G21499.13(-)', ('perfect', 'perfect')], [220, 'G21499;G21499.2(-)', ('perfect', 'perfect')]]


