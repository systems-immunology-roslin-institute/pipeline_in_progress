# NGS-Graph-Generator2
## Pre-requisites
### Required tools:
- Samtools (version 1.4)
- R (version 3.2.5)
- Sambamba (version 0.6.6)
- BEDTools (version 2.25.0)
- Blast+ (version 2.6.0)

Please install/add the above tools prior to running this pipeline. The tool versions used in the developing this pipeline are indicated in the bracket. Please modify line 2 to 9 in generate_network.sh to the paths/commands for these tools on your system. Once the modifications are done, run generate_network.sh to create alternative-splicing network graphs.

## Input options
### Required:
- b	: 	Folder with BAM files 
- g	: 	GTF file (matching the genome build for bam alignment)
- n 	: 	A text file with gene names in each line

### Optional:
- o	: 	The output folder (default “output”)
- c	: 	The cache directory (default “cache”)
- p	: 	The percentage similarity value (default 85) 
- l	: 	The percentage coverage value (default 30)
- u	:	The method used to unify reads. Options are:
       by_seq (by 100% sequence identity) 
       by_coord (by same start and end positions) (default)
       no_blast (by binning the chromosome positions into specific intervals, require the -r setting)
-a 	:	If sample annotation is supplied, node/edge are averaged within specified groups 
-s 	:	Folder with sorted input bam files. 

> Example 1:
> bash generate_network.sh \
> -b input_bam \
> -n gene_list.txt \
> -g input_gtf/Mus_musculus.GRCm38.91.chr.gtf 

> Example 2:
> bash generate_network.sh \
> -b input_bam \
> -n gene_list.txt \
> -g input_gtf/Mus_musculus.GRCm38.91.chr.gtf \
> -a sample_annotation.txt \
> -o output \
> -u no_blast \
> -s input_bam
