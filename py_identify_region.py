#### SCRIPT DESCRIPTION ####
# This script identify the chromosomal region of the gene/genes of interest
# It saves the gtf file for each gene into a separate file for later use

#### INPUT / OUPUT ####
# INPUT (from .sh):
# - gene annotation (.gtf/.gff) filepath (this script has only been optimalised to recognise .gtf or .gff3 file format)
# - cache folder filepath
# - gene of interest
# OUTPUT (to .sh):
# - chromomsome regions for the gene of interest
# OUTPUT (to file):
# - gene specific gtf (under cache folder)
# - chromosome regions for the gene of interest

#### INPUT BACKGROUND ####
# GTF/GFF file format (downloadable from Ensembl. Needs to be the one used for sorting the .bam file)
# Although the gene annotation file can be .gtf or gff, throughout the script, it is referred to as "GTF"
# Both gtf and gff files are tabl delimited files with the attributes in the following columns:
#[0] seqname : name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
#[1] source : name of the program that generated this feature, or the data source (database or project name)
#[2] feature : feature type name, e.g. Gene, Variation, Similarity
#[3] start : Start position of the feature, with sequence numbering starting at 1.
#[4] end : End position of the feature, with sequence numbering starting at 1.
#[5] score : A floating point value.
#[6] strand : defined as + (forward) or - (reverse).
#[7] frame : One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
#[8] attribute : A semicolon-separated list of tag-value pairs, providing additional information about each feature.
# The script has only been optimalised to understand .gtf gene/exon annotation formats

#### STEP 0 : SETTING ####
setting_warning_message_multichromosome = "Gene found on multiple chromosomes"
setting_warning_message_gene_not_found = "Gene not found in supplied gtf/gff file"
setting_warning_message_scripterror = "Error when identifying chromosome region for the gene of interest"
setting_col_chr = 0
setting_col_chr_start = 3
setting_col_chr_end = 4
setting_col_strand = 6


#### Import library ####
import sys, re, os, datetime, subprocess
import numpy as np

#### Read in inputs ####
in_fp_gtf = sys.argv[1] # input from shell script #'input_gtf/mm10_chr19-1-20000000_Ensembl.gtf' #
in_fp_cache = sys.argv[2] # input from shell script # 'cache' #
in_gene = sys.argv[3] # input from shell script #'IGHMBP2' 
in_gene = in_gene[0:(len(in_gene)-1)].split('|')

# find out the GTF file name and extension. This is for identifing if it is a .gtf or .gff file, and for naming the cache-gene-specific gtf/gff file
if re.match(r'^(.*?)/(.*?)', in_fp_gtf, re.IGNORECASE):
    # if the input gtf filepath is in a subdirectory, strip away the folder names and retain only the file name
    gtf_filename_extension = in_fp_gtf.split("/")[-1]
else:
    gtf_filename_extension = in_fp_gtf
gtf_extension = (gtf_filename_extension.split("."))[-1].lower()

#### Organise the outputs ####
# Create the output filepaths and check if the subfolders exists. Create them if not.
out_region_fp = in_fp_cache + '/chromosome_region.txt'
out_region_warning_fp = in_fp_cache + '/chromosome_region_error.txt'
#### STEP 1: Go through the gtf file ####
# Identify the starting/ending position of the gene of interest, and wether if the gene is present in multiple chromosomes
# and print out the gtf in cache for the gene(s) of interest for use in later steps (annotating the exon ID)
# Do these steps at a per-gene basis
merged_gene_name = "|".join(in_gene)
match_regex= r'^.*?gene_name \"(%s)\";.*?$' %(merged_gene_name) # Find lines in the gtf files describing the the gene of interest
#match_regex= r'^.*?gene_id \"(%s)\";.*?$' %(merged_gene_name) # Find lines in the gtf files describing the the gene of interest
if not os.path.exists(os.path.dirname(out_region_fp)):
    os.makedirs(os.path.dirname(out_region_fp))
out_region_warning = open(out_region_warning_fp, 'w+')
in_gtf = open(in_fp_gtf, "r") # read in input gtf
all_gene_info = {}
previous_gene=None
for line in in_gtf:
    line = line.rstrip("\r\n")
    # only regex the line if it is labelled as gene
    matched_values=re.match(match_regex, line, re.IGNORECASE)
    if matched_values:
        c_gene = matched_values.group(1).upper()
        # If this gene is different from the previous gene, close the file for the previous gene, unless it's None
        if previous_gene != c_gene:
            out_gtf.close() if previous_gene else 0 # close the file
            # Open the gtf file for writing for the current gene
            out_gtf_fp = in_fp_cache + '/' + c_gene +'/gtf/' + c_gene + '.' + gtf_extension # label the gene-specific gtf with the original file name and gene name
            if not os.path.exists(os.path.dirname(out_gtf_fp)): # make the file/folder if it does not exist
                os.makedirs(os.path.dirname(out_gtf_fp))
            out_gtf = open(out_gtf_fp, 'a') # open gtf file location for appending the matched line
            previous_gene = c_gene
        out_gtf.write("%s\n" %(line))
        # column 0 is chromosome id, column 3 is chromosome start, column 4 is chromosome end
        values = re.split('\t', line)
        # If this line defines gene, store the gene information
        c_chr = str(values[setting_col_chr])
        c_start = int(values[setting_col_chr_start])
        c_end = int(values[setting_col_chr_end])
        c_strand =values[setting_col_strand]
        if c_gene not in all_gene_info:
            all_gene_info[c_gene]={}
            all_gene_info[c_gene]['chr_num'] = c_chr
            all_gene_info[c_gene]['chr_start'] = c_start
            all_gene_info[c_gene]['chr_end'] = c_end
            all_gene_info[c_gene]['strand'] = c_strand
        else:
            if (c_chr != all_gene_info[c_gene]['chr_num']) or (c_strand != all_gene_info[c_gene]['strand']):
                out_region_warning.write("%s found on multiple chromosomes "  %(c_gene) + ':'.join([c_chr + c_strand, all_gene_info[c_gene]['chr_num'] + all_gene_info[c_gene]['strand']]))
                print(setting_warning_message_multichromosome)
            else:
                if all_gene_info[c_gene]['chr_start'] > c_start:
                    all_gene_info[c_gene]['chr_start']=c_start
                if all_gene_info[c_gene]['chr_end'] < c_end:
                    all_gene_info[c_gene]['chr_end']=c_end 

in_gtf.close()
out_gtf.close()


#### STEP 2: Output ####
# Write the chromosome region into the cache folder. This is not used for later, but more for trouble shooting in case of 1) gene not found or 2) gene found in multiple chromosome locations

# write down the chromosome region
out_region =open(out_region_fp, 'w+')
bed_region=os.path.join(in_fp_cache, "bed_chr_regions.bed")
open_bed_region=open(bed_region,'w+')
gene_ranges={}
gene_names={}
# if the region is unique in 1 chromosomal location, and the forward/reverse strand were found
# save the infomration in a .txt file, which would be used as input for .sh later
# Make a bed file from the chromosome regions for the genes of interest so it can be intersected with merged gtf
for gene in in_gene:
    if gene not in all_gene_info:
        out_region_warning.write("%s not found : " %(gene) + setting_warning_message_gene_not_found)
    else:
        gene_info = all_gene_info[gene]
        gene_chr = gene_info['chr_num'] 
        gene_start = gene_info['chr_start']
        gene_end = gene_info['chr_end']
        gene_strand = gene_info['strand']
        out_region_str = "%s:%d-%d" %(gene_chr,gene_start,gene_end)
        out_strand_str="f" if gene_info['strand'] == '+' else "r" # the input for samtools is -F for forward strand and -f for reverse strand   
        out_region.write('%s,' %(gene))
        out_region.write(out_region_str) # chromosome location. Input for samtool/sambamba
        out_region.write(',%s\n' %(out_strand_str)) # strand. Input for samtool/sambamba
        bed_outline = "%s\t%s\t%s\t%s\n"%(gene_chr, gene_start, gene_end, gene)        
        open_bed_region.write(bed_outline)
        gene_ranges.update({gene_chr:[]}) if gene_chr not in gene_ranges else 0
        gene_names.update({gene_chr:[]}) if gene_chr not in gene_names else 0
        gene_ranges[gene_chr].append([gene_start,gene_end])
        gene_names[gene_chr].append(gene)
out_region.close()
out_region_warning.close()
open_bed_region.close()

