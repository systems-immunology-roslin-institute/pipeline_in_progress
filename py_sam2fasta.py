#### SCRIPT DESCRIPTION ####
# script description needs rewriting
####

import sys
import os
#import numpy as np
import py_class
from py_class.unit_gtf import gtf_info
from py_class.unit_sam import read_sam
from py_class.unit_reads import read_unification
from py_class.unit_reads import read_unification_no_blast
from py_class.unit_sample_annotation import sample_annotation
from py_functions import coord_str2lst

in_fp_cache =  sys.argv[1] # cache directory
in_fp_output = sys.argv[2] # output folder
in_unification =  sys.argv[3]
in_gene_name = sys.argv[4]
in_strand = sys.argv[5]
in_coord_range_interval = sys.argv[6]
in_sample_annotation = sys.argv[7]
in_number_pair2 = int(sys.argv[8])

#### STEP 0 : Use the .sys input to determine the fp for required files and outputs
in_fp_gtf = in_fp_cache + '/' + in_gene_name + '/gtf/' + in_gene_name + '.gtf'
in_fp_sam = in_fp_cache + '/' + in_gene_name + '/sam'
in_fp_bed = in_fp_cache + '/' + in_gene_name + '/bed'
out_fp_gml = in_fp_cache + '/' + in_gene_name + '/gml'
in_fp_merged_gtf = in_fp_cache + '/' + in_gene_name + '/gtf/' + in_gene_name + '_merged.gtf'
if not os.path.exists(out_fp_gml):
    os.makedirs(out_fp_gml)
# other input setting from sys.argv
chromosome=in_coord_range_interval.split(":")[0]
paired_read = True if in_number_pair2 > 0 else False


#### STEP 1: read in and organise exon annotation from the .gtf/.gff3 file ####
# Characterise exon definitions for transcript models in the gtf file
gtf_file_extension = (in_fp_gtf.split('.')[-1]).lower()
if gtf_file_extension not in ['gtf', 'gff3']:
    print("Error: the input gtf/gff3 has incoorect extension. Please check the format.")
gtf_organised = gtf_info(strand = in_strand)
with open(in_fp_gtf, "r") as gtf:
    for line in gtf:
        gtf_organised.read_gtf(line)
gtf_organised.position_annotation_exon() # reorganise the dictionary so it is keyed by chromosome loci

#### STEP 2 : organise the sam data ####
sample_annotation_obj=sample_annotation(sam_fp=in_fp_sam, bed_fp=in_fp_bed,sample_annotation_fp=in_sample_annotation)
if in_unification == 'no_blast':
    read_summary=read_unification_no_blast(sample_annotation_obj=sample_annotation_obj, genename=in_gene_name, strand=in_strand, chromosome=chromosome, paired_read=paired_read, unification=in_unification)
else:
    read_summary=read_unification(sample_annotation_obj=sample_annotation_obj, genename=in_gene_name, strand=in_strand, chromosome=chromosome, paired_read=paired_read, unification=in_unification)
read_summary.organise_readid_info()
gtf_organised.add_read_defined_segments(read_summary.read_defined_segments(gtf_organised)) # add read_defined gene annotation into the gtf object

# use the junctions in the read summary to remake the all transcript models in the gtf file that contains the string "merged_gtf" in their name
# Then, annotate the reads
read_summary.print_output(gtf_obj=gtf_organised, cache_dir = in_fp_cache, out_dir=in_fp_output)


