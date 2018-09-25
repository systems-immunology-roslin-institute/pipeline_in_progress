#### THIS CLASS ORGANISES THE SAM INFORMATION #####
# It contains the details on how sam file input should be read 
# .sam file format are as follows:
# [0] : Read ID (forward and reverse have the same ID)
# [1] : FLAG (check Picard SAM flags for detail https://broadinstitute.github.io/picard/explain-flags.html)
# [2] : Reference sequence NAME
# [3] : 1-based leftmost position
# [4] : Mapping quality - this score system varies drastically between different alignment programs and is therefore not used
# [5] : CIGAR string
# [6] : Ref.name of the next mate/read
# [7] : Position of the mate/next read (i.e. column 3 for the paired mate)
# [8] : Observed tempate length (this length is combined with the paired mate)
# [9] : Segment seqence
# [10] : ASCII of Phred-scaled base quality +33
# [11] : Optional information  (not used, but described here in case of future use)
# Other additional inforamtion about .sam and .gtf/.gff3 files
# Depending on the gene, some are on forward (+) strand, some on reverse (-) [within .gtf/.gff3 file]
# Similarily, in data from paired read .bam files, each read is treated as an indiviual identity. 
# If it is a reverse strand, the chromosome position would need to be reversed
# (i.e. the sequence goes from higher coordinate to lower coordinate)
# When treating individual read (forward and reversse separately), the actual coverage by the read is calculated by using the QC alignment
# The CIGAR score format is as follow:
# M : match or mismatch (can be a match or mismatch. The nucleotide is present in the reference)
# I : insertion (extra in the read)
# D : deletion (missing in th read)
# N : skipped bases on reference (skipped region; a region of nucleotides that is not present in the read) (not to be confused with when the sequence contains N)
# S : soft clipping (the clipped nucleotids are present in the read) - soft-clipped: bases in 5' and 3' of the read are NOT part of the alignment.
# H : hard clipping (the clipped nucleotides are not present in the read) - hard-clipped: bases in 5' and 3' of the read are NOT part of the alignment AND those bases have been removed from the read sequence in the BAM file. 
#     The 'real' sequence length would be length(SEQ)+ count-of-hard-clipped-bases
# P : padding (padded area in the read and not in the reference)
# X : read mismatch (the nucleotide is present in the reference)
# = : read match (the nucleotide is present in the reference)
import re
class read_sam(object):
    def __init__(self, sam_line):
        self.cigar_score_cut = ('S')
        self.cigar_score_in_read = ('M', 'I', 'S', 'X', '=', 'EQ') 
        self.cigar_score_in_read_only = ('I', 'S', 'P') 
        sam_line = sam_line.rstrip("\n\r")
        sam_values = sam_line.split('\t')
        self.readid=sam_values[0]
        self.flag=int(sam_values[1])
        self.cigar=sam_values[5]
        self.seq=sam_values[9]
        self.readid_out=None
        self.cigar_walk() if 'S' in self.cigar else 0 # clean up the sequence if there has been soft clipping
    # This function uses cigar string to split the matching segments based on the cigar score
    # Cutting off bases with a cigar score of S, which are bases present in the read sequence but not used for alignment 
    # This function is only used for when blast+ is required downstream
    def cigar_walk(self):
        cigar_split =  re.findall('\d*\D+', self.cigar)
        cigar_split = [re.findall('\d+|\D+', score) for score in cigar_split]
        read_walk = 0
        trimmed_seq = ''
        for score_pair in cigar_split:
            if score_pair[1] in self.cigar_score_in_read:
                read_walk_start = read_walk
                read_walk = read_walk+ int(score_pair[0]) 
                if score_pair[1] not in self.cigar_score_cut:                      
                    trimmed_seq = trimmed_seq + self.seq[read_walk_start:read_walk]
        self.seq=trimmed_seq
    # This function prints out the coordinates from the bed_unit
    def read_seq_coord(self, bed_unit):
        # Only print out the results if the readid_out (i.e. the equalivant of readid in the bed file) has been defined. 
        # The step filters out the reads that are not mapped to the strand of interest in paired end data
        # It takes in a bed unit and return the corresponding coordinates and sequences
        # check if the readid_out is in the bed unit            
        read_details={}    
        if self.readid_out and self.readid_out in bed_unit.readid2position:
            coord=bed_unit.readid2position[self.readid_out]
            read_details={'readid':self.readid_out,
                          'coord':coord,
                          'seq':self.seq}            
        return(read_details)

            
# This child object class is for single end RNAseq
# This child object class inherit from the parent object, which has all that cigar score details
class read_sam_single(read_sam):
    def __init__(self, *args, **kwargs):
        super(read_sam_single,self).__init__(*args, **kwargs)
        self.flag_pos=(0) 
        self.flag_neg=(16)
        self.readid_out=self.readid


# This child object class is for paired end RNAseq
# This child object class inherit from the parent object, which has all that cigar score details
# Has the additional input of strand
class read_sam_paired(read_sam):
    def __init__(self, sam_line, strand, *args, **kwargs):
        self.cigar_score_cut = ('S')
        self.cigar_score_in_read = ('M', 'I', 'S', 'X', '=', 'EQ') 
        self.cigar_score_in_read_only = ('I', 'S', 'P') 
        sam_line = sam_line.rstrip("\n\r")
        sam_values = sam_line.split('\t')
        self.readid=sam_values[0]
        self.flag=int(sam_values[1])
        self.cigar=sam_values[5]
        self.seq=sam_values[9]
        self.cigar_walk() if 'S' in self.cigar else 0 # clean up the sequence if there has been soft clipping
        self.flag_pair1 = (73, 89, 121, 101, 117, 69, 99, 83, 131, 179, 81, 97, 65, 113) # 131, 179 are for pair1 with wrong orientation
        self.flag_pair2 = (133, 165, 181, 153, 185, 137, 147, 163, 67, 115, 161, 145, 129, 177) # 67 and 115 are pair2 with wrong orientation
        self.flag_pos = (73, 89, 99, 97, 0, 137, 163, 161, 65) # 73, 89, 99, 97, 0 for pair/1, and 137, 163, 161 are for pair/2
        self.flag_neg = (121, 83, 81, 16, 153, 185, 147, 145, 129) # 121, 83, 81, 16 for pair/1, and 153, 185, 147, 145 are for pair/2
        self.strand=strand
        self.readid_out=None
        self.paired_readid()
    # This function is specifically designed for pair end data
    # To define the correspnding read_id in the bed file
    def paired_readid(self):
        record=False
        if self.strand > 0: # the gene is on the forward strand
            record = True if self.flag in self.flag_pair1 else False
        elif self.strand < 0:
            record = True if self.flag in self.flag_pair2 else False
        if record:
            if self.flag in self.flag_pos:
                pair_num = '/1'
                self.readid_out = self.readid + pair_num
            elif self.flag in self.flag_neg:
                pair_num = '/2'
                self.readid_out = self.readid + pair_num
            else:
                pair_num = ""
        
    