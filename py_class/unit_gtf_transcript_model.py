#### This class object stores information for transcript models ####
# Once the full transcript model has been added in , .make_exon_dict and .make_other_dict are used to summarise all info
# .make_exon_dict and .make_other_dict are a dictionary of coordinates, and the exon annotation for these coordinates (for instance 12323:'exon3', 12324:'stop_codon')
import operator
class transcript_model(object):
    def __init__(self, other_overwrite_order = None):
        self.exons = []
        self.other = {}
        self.exon_dict = {}
        self.other_dict = {}
        self.merged_dict = {}
        self.strand=None # this strand information is stored as -1 or 1 for reverse and forward
        self.non_overlapping_threshold=3 # As some of the starts/ends of exons may be fuzzy, I have enforced a required threshold for an intron/exon junction to be defined as separate exons (see self.non_overlapping_coordset for details)
        # In order to minimalise the number of "tracks" in annotation there are certain orders that permits overwriteing of annotation. 
        # For instance, "start codon" overwrites "CDS" (protein coding sequence). The smaller numbers indicates higher priority
        # Ignore all other annotation that are not listed here (such as "gene", "transcript", "Selenocysteine")
        if other_overwrite_order is None:
            self.other_overwrite_order={'start_codon':1,
                                        'stop_codon':2,
                                        'five_prime_utr':3,
                                        '5utr':3,
                                        'three_prime_utr':4,
                                        '3utr':4,
                                        'utr':5,
                                        'cds':6,
                                        'inter':7, #  intergenic region, one which is by almost all accounts not transcribed
                                        'inter_cns':8, #intergenic conserved noncoding sequence region
                                        'intron_cns':9, #conserved noncoding sequence region within an intron of a transcript
                                        'selenocysteine':10}
    def add_annotation(self, transcript_info):
        if not self.strand: # find out the strandedness of the current transcript model
            self.strand=-1 if transcript_info['strand'] == '-' else 1 
        if transcript_info['annotation_type'] == "exon":
            self.exons.append(transcript_info['coord'])
        # only record pre-set values
        elif transcript_info['annotation_type'] in self.other_overwrite_order:
            if transcript_info['annotation_type'] not in self.other:
                self.other.update({transcript_info['annotation_type']:[transcript_info['coord']]})
            else:
                self.other[transcript_info['annotation_type']].append(transcript_info['coord'])
        elif transcript_info['annotation_type'] not in self.other_overwrite_order:
            pass
    # Sort the exons so the returned exon numbers would be labelled sequentially from 5' to 3'
    # Create a searchable hash for each bp with an exon annotation
    def make_exon_dict(self):
        # As there may be overlapping exon annotations (due to merged file)
        # each coordinate range is walked through, with ends recorded for breakage
        if len(self.exons) > 0:
            coord_list = self.non_overlapping_coordset(self.exons)
            coord_list = self.exons
            if self.strand>0:
                coord_list.sort()
            elif self.strand<0:
                coord_list.sort(reverse=True)
            else:
                coord_list.sort()
            for exon_id, coord_set in enumerate(coord_list):
                current_coord_set = sorted(coord_set)
                for coord in range(current_coord_set[0], current_coord_set[1] + 1):
                    self.exon_dict[coord] = "e%02d"%(exon_id+1) # plus one because the counting in enumerate is 0-index'd

        self.exons = [] # set this to null as it is no longer required
    # Create a serachable hash for each bp with an "other" annotation
    def make_other_dict(self):
        if len(self.other) > 0:
            for annotation_type, coord_set_lst in self.other.items():
                current_annot_rank = self.other_overwrite_order[annotation_type]
                for coord_set in coord_set_lst:               
                    for coord in range(coord_set[0], coord_set[1]+1):
                        if coord in self.other_dict:
                            self.other_dict[coord] = annotation_type if current_annot_rank < self.other_overwrite_order[self.other_dict[coord]] else self.other_dict[coord]
                        else:
                            self.other_dict[coord] = annotation_type
        self.other={} # set this to null as it is no longer required
    # This function makes a merged dictionary from self.exon_dict and self.other_dict
    # This is then used for querying against genes
    def merge_exon_other(self):
        self.merged_dict = self.exon_dict
        for coord, annotation in self.other_dict.items():
            if coord not in self.merged_dict:
                self.merged_dict[coord] = annotation
            # overwrite the exon annotation if it is any of the following: start_codon':1, 'stop_codon':2, 'five_prime_utr':3, '5utr':3, 'three_prime_utr':4, '3utr':4, 'utr':5, 
            elif self.other_overwrite_order[annotation] <= 5:
                self.merged_dict[coord] = annotation
            else:
                pass
    # This function takes in a list of coordinate list (in the form of [start,end])
    # organise it and return non-redundent segments
    def non_overlapping_coordset(self, coord_list):
        coord_list.sort()
        min_coord=coord_list[0][0]
        max_coord=coord_list[-1][-1]
        coord_covered={}
        coord_start={} # this indicate that this coordinate has been marked as the start of an exon
        coord_end={} # this indicate that this coordinate has been marked as the end of an exon
        # mark all bases covered by at least one exon range  
        for coord_set in self.exons:
            coord_set.sort()
            min_coord = coord_set[0] if coord_set[0] < min_coord else min_coord
            max_coord = coord_set[1] if coord_set[1] > max_coord else max_coord
            coord_start[coord_set[0]]=1
            coord_end[coord_set[1]]=1
            for i in range(coord_set[0], coord_set[1]+1):
                coord_covered[i]=1
        sorted_coord_list=[]
        current_chunk=[min_coord,None]
        for coord in range(min_coord+1, max_coord+1): 
            if coord in coord_covered:
                if coord in coord_end and current_chunk[0]:
                    # if the start coordinate of chunk has already been defined, append this to the list
                    if current_chunk[0]: 
                        current_chunk[1] = coord
                        sorted_coord_list.append(current_chunk)
                        current_chunk=[None,None]
                    # otherwise the current coord is the start of the current_chunk
                elif coord in coord_start and not current_chunk[0]:
                    current_chunk[0]=coord
                else:
                    pass
        return(sorted_coord_list)
    # This function allows a list of coordinates to be added in as exon_list
    def defined_by_segment_list(self, segment_list, strand):
        self.strand = strand
        if self.strand>0:
            segment_list.sort()
        elif self.strand<0:
            segment_list.sort(reverse=True)
        else:
            segment_list.sort()
        for segment_id, current_segment in enumerate(segment_list):
            current_segment = sorted(current_segment)
            for coord in range(current_segment[0], current_segment[1] + 1):
                self.exon_dict[coord] = "e%02d"%(segment_id+1) # plus one because the counting in enumerate is 0-index'd
        return(segment_list)
    