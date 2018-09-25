#### THIS CLASS ORGANISES THE SAMPLE ANNOTATION INFORMATION #####

import re, os

class sample_annotation(object):
    def __init__(self, sam_fp, bed_fp, sample_annotation_fp):
        self.sam_fp=sam_fp
        self.bed_fp=bed_fp
        self.annotation_fp=sample_annotation_fp
        self.sample_info={}
        self.group_size={}
        self.annotation_other={}
        self.annotation_other_header=[]
        self.special_char='[^\w -]'
        self.group_order_idx={}
        self.organise_annotation()
    # This function assumes that the basenames for sam, bed and first column of sample annotation would be identical (in the absence of file extension)
    def organise_annotation(self):
        sample_dict=self.get_sample_annotation()
        sam_dict={}
        bed_dict={}
        for filename in os.listdir(self.sam_fp):
            sam_dict.update({os.path.splitext(filename)[0]:os.path.join(self.sam_fp, filename)})
        for filename in os.listdir(self.bed_fp):
            bed_dict.update({os.path.splitext(filename)[0]:os.path.join(self.bed_fp, filename)})
        sample_names=list(set(sample_dict.keys()).intersection(set(sam_dict.keys())).intersection(set(bed_dict.keys())))
        sample_names.sort()
        for sample_name in sample_names:
            self.sample_info.update({sample_name:{'annotation': sample_dict[sample_name], 'sam_fp':sam_dict[sample_name], 'bed_fp':bed_dict[sample_name]}})
            if self.annotation_fp != " ":
                self.group_size.update({sample_dict[sample_name]:0}) if sample_dict[sample_name] not in self.group_size else 0
                self.group_size[sample_dict[sample_name]] +=1
            else:
                self.group_size.update({sample_name:1})
        all_group_names=list(set(self.group_size.keys()))
        all_group_names.sort()
        for i in range(0,len(all_group_names)):
            self.group_order_idx[all_group_names[i]]=i
        
    # This function get the sample annotation from the annotation file, if the file path does not equal to " "
    def get_sample_annotation(self):
        sample_group={}
        if self.annotation_fp != " ":
            with open(self.annotation_fp, 'r') as open_sample_annotation:
                for line_num, line in enumerate(open_sample_annotation):
                    line = line.rstrip()
                    line = line.split("\t")
                    if line_num == 0:
                        self.annotation_other_header = [re.sub(self.special_char, '_', val) for val in line[2:len(line)]]
                    else:
                        sample=line[0]
                        # remove special characters
                        anontation=re.sub(self.special_char, '_', line[1])# remove any special characters in the group name
                        sample_group.update({sample:anontation})
                        if len(line) > 2:
                            self.annotation_other[sample]= [re.sub(self.special_char, '_', val) for val in line[2:len(line)]]
        # If no sample annotatino was given (i.e. sample_annotation_fp == " "), use the sample name as the annotation
        else:
            sam_list = os.listdir(self.sam_fp)
            for i in range(0,len(sam_list)):
                sample=os.path.basename(sam_list[i])
                sample=sample.split(".")[0]     
                sample_group.update({sample:sample})
        return(sample_group)




