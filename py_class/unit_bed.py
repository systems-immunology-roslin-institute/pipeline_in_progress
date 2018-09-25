#### THIS CLASS ORGANISES THE BEDFILE INFORMATION #####
## .bed file format are as follows:
# note that bed format starts are zero based and ends are one based, meaning that it is easier for calculating coordinates (without needing to keep on subtracting 1 when calculating ranges)
# However, this means that for the actual coordinates, one needs to be added to the start position 
# [0] : chromosome
# [1] : coordinate start
# [2] : coordinate end
# [3] : read_id/1 or /2 for paired reads. These reads are treated as separate instances
# [4] : score
# [5] : strand
# [6] : coordinate start
# [7] : coordinate end
# [8] :
# [9] : number of blocks
# [10] : size of each block 
# [11] : starting position of each block

import re

class bed_unit(object):
    def __init__(self, bed_fp):
        self.bed_fp=bed_fp
        self.readid2position = {}
        self.read_bed()
    # This function is used to read a line from a bed file and store it in self.readid2position
    def read_bed(self):
        with open(self.bed_fp, 'r') as open_bed:
            for line in open_bed:
                line = line.rstrip("\n\r")
                values = re.split('\t', line)
                read_name = values[3]
                block_size = [int(x) for x in values[10].split(',')]   
                block_start = [int(x) for x in values[11].split(',')]
                coord_start = int(values[1])
                block_coord = []
                for block_id in range(0, len(block_size)):
                    current_block_coord=[coord_start + block_start[block_id] , coord_start + block_start[block_id]+block_size[block_id]] # bedfile start position is 0 based, and end position is 1 based to allows for simplier calculation for length, so for the true start position, one needs to be added to the start position only
                    current_block_coord='-'.join([str(x) for x in current_block_coord])
                    block_coord.append(current_block_coord)
                read_address='|'.join(block_coord) if len(block_coord)>1 else block_coord[0]
                if read_name in self.readid2position:
                    if read_address not in self.readid2position[read_name]:
                        self.readid2position[read_name].append(read_address)
                else:
                    self.readid2position[read_name]=[read_address]
    
    
