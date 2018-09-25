##### This file contains the different functions for comparing/organising coordinates ####

# This function converts the string input of coordinates to a list of numeric values for the coordinates 
def coord_str2lst(in_str):
    coord_lst = in_str.split("|")
    coord_lst = [coord_set.split("-") for coord_set in coord_lst]
    coord_lst_sorted = sorted([sorted([int(coord) for coord in coord_set]) for coord_set in coord_lst])
    return(coord_lst_sorted)
# This function find out the overlap between two sets of coordinates
# a is coordinates for read, b is coordinates for exon
def overlap(a,b):  
    overlap_len=None    
    if a[0] <= b[0] <= a[1] :
        overlap_len=b[1]-b[0]+1 if a[1]>=b[1] else a[1]-b[0]+1
    elif b[0] <= a[0] <= b[1]:
        overlap_len=a[1]-a[0]+1 if b[1]>=a[1] else b[1]-a[0]+1
    return(overlap_len)

# check if the value is a float
def is_number(n):
    try:
        float(n)   # Type-casting the string to `float`.
                   # If string is not a valid `float`, 
                   # it'll raise `ValueError` exception
    except ValueError:
        return False
    return True