def coord_str2lst(self,in_str):
    coord_lst = in_str.split("|")
    coord_lst = [coord_set.split("-") for coord_set in coord_lst]
    coord_lst_sorted = sorted([sorted([int(coord) for coord in coord_set]) for coord_set in coord_lst])
    return(coord_lst_sorted)
