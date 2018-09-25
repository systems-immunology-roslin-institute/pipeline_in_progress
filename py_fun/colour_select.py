import hashlib
import colorsys
def graphia_colour(in_str):
    md5_obj=hashlib.md5(in_str.encode('utf-8'))
    converted=md5_obj.digest()
    hue=0
    lightness =0
    
    for i in range(0,len(converted)):
        current_byte=converted[i] + 128;
        if i < len(converted)/2:
            hue=(hue+current_byte)%255
        else:
            lightness = (lightness+current_byte)%127
    lightness += 64
    rgb_val=colorsys.hls_to_rgb(hue/360,lightness/255,255/255)
    return(rgb_val)

