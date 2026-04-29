reinitialize
fetch 7MCA

run C:/Users/Luis/Dropbox (MIT)/Lab/Projects/automate-the-boring-stuff/script/niceify2.pml

set cartoon_side_chain_helper, on
color orange, c. I and 7MCA
util.cnc 

sele aaafolds, (c. A and (resi 471-628)) or (c. B and (resi 300-489)) or (c. C and (resi 92-273)) or (c. D and (resi 66-278)) or (c. E and (resi 31-168)) or (c. I and (resi 102-263)) 
sele dnainchannel, (c. G and (resi 10-31)) or (c. H and (resi 56-76))
hide everything
show cartoon, aaafolds
show cartoon, dnainchannel
show sticks, lig

#Use the following lines to get the view below. 
#centerofmass aaafolds
#origin position=[ 107.580, 120.431, 127.586]
#center origin

set_view (\
    -0.527550519,   -0.201993212,    0.825157940,\
    -0.841823280,    0.254780442,   -0.475835592,\
    -0.114118114,   -0.945666373,   -0.304451793,\
    -0.000000000,    0.000000000, -323.199737549,\
   107.580001831,  120.430999756,  127.585998535,\
   254.813293457,  391.586181641,  -20.000000000 )

set ray_trace_fog,0
set ray_shadows,0
unset depth_cue
bg_color white
set antialias,2
set hash_max, 300

ray 2400, 2400
png C:/Users/Luis/Dropbox (MIT)/Lab/Projects/automate-the-boring-stuff/output/AAA+folds.png, 0, 0, -1, ray=0


sele c. G or c. H
ray 2400, 2400
png C:/Users/Luis/Dropbox (MIT)/Lab/Projects/automate-the-boring-stuff/output/AAA+folds.png, 0, 0, -1, ray=0

#convert  -compress lzw  raytraced.png raytraced.tiff
