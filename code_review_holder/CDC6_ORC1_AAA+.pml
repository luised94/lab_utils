reinitialize
fetch 7MCA

run C:/Users/Luis/Dropbox (MIT)/Lab/Projects/automate-the-boring-stuff/script/niceify2.pml

set cartoon_side_chain_helper, on
color orange, c. I and 7MCA
util.cnc 

sele c6o1, (c. A and (resi 616)) or (c. I and (resi 114 or resi 224)) or (lig and c. I)
center c6o1
orient c6o1
show sticks, c6o1


#sele intheway, (c. A and resi 478,479,481,482,695,696,699,694,698,697,700,701,703,704,705,483,628,627,626,625)
#hide everything, intheway
set cartoon_transparency, 0.3, 7MCA

set_view (\
     0.329257607,   -0.103709482,   -0.938521206,\
     0.927203953,    0.223447070,    0.300590754,\
     0.178538367,   -0.969175398,    0.169733956,\
     0.000000000,    0.000000000,  -53.428050995,\
   127.544395447,  152.950195312,  131.362655640,\
    42.123081207,   64.733016968,  -20.000000000 )

set ray_trace_fog,0
set ray_shadows,0
unset depth_cue
bg_color white
set antialias,2
set hash_max, 300

ray 2400, 2400
png C:/Users/Luis/Dropbox (MIT)/Lab/Projects/automate-the-boring-stuff/output/CDC6_ORC1_AAA+.png, 0, 0, -1, ray=0

#convert  -compress lzw  raytraced.png raytraced.tiff
