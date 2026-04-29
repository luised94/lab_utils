reinitialize
fetch 5BK4

run C:/Users/Luis/Dropbox (MIT)/Lab/Projects/automate-the-boring-stuff/script/niceify2.pml

set cartoon_side_chain_helper, on
color orange, c. I and 7MCA
util.cnc 

alter c. S, resi=int(resi)+61
sele aaafolds, (c. A and (resi 493-700)) or (c. B and (resi 359-566)) or (c. C and (resi 518-725)) or (c. D and (resi 366-573)) or (c. E and (resi 525-732)) or (c. F and (resi 410-617)) 
sele dnainchannel, (c. S and (resi 46-60)) or (c. O and (resi 1-16))

hide everything
show cartoon, aaafolds
show cartoon, dnainchannel
show sticks, lig

#Use the following lines to get the view below. 
#centerofmass aaafolds
#origin position=[ 139.621, 209.389, 139.520]
#center origin

set_view (\
    -0.761089265,    0.639125884,    0.110753059,\
    -0.066383310,   -0.246594980,    0.966842890,\
     0.645244181,    0.728501499,    0.230107248,\
     0.000000000,    0.000000000, -331.004669189,\
   139.611953735,  209.416946411,  139.502182007,\
   260.966613770,  401.042724609,  -20.000000000 )

set ray_trace_fog,0
set ray_shadows,0
unset depth_cue
bg_color white
set antialias,2
set hash_max, 300

ray 2400, 2400
png C:/Users/Luis/Dropbox (MIT)/Lab/Projects/automate-the-boring-stuff/output/MCM_AAA+folds.png, 0, 0, -1, ray=0

#convert  -compress lzw  raytraced.png raytraced.tiff
