reinitialize
fetch 5zr1

run C:/Users/Luis/Dropbox (MIT)/Lab/Projects/automate-the-boring-stuff/script/niceify2.pml

set cartoon_side_chain_helper, on

sele o14, (c. A and (resi 485 or resi 567)) or (c. D and (resi 267)) or (lig and c. A)
center o14
orient o14
show sticks, o14

distance Rlig, /5zr1/D/D/ARG`267/CZ, /lig/I/A/AGS`2001/PG
distance Elig, /5zr1/A/A/GLU`567/CD, /lig/I/A/AGS`2001/PG
distance Klig, /5zr1/A/A/LYS`485/NZ, /lig/I/A/AGS`2001/PG
# hide just the distance labels
hide labels
set dash_gap, 0.5
set dash_radius, 0.1


sele intheway, (c. A and resi 478,479,481,482,695,696,699,694,698,697,700,701,703,704,705,483,628,627,626,625)
hide everything, intheway
set cartoon_transparency, 0.5, 5zr1

set_view (\
     0.653128624,   -0.188511506,    0.733406842,\
    -0.033014823,    0.960509539,    0.276285082,\
    -0.756525755,   -0.204663277,    0.621111512,\
     0.000000000,    0.000000000,  -63.208839417,\
   146.743408203,  113.786521912,   90.849098206,\
    49.977756500,   76.439926147,  -20.000000000)

set ray_trace_fog,0
set ray_shadows,0
unset depth_cue
bg_color white
set antialias,2
set hash_max, 300

ray 2400, 2400
png C:/Users/Luis/Dropbox (MIT)/Lab/Projects/automate-the-boring-stuff/output/ORC1_4_AAA+.png, 0, 0, -1, ray=0

#convert  -compress lzw  raytraced.png raytraced.tiff
