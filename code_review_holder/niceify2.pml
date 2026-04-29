#Reference: https://www.youtube.com/watch?v=Ztz0QCX1oF8
#set background color
bg_color white

#make selections and objects
sele lig, organic
select bb, bb.
#select sc, sc. + n. ca + PRO/N

#clean and show
remove hydrogen + solvent
hide all 
as cartoon
#show sticks, sc
show sticks, lig

#color by chain
util.cbc 
color black, lig
#color carbon by non-carbon colors
util.cnc 
#orient around backbone
orient bb 
select none


#get_view to get particular orientation.
#ray 2400,2400
#set hash_max, 200
#reinitialize