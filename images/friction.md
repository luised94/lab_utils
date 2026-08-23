# Recording thoughts and friction encountered during gel shift imaging
date: 2026-08-22

Annoying but not a big deal. Open files in imagej. I think rotating and saving does duplicate data unnecessarily but not super important. Solution could be to record how to rotate the original.
    > Requires more data input. Wonder if it can be back calculated from the rotated image.
I can assign the lanes while rotating right? What is better? Close and reopen, make sure image is fine? Use assembly line method for each step.

## 2026-08-23
Did not do all of the gels. Went home and came back.
Have a loading gel that I have to analyze.
Got distracted trying to back up dropbox. Just noting the context switch.
Rotated the image a bit. nothing annoying there. Wonder if there is a more sophisticated and mathematical real time way to see if there is an optimal rotation angle (0 to 180 i suppose?) I know the lane profiles can be collapsed along the two axis (loadign and migration axis). At some point those profiles should reach a maximum of signal and distinctiveness right? Probably use as realtime feedback in imagej would be cool to handle that but python also possible and then it applies the transformation. Lane definitions would still be manual. Probably deferred idea.
Holy moly. the gel is so clean this time that I can see all of the lanes. May have to adjust that. Could have been a wrong assumption all along.Best I can do is to space evenly or something.
Trying the distribute selection for the first time. Kind of annoying that I have to find the file in wsl. Not sure what to do about that now. Probably low prio.
Some of the lanes were smiling a bit and if they have to much protein, I can draw the lanes cleanly. Consider best practices for salvaging in imagej and mostly a experimental issue. Right now, I try to fit as tight as possible to the hardest one. If I can decide band versus region easily then this may not be a big issue. I wonder if instead of regions, I just draw multiple horizontal region for each band I see. Region is then that. Risky for tight clusters.
Messed up the distribute selection immediately. I think the current options is too big. They went out of the image bounds and the rest got drawn to the center. Should I just skip? Tried one more time. Unchecked the box to apply in pixels. Some of the options dont denote units. The script also duplicates the current selection as opposed to just drawing the rectangle. I was not clear enough that the rectangle shape and region is what is duplicated. woops. Decide to continue manually. Would be nice to fix but its not a big deal at this moment.
Distracted by archivign dropbox. Important thing I have to get done .
Cant see lanes when pressing t. Have to click show all in ROI manager. Ctrl + 1 is unneccessary. Drew 14 lanes instead of 15. Did it twice. With this gel its very hard to figure out how many lanes there are (this is actually optimal. Very low background.) May have to reduce the guard.
Yup. the export script blocks me. Had to choose other so pretty ok. Maybe a better clearer message.
All outputs good for this gel. Continue with defining the lanes in the other gels.
Had to click a lot to get to the imagej script but runs even though its on wsl side.
Good. Did one of the gel shift experiments. Lanes have to be balanced between thin lanes and bulckier lanes.
Can copy paste the imagej script and run using Plugins > New > Macro and pasting the file into the box then hit run instead of Run and go through all of the clicks and folders. After first gel, can leave the macro window open, or just open at the beginning.
Huh. If i deselect the window, something happens to the selections. Basically need to start over. They dont stay in place. Not sure how to describe Had to delete, reselect and then start over.
Ok. Did all three. The process is not too terrible to be honest. Moving the rectangle selection is kind of slow. if that was faster or tunable, that would be the best.
At this time, I think I need to export sheets but I did forget how to use session template or if the scripts can be done next. Let me read PROTOCOL.md (actually HOWTO.md forgot, also press gf on top and I go to file. Very nice.)
Yup. Can use gf to navigate. Very nice. Thanks nvim.
 Made some modifcations to HOWTO.md with markers of edit. Maybe I should delete them or mark but then get the diff to make sure I know. Maybe there is a better protocol or procedure for annotations/tracking changes/git for writing and editing in nvim.
The sample sheet has 10 wells but I exported nine. Will that cause an error? Think the error should be removed. If I have option to select a different number of lanes, why would I block it?
Got kind of hungry. Left to eat something.
Came back after eating. Had quick beer. Kinda late now (7pm). Taking me a while. Maybe its distractions. I shall see.
Now is the export sheets. annoying if the multiple gels are in different folders but not too bad.
Forgot to fix some of the labels but editing in excel is good.
Encountered the windows path character limit while exporting sample sheets.
Should I incldue a small csv export automator in the session template? Error during export. file cant be accessed. Thought it was the path but still getting the pop up window. Maybe I need to fully restart.
Had to copy the manual profiles csv ctrl a delete and then paste the sample sheet.
Got to the session run. Not sure what regions to pick and it doesnt work for gels from different experiments. Damn.
Went back to the serve gel picker (should be in the howto. think its purpose was mistaken at some point.)
Had to go to the directory in wsl and use the printf line.
Error using serve gel picker. The parent directory gets set by tif file so the script looks in that parent directory when the file is in the <stem>_gel_analysis directory. I think files should have the gel as point of entry. Everything is then derived from that.
