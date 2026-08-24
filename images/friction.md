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

Dealt with the serve gel picker but forgot that I have to run the other scripts before.
Had to do manually.
usb[ ]luis@Luis:~/personal_repos/lab_utils/images$ GEL_PATH="/mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.22_20.11.54_rotated_LM-0008_s0015_SDSPAGE_input-double-check-of-ORC-preps_Fl-Green.tif"
Went back to the HOWTO.md to see order. serve gel picker not there.
Had to patch the first three files of the pipeline as well. then got the error as predicted:
usb[ ]luis@Luis:~/personal_repos/lab_utils/images$ uv run validate_sample_sheet.py "$GEL_PATH"
[input] gel analysis directory: /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.22_20.11.54_rotated_LM-0008_s0015_SDSPAGE_input-double-check-of-ORC-preps_Fl-Green_gel_analysis
[input] using sample sheet: /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.22_20.11.54_rotated_LM-0008_s0015_SDSPAGE_input-double-check-of-ORC-preps_Fl-Green_gel_analysis/sample_sheet.csv
[input] using profile CSV: /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.22_20.11.54_rotated_LM-0008_s0015_SDSPAGE_input-double-check-of-ORC-preps_Fl-Green_gel_analysis/manual_lane_profiles.csv
[profile] profile has 14 lanes, comb_well_count=14
[read] sample sheet has 15 rows and columns: lane_index, well_number, lane_content, analysis_role, condition_type, sample_label, prep_source
[identity] lane-to-well relationship: undetermined
[checks] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.22_20.11.54_rotated_LM-0008_s0015_SDSPAGE_input-double-check-of-ORC-preps_Fl-Green_gel_analysis/sample_sheet_checks.json; overall fail
[checks] WARNING positive_control_present: no lane marked condition_type=positive_control
[checks] WARNING negative_control_present: no lane marked condition_type=negative_control
[checks] FAIL row_count_matches_comb: sheet rows=15, profile comb_well_count=14
[checks] FAIL lane_index_is_full_set: lane_index set = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], expected 1..14 each once
[checks] FAIL well_number_is_permutation: well_number set = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], expected a permutation of 1..14
[checks] FAIL lane_index_matches_profile: sheet lanes [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15] vs profile lanes [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
[checks] FAIL label_and_prep_source_explicit: samples need real label/prep; non-samples need not_applicable (never blank); lane 1 (ladder): sample_label must be not_applicable, not 'ladder'; lane 2: sample needs a real prep_source; lane 11 (empty): sample_label must be not_applicable, not 'empty'; lane 12 (empty): sample_label must be not_applicable, not 'empty'; lane 13 (empty): sample_label must be not_applicable, not 'empty'; lane 14 (empty): sample_label must be not_applicable, not 'empty'; lane 15 (empty): sample_label must be not_applicable, not 'empty'
[checks] ERROR: 5 hard failure(s); sample sheet is not usable. See /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.22_20.11.54_rotated_LM-0008_s0015_SDSPAGE_input-double-check-of-ORC-preps_Fl-Green_gel_analysis/sample_sheet_checks.json.
Maybe its better to do more than 10 but less than 15 as ok? Not sure how to deal with the issue of the lanes I cant see. So maybe its just incorrect.
Deleted the lane 15 entry as simple fix.
[checks] FAIL label_and_prep_source_explicit: samples need real label/prep; non-samples need not_applicable (never blank); lane 1 (ladder): sample_label must be not_applicable, not 'ladder'; lane 2: sample needs a real prep_source; lane 11 (empty): sample_label must be not_applicable, not 'empty'; lane 12 (empty): sample_label must be not_applicable, not 'empty'; lane 13 (empty): sample_label must be not_applicable, not 'empty'; lane 14 (empty): sample_label must be not_applicable, not 'empty'
[checks] ERROR: 1 hard failure(s); sample sheet is not usable. See /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.22_20.11.54_rotated_LM-0008_s0015_SDSPAGE_input-double-check-of-ORC-preps_Fl-Green_gel_analysis/sample_sheet_checks.json
[checks] FAIL label_and_prep_source_explicit: samples need real label/prep; non-samples need not_applicable (never blank); lane 1 (ladder): sample_label must be not_applicable, not 'ladder'; lane 2: sample needs a real prep_source
[checks] ERROR: 1 hard failure(s); sample sheet is not usable. See /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.22_20.11.54_rotated_LM-0008_s0015_SDSPAGE_input-double-check-of-ORC-preps_Fl-Green_gel_analysis/sample_sheet_checks.json.
[checks] FAIL label_and_prep_source_explicit: samples need real label/prep; non-samples need not_applicable (never blank); lane 2: sample needs a real prep_source
[checks] ERROR: 1 hard failure(s); sample sheet is not usable. See /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.22_20.11.54_rotated_LM-0008_s0015_SDSPAGE_input-double-check-of-ORC-preps_Fl-Green_gel_analysis/sample_sheet_checks.json.
usb[ ]luis@Luis:~/personal_repos/lab_utils/images$ uv run validate_sample_sheet.py "$GEL_PATH"
[input] gel analysis directory: /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.22_20.11.54_rotated_LM-0008_s0015_SDSPAGE_input-double-check-of-ORC-preps_Fl-Green_gel_analysis
[input] using sample sheet: /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.22_20.11.54_rotated_LM-0008_s0015_SDSPAGE_input-double-check-of-ORC-preps_Fl-Green_gel_analysis/sample_sheet.csv
[input] using profile CSV: /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.22_20.11.54_rotated_LM-0008_s0015_SDSPAGE_input-double-check-of-ORC-preps_Fl-Green_gel_analysis/manual_lane_profiles.csv
[profile] profile has 14 lanes, comb_well_count=14
[read] sample sheet has 14 rows and columns: lane_index, well_number, lane_content, analysis_role, condition_type, sample_label, prep_source
[identity] lane-to-well relationship: identity (well_number == lane_index)
[checks] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.22_20.11.54_rotated_LM-0008_s0015_SDSPAGE_input-double-check-of-ORC-preps_Fl-Green_gel_analysis/sample_sheet_checks.json; overall warning
[checks] WARNING positive_control_present: no lane marked condition_type=positive_control
[checks] WARNING negative_control_present: no lane marked condition_type=negative_control
[done] sample sheet valid; relationship is identity (well_number == lane_index); 2 warning(s).
The checks were rough but I was able to fix them at least. Could be clearer.
Ran the analyze gel
[reference] WARNING lane 8 (well 8) scores clearly higher as reference than your pick, lane 3 (well 3)
[measure] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.22_20.11.54_rotated_LM-0008_s0015_SDSPAGE_input-double-check-of-ORC-preps_Fl-Green_gel_analysis/band_measurement_overlay.png
[done] measured 12 bands x 8 lanes; 0 saturated, 38 fragile; 1 reference warning(s). Inspect band_measurement_overlay.png, band_measurements.csv, and reference_quality.json.
usb[ ]luis@Luis:~/personal_repos/lab_utils/images$ uv run extract_lane_values.py /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.22_20.11.54_rotated_LM-0008_s0015_SDSPAGE_input-double-check-of-ORC-preps_Fl-Green_gel_analysis --region 62.2 65.0
[output] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.22_20.11.54_rotated_LM-0008_s0015_SDSPAGE_input-double-check-of-ORC-preps_Fl-Green_gel_analysis/extract_region_62.2-65mm_none.csv
[output] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.22_20.11.54_rotated_LM-0008_s0015_SDSPAGE_input-double-check-of-ORC-preps_Fl-Green_gel_analysis/extract_region_62.2-65mm_none_checks.json
[summary] gel 2026.08.22_20.11.54_rotated_LM-0008_s0015_SDSPAGE_input-double-check-of-ORC-preps_Fl-Green_gel_analysis
[summary] selection: {"mode": "region", "start_millimetres": 62.2, "end_millimetres": 65.0, "net_baseline": "none"}
[value] lane  1 well    role not_in_measurements                 value 2953309.0
[value] lane  2 well  2 role measured             wt             value 18575366.0
[value] lane  3 well  3 role reference            4r             value 12069448.0
[value] lane  4 well  4 role measured             1ek            value 11914628.0
[value] lane  5 well  5 role measured             3pl            value 13486116.0
[value] lane  6 well  6 role measured             4ps            value 8402105.0
[value] lane  7 well  7 role measured             5ek            value 17371812.0
[value] lane  8 well  8 role measured             6ek            value 7892439.0
[value] lane  9 well  9 role measured             unknown        value 17030184.0
[value] lane 10 well    role not_in_measurements                 value 555891.0
[value] lane 11 well    role not_in_measurements                 value 415570.0
[value] lane 12 well    role not_in_measurements                 value 429134.0
[value] lane 13 well    role not_in_measurements                 value 447964.0
[value] lane 14 well    role not_in_measurements                 value 324987.0
[summary] n=14 min=324987.0 median=8147272.0 max=18575366.0 CV=0.8883
gel picker was probably the coolest idea. Life realtime feedback with rich visual data. Bret victor would approve (even though its kinda trash and not impreesive at the moment, I think its very cool.)
After this, manually processed the excel file to calculate the updated concentrations. Next steps for that is to perform more experiments. Now I can move on to the gel shift files.
I really do just overengineer.
Good. After fixing the sample sheet, running the pipeline was relatively easy.
Did two regions per gel for gel shift assays. Pretty nice. Also proper metric is top / top + bottom region/band. Avoids dividing or denominator as zero or small number.
Made the calculate_gel_shift_ratio.py
