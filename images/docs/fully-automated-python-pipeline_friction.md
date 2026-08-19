# friction
date: 2026-08-18

## Stage 0
Document issues found while trying to quantify LM-0008_s0012 to measure MCM concentration and finally do the last two experiments.
Have three tmux panes. One at images directory, another with protocol and friction files open in nvim, and last one in the experiment 0008 directory.

Right now, know I can open the PROTOCOL.md document but only from memory. Set GEL to the target image. Use the "printf '%s\n' "$(pwd)"/*.tif". Copy paste from the output then set the GEL variable.

GEL="/mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green.tif"

IJ_EXE="C:\Users\Luised94\Desktop\ij154-win-java8\ImageJ\ImageJ.exe"

Open the file from the command line.
cmd.exe /c "$IJ_EXE" "$(wslpath -w "$GEL")"

The protocol has a lot of info accounting for the orientation. We shall see if that matters as much. The current gel is inverted along the loading axis but migration axis is correct.

```bash
usb[ ]luis@Luis:~/personal_repos/lab_utils/images$ uv run inventory_tiff_tags.py "$GEL" > tags.txt
[path] normalized to /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green.tif
[readability] opened and read 1 byte; reported size is 11535548 bytes
[inventory] complete; 1 page(s) inventoried, nothing interpreted
```

Ran uv run prepare_gel_image.py "$GEL" before preparing the sidecar.

```bash
usb[ ]luis@Luis:~/personal_repos/lab_utils/images$ uv run prepare_gel_image.py "$GEL"
[arguments] input image path as given: '/mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green.tif'
[input_image_path] normalized to /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green.tif
[stage 0] path normalization complete
[inf sidecar] WARNING: no sidecar at /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green.inf. Pixel size, orientation, ScaleType, Invert and PMT voltage are then unavailable and must come from the container's own tags or be recorded as unknown.
[readability] input image: opened and read 1 byte; reported size is 11535548 bytes
[stage 0] input checks complete
[output directory] created and confirmed writable: /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis
input_image_absolute_path       /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green.tif
input_image_physical_path       /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green.tif
inf_sidecar_absolute_path       /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green.inf
inf_sidecar_is_present  no
output_directory_path   /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis
[stage 0] complete; no image data was read and 7 run log lines accumulated for the stage 2 provenance JSON
```
Should I split computations into the ones that require the sidecar and the ones that dont? Clear that sidecar is not present. What is the deal with the inf sidecar? That should not be required.
Unclear what I should copy and how to fill it at this stage (point to PROTOCOL.md?). Go back to protocol file.
Copy the template into notepad (windows), paste, save as. Minor annoyance. Had to create new tmux window.
Then ran:
```bash
usb[ ]luis@Luis:/mnt/c/Users/Luised94$ cd Desktop/lab/experiments/*0008*
usb[ ]luis@Luis:/mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps$ nvim *.txt
```
The template has a lot of comments. Wonder if they should be reduced. Do I read protocol or the template at this point?
Continued with protocol and then looked at the template text.
Record identity good until gel_migration_axis. Protocol was vague but template was good.
Drew the line following the protocol instructions. Angle was around 179. Left the raw values in the file as comments just in case.
Angle seems high but dont think the gel has to be rotated. Wells maybe are smiling a bit.
Ran macro. pretty straightforward.
Instruction at line 157-158 is ambiguous. File references _x and _y. So unclear. I just put in the order that they were printed.
Drew the rectangle. Output in centimeters. Had to go to Analyze > Set Scale and set to zero then remeasure. 
Had to go to https://www.perplexity.ai/search/1ea59212-df86-4648-acc9-207b4951d02d to get the info again.
Use bx to see what the values are and then fill from there in order.
Rotation provenance seemed like a lot of text and there was still a ton of work to do. Should probably move to top to be honest since I will do rotation, write down, then fill out preprocess sidecar. So everything should be clearer up front. No mention of rotation provenance in PROTOCOL.md
Closed the preprocess sidecar file as I thought there was nothing else to do.

## Stage 1
Ran immediately. No warnings, overviews or anything so I assumed that I should just run.
```bash
usb[ ]luis@Luis:~/personal_repos/lab_utils/images$ uv run validate_gel_image.py "$GEL"
[arguments] input tiff path as given: '/mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green.tif'
[path] normalized to /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green.tif
[readability] opened and read 1 byte; reported size is 11535548 bytes
[output directory] writing into /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis
[tiff] opened; byte order <, bigtiff False, 1 page(s)
[tiff] reading page 0
[tiff] inventoried 21 standard tags and 0 private tags at or above 65000
[statistics] min 958, max 65535, mean 1296.27, median 1068, sd 1709.68, distinct 25560, stride 1, isolated extremes 371
[histogram] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/input_file_histogram.png
[crop preview] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/input_file_crop_preview.png at stride 3
[report] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/input_file_validation_report.json
overall_status  fail
failed_hard_stop_count  1
warning_count   1
validation_report_path  /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/input_file_validation_report.json
histogram_image_path    /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/input_file_histogram.png
crop_preview_image_path /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/input_file_crop_preview.png
[finding] WARNING whole_image_saturation: maximum observed value 65535 held by 462 pixels; 462 pixels at the container maximum 65535; spike ratio at the maximum is 416.22 against the 100 present values below it. This is evidence of clipping. Stage 1 reports it and does not adjudicate: per-band at-ceiling counts and the plateau statistic decide whether any particular measurement is affected, and the plateau statistic is the primary detector because the detector can saturate below the numeric ceiling.
[finding] FAIL preprocess_sidecar_filename_matches_input: sidecar was measured against '/mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green.tif' and the input is '2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green.tif'
[stage 1] complete; overall status fail, 1 failed hard stop(s), 1 warning(s). The report was written either way, so a failure is diagnosable from the file rather than only from this scrollback.
```
No mention that the path is relative then? Good constrain but not mentioned. Went into the preprocess sidecar and modified.
```bash
usb[ ]luis@Luis:~/personal_repos/lab_utils/images$ uv run validate_gel_image.py "$GEL"
[arguments] input tiff path as given: '/mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green.tif'
[path] normalized to /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green.tif
[readability] opened and read 1 byte; reported size is 11535548 bytes
[output directory] writing into /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis
[tiff] opened; byte order <, bigtiff False, 1 page(s)
[tiff] reading page 0
[tiff] inventoried 21 standard tags and 0 private tags at or above 65000
[statistics] min 958, max 65535, mean 1296.27, median 1068, sd 1709.68, distinct 25560, stride 1, isolated extremes 371
[histogram] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/input_file_histogram.png
[crop preview] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/input_file_crop_preview.png at stride 3
[report] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/input_file_validation_report.json
overall_status  pass
failed_hard_stop_count  0
warning_count   1
validation_report_path  /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/input_file_validation_report.json
histogram_image_path    /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/input_file_histogram.png
crop_preview_image_path /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/input_file_crop_preview.png
[finding] WARNING whole_image_saturation: maximum observed value 65535 held by 462 pixels; 462 pixels at the container maximum 65535; spike ratio at the maximum is 416.22 against the 100 present values below it. This is evidence of clipping. Stage 1 reports it and does not adjudicate: per-band at-ceiling counts and the plateau statistic decide whether any particular measurement is affected, and the plateau statistic is the primary detector because the detector can saturate below the numeric ceiling.
[stage 1] complete; overall status pass, 0 failed hard stop(s), 1 warning(s). The report was written either way, so a failure is diagnosable from the file rather than only from this scrollback.
```
Inspected the files. Looks good. JSON is hard to read per se (json issue but just saying that maybe I should write down what to check or something)
Crop file looks good.
The warnings next. Breaks flow I think and perhaps they should become active checks and go into gel imaging protocol, not analysis.
Went straight to stage 2 after reading.

## Stage 2

```bash
usb[ ]luis@Luis:~/personal_repos/lab_utils/images$ uv run measure_gel.py "$GEL"
[provenance] input SHA-256 4e18597d7900673e98be803076e523248a15c36f76f5058dbdc27e50db64f421 over 11535548 bytes
[overlay] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/lane_grid_overlay.png
[profiles] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/lane_profiles.png
[profiles] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/lane_profiles.csv
[profiles] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/lane_migration_profiles.csv
[bands] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/band_definitions.csv
[csv] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/band_measurements.csv
[map] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/consensus_center_map.png
[overlay] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/overlay_consensus_windows_by_occupancy.png
[overlay] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/overlay_band_windows_labeled.png
[profiles] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/lane_migration_profiles.png
[csv] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/band_centers.csv
[saturation] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/saturation_derivation.png
[baseline] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/baseline_comparison.png
[figure] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/figure_strip_trace.png
[figure] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/figure_metric_compare.png
[figure] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/figure_net_vs_apex.png
[figure] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/overlay_reported_values.png
[report] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/2026.08.18_11.57.04_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/band_detection_report.json; overall warning
```
No errors. Some warnings. Lets see the files. Maybe files should be ordered in the inspection order and also added a label of what they are (inspect, quality control, quantification/results)

Overall ok. Seems like the two edge lanes were not flagged or identified correctly. Seems like two of the lanes overlapped into one.
The ones that were flagged seem pretty good.
After this, there is no indication of what to do next. The csvs are there but they are kind of incomprehensible just opening. Lots of data.
10 lanes with two ladders. So 12 lanes in a 15 well gel. Only eight properly identified. The line down the lane does show that the lanes should be rotated I think. Try either fiji manual method or rotation?
Cant continue to next step since the identification did not work.
overlay reported values is not interpretable. the letters are too small.

Seems like rotation should be the next step but there is no instructions for that at this stage. No convention for denoting rotated tif file.
Copied the preprocess file and renamed. Reopened the tif file.
otated the file -1.00, grid line 1,bilinear. Transfered the info to preprocess file.
Save as, replaced the timestamp with manual 20280818_rotated prefix. 
Reopened the file in imagej and the preprocess file in notepad. 
Had to go back to protocol just to double check the process.
Did the set scale first.
measured against image remained unchanged.
Updated everything. Pretty straightforward the second time (even with the rotation provenance section added)
Tmux pain that the cmd.exe remains busy until imagej is closed.
Closed imagej now.
Set GEL variable again.
```bash
GEL="/mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green.tif"
usb[ ]luis@Luis:~/personal_repos/lab_utils/images$ uv run validate_gel_image.py "$GEL"
[arguments] input tiff path as given: '/mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green.tif'
[path] normalized to /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green.tif
[readability] opened and read 1 byte; reported size is 11537970 bytes
[output directory] ERROR: /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis does not exist. Run prepare_gel_image.py on this file first; stage 1 does not create the output directory, so that it cannot run against a directory stage 0 never validated.
usb[ ]luis@Luis:~/personal_repos/lab_utils/images$ uv run prepare_gel_image.py "$GEL"
[arguments] input image path as given: '/mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green.tif'
[input_image_path] normalized to /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green.tif
[stage 0] path normalization complete
[inf sidecar] WARNING: no sidecar at /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green.inf. Pixel size, orientation, ScaleType, Invert and PMT voltage are then unavailable and must come from the container's own tags or be recorded as unknown.
[readability] input image: opened and read 1 byte; reported size is 11537970 bytes
[stage 0] input checks complete
[output directory] created and confirmed writable: /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis
input_image_absolute_path       /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green.tif
input_image_physical_path       /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green.tif
inf_sidecar_absolute_path       /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green.inf
inf_sidecar_is_present  no
output_directory_path   /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis
[stage 0] complete; no image data was read and 7 run log lines accumulated for the stage 2 provenance JSON
usb[ ]luis@Luis:~/personal_repos/lab_utils/images$ uv run validate_gel_image.py "$GEL"
[arguments] input tiff path as given: '/mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green.tif'
[path] normalized to /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green.tif
[readability] opened and read 1 byte; reported size is 11537970 bytes
[output directory] writing into /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis
[tiff] opened; byte order >, bigtiff False, 1 page(s)
[tiff] reading page 0
[tiff] inventoried 17 standard tags and 0 private tags at or above 65000
[statistics] min 0, max 65535, mean 1285.76, median 1067, sd 1708.67, distinct 25503, stride 1, isolated extremes 167
[histogram] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/input_file_histogram.png
[crop preview] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/input_file_crop_preview.png at stride 3
[report] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/input_file_validation_report.json
overall_status  pass
failed_hard_stop_count  0
warning_count   4
validation_report_path  /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/input_file_validation_report.json
histogram_image_path    /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/input_file_histogram.png
crop_preview_image_path /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/input_file_crop_preview.png
[finding] WARNING input_file_provenance: ImageDescription carries an ImageJ resave stamp, so this file has been resampled and its instrument tags may be stripped; its pixels may no longer be raw counts. Point the pipeline at the pristine instrument scan. See FINDINGS.md section 1.
[finding] WARNING whole_image_saturation: maximum observed value 65535 held by 376 pixels; 376 pixels at the container maximum 65535; spike ratio at the maximum is 372.28 against the 100 present values below it. This is evidence of clipping. Stage 1 reports it and does not adjudicate: per-band at-ceiling counts and the plateau statistic decide whether any particular measurement is affected, and the plateau statistic is the primary detector because the detector can saturate below the numeric ceiling.
[finding] WARNING orientation_sources_agree: no orientation statement in tag 274, in ImageDescription, or in a .inf. Row order is then unverified and must be confirmed by eye on the first overlay before any band label is trusted.
[finding] WARNING sidecar_rotation_interpolation: sidecar records a bilinear rotation, so the image was resampled by neighbour blending and its pixel values are no longer raw counts. This stacks on encoding_verified: false. Prefer placing the gel square and recording rotation_interpolation_method=none.
[stage 1] complete; overall status pass, 0 failed hard stop(s), 4 warning(s). The report was written either way, so a failure is diagnosable from the file rather than only from this scrollback.
usb[ ]luis@Luis:~/personal_repos/lab_utils/images$ uv run measure_gel_image.py "$GEL"

error: Failed to spawn: `measure_gel_image.py`
  Caused by: No such file or directory (os error 2)
usb[ ]luis@Luis:~/personal_repos/lab_utils/images$ uv run measure_gel.py "$GEL"
[provenance] input SHA-256 77437c8672055d0f6e5d6a26ee5688f8ee6b4c3bf4b148bbb3c13adf26fc9e33 over 11537970 bytes
[overlay] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/lane_grid_overlay.png
[profiles] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/lane_profiles.png
[profiles] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/lane_profiles.csv
[profiles] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/lane_migration_profiles.csv
[bands] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/band_definitions.csv
[csv] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/band_measurements.csv
[map] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/consensus_center_map.png
[overlay] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/overlay_consensus_windows_by_occupancy.png
[overlay] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/overlay_band_windows_labeled.png
[profiles] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/lane_migration_profiles.png
[csv] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/band_centers.csv
[saturation] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/saturation_derivation.png
[baseline] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/baseline_comparison.png
[figure] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/figure_strip_trace.png
[figure] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/figure_metric_compare.png
[figure] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/figure_net_vs_apex.png
[figure] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/overlay_reported_values.png
[report] wrote /mnt/c/Users/Luised94/Desktop/lab/experiments/20260708_LM-0008_load_repeat2-KGlut300mM-allsupps/20260818_rotated_LM-0008_s0012_SDSPAGE_run-old-MCM-aliquots_Fl-Green_gel_analysis/band_detection_report.json; overall warning
```
Looked at the output. Missed lanes were still not detected. Very surprising. If the smiling could be corrected, then I think it would work. Start loading all lanes with buffer just in case.
My sense is that I should do the manual method wiht imagej and that certain macros for saving the intermediates would be all I need to really automate. Let me see. The images investment would have been bad but still. Some progress and never would have known. Still chance that the python approach can be fixed but probably not.

## Manual imagej
I suspect I can do the manual protocol and just need to create a few macros or extract part of the python file to run on the lane profiles.
Maybe do as single pass and write down transformations. Only go with crop for figures, not for analysis.
Draw rectangle > Select Gel one > move using keyboard and then select another lane (Ctrl 2) > Plot (Ctrl 3). This is where the script can come in.
The process was actually pretty fast, skipping the saving of any intermediates.
So what no? Can the plots be converted to csv or something?
Started from the already rotated image.

