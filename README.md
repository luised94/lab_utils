# lab_utils
Code used for laboratory analysis

# STICKY NOTES
This note is just a section of thoughts and comments that documents though process and some notes. 
I have separated it into dates. 

## 2024 04 21
- Most of the bash code I wrote for my CHIP-seq analysis previously served only to run R code. 
- Need to see how I read feature files.
- Can definitely the quality control with bash (deeptools or other)
- Think I finally understand how I was coding. I was essentially debugging and testing but leaving it in that state.
- Furthermore, I think I was using weird unintuitive approaches and building on top of it if it worked.
- Need to output environment info after I finish analysis with git info. 
- Installed all packages in 1-package-installation.r. I am reworking this file so I will just have to remember to document installed and used pakcages.
- Was able to replace most the initial part of the pipeline with bash instead of R. 
- Copied rsync files from previous code collecition  because it had useful parts to extract.

# TODO

1. MED- See Random R questions thread to create a list of packages that have been used and document use in scripts.
2. HIGH- Create genome track plotting code.
3. MED- Create rsync tool for plots of experiments
4. HIGH- Download Eaton Data for controls.
5. LOW- Output environment info and git info to NGS data directory. HIGH once done. 
