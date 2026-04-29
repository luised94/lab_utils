---
Title: Install apbs for PyMOL3
Date: 2025-11-26
---

# Introduction
APBS is a program/pymol that can calculate electrostatic surface potential of a molecule chain in pymol.
I use it to generate pictures for displaying different characteristics of proteins that I am interested in. 
I encountered some problems recently. Running the APBS via the gui resulted in error 3221225781 (apparently a DLL not found error). I must have installed apbs a long time ago on my old computer. After checking the executable (in the AppDAta folder), the binary did not execute correctly but the .in and .pqr file were created succesfully. Likely need to install the apbs for this device.
This document writes down how I installed apbs on the new device.

# Installation
## Errors encountered:
1. Had to disable windows defender real-time protection and unblock the file (via Properties).
2. Tried running the bin apbs.exe --version and no output was found.
3. Tried apbsplugin.py but that doesnt work (not super sure what it does).
4. Clicked the exe and that reported that python39.dll was missing.
5. Tried to install python via powershell and move dll to the working directory and that did not work.
6. Uninstall PyMOL3 to ensure clean install.
7. Install Visual C++ Redistributable.
8. Install pdb2pqr via pip (the python 3.9 we used).
9. Download Pymol3 and the license again.
10. Restart computer. Still no progress.
11. Had to install miniconda and Redistributables from 2013 then point to the location of those exe.
    > Need to include schrodinger to find the apbs package
    > conda create -n apbs-env python=3.9 apbs pdb2pqr -c schrodinger -c conda-forge

## Succesful steps
1. Download [apbs](https://github.com/Electrostatics/apbs/releases/download/v3.4.1/APBS-3.4.1.Windows.zip) from Github.
    > Had to disable virus protection...
2. Extract to a path of choice.


# Reference
1. https://www.perplexity.ai/search/i-am-getting-an-error-using-ap-rRp1KCO0RXuULk28FGbNXw
2. https://github.com/Electrostatics/apbs
3. https://www.perplexity.ai/search/install-apbs-on-pymol-version-ciXT3zOsTq.mSsaW6K.0Kw
