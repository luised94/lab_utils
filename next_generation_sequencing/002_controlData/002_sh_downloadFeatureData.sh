#!/bin/bash

download_dir=$HOME/data/features_files

git clone --depth=1 https://github.com/CEGRcode/2021-Rossi_Nature.git "$download_dir"
#wget -O "$download_dir"/02_References_and_Features_Files/hawkins-origin-timing.xlsx "https://ars.eln-cdn.com/content/image/1-s2.0-S2211124713005834-mmc2.xlsx"
#"https://www.cell.com/cms/10.1016/j.cellrep.2013.10.014/attachment/a16ab8a7-4e28-4ca6-aafe-e8857b760c80/mmc2.xlsx 
#"https://ars.eln-cdn.com/content/image/1-s2.0-S2211124713005834-mmc2.xlsx" 
#mv "$download_dir/CEGRcode/2021-Rossi_Nature/02_References_and_Features_Files" "$download_dir"
#curl --output $HOME/data/playground/test.xlsx "https://ars.eln-cdn.com/content/image/1-s2.0-S2211124713005834-mmc2.xlsx"


