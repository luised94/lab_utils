
# Description: Download hawkins origin timing data for categorical analysis and track plotting
# Usage: $Rscript 003_downloadGawkinsTimingData.R
#TODO: Test that the script downloads the data.
#TODO: Process the files into a single format.
hawkins_timing_url <- "https://ars.els-cdn.com/content/image/1-s2.0-S2211124713005834-mmc2.xlsx"
dir.create
hawkins_timing_path <- paste(Sys.getenv("HOME"), "data", "feature_files", "hawkins-origins-timing.xlsx", sep = "/")

curl::curl_download(hawkins_timing_url, hawkins_timing_path)

eaton_peaks_url <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM424nnn/GSM24494/suppl/GSM424494_wt_G2_orc_chip_combiend.bed.gz"
curl::curl_download(eaton_peaks_url, feature_folder)
