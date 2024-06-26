hawkins_timing_url <- "https://ars.els-cdn.com/content/image/1-s2.0-S2211124713005834-mmc2.xlsx"
hawkins_timing <- paste(Sys.getenv("HOME"), "data", "feature_files", "hawkins-origins-timing.xlsx", sep = "/")

curl::curl_download(hawkins_timing_url, hawkins_timing)
#readxl::read_excel(hawkins_timing)
#df_originTiming <- as.data.frame(readxl::read_excel(hawkins_timing)) %>% dplyr::select(1:7) %>% filter(!is.na(Chromosome))




df_names <- c("ChExMix","MEME","XUTs","CUTs","ORF","Nucleosome", "timing", "Rhee" , "SGD_features", "SUTs", "G2_orc")

for(i in 1:length(df_names)){
  print(feature_files[grepl(df_names[i], feature_files)])
  filepath <- file.path(feature_folder, feature_files[grepl(df_names[i], feature_files)])
  if (grepl("xls|xlsx", tools::file_ext(filepath))){
    assign(df_names[i], readxl::read_excel(filepath))
    if(df_names[i] == "timing"){
      assign(df_names[i], get(df_names[i]) %>% as.data.frame() %>% dplyr::select(1:7) %>% filter(!is.na(Chromosome)))
    }
  } else if(tools::file_ext(filepath) == "bed"){
    assign(df_names[i], readBed(filepath))
    print(seqlevelsStyle(eval(parse(text = df_names[i]))))
  } else if(tools::file_ext(filepath) == "gff3"){
    assign(df_names[i], toGRanges(makeTxDbFromGFF(filepath), format = 'gene'))
    print(seqlevelsStyle(eval(parse(text = df_names[i]))))
  } else if(tools::file_ext(filepath) == "tab"){
    assign(df_names[i], read.delim(filepath, header = FALSE))
    if(df_names[i] == "SGD_features"){
      assign(df_names[i], get(df_names[i]) %>% as.data.frame() %>% dplyr::select(c(9,10,11,12,2,4,5)))
    }
  }
}














