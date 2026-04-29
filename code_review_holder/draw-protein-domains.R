library(ggplot2)
tsv_files <- list.files("C:\\Users\\Luis\\Dropbox (MIT)\\Lab\\Projects\\automate-the-boring-stuff\\reference_data\\Interpro-tsv\\", pattern = "*.tsv$", full.names = TRUE)
data_list <- lapply(tsv_files, read.delim)
df <- do.call(rbind, data_list)
df <- unique(df)
df_names <- names(df)

begin_end <- strsplit(df$Matches, "..", fixed = TRUE)
result_df <- data.frame(
  begin = unlist(lapply(begin_end, function(x) as.numeric(x[[1]]))),
  end = unlist(lapply(begin_end, function(x) as.numeric(x[[2]])))
)

result_df$length <- result_df$end - result_df$begin

columns_to_exclude <- c("Integrated", "GO", "Source","Accession", "Matches", "length")
index_to_exclude <- unlist(lapply(columns_to_exclude, function(x){
  grep(x, names(df))
}))
names(df)[grep("Protein.Accession", names(df))] <- "accession"
df$accession <- toupper(df$accession)
features_df <- drawProteins::get_features(paste(unique(df$accession), collapse = " "))
features_df <- drawProteins::feature_to_dataframe(features_df)
features_df
# output_file <- "C:\\Users\\Luis\\Dropbox (MIT)\\Lab\\Projects\\automate-the-boring-stuff\\reference_data\\Interpro-tsv\\features_example.csv"
# write.csv(features_df, file = output_file)

# df[grep("P54784", df$accession),c(2,4,8,10)]
# write.csv(features_df[grep("DOMAIN", features_df$type),], file = paste(dirname(output_file), "/", "domain_info.csv", sep = ""))

features_df <- readxl::read_xlsx("C:\\Users\\Luis\\Dropbox (MIT)\\Lab\\Projects\\automate-the-boring-stuff\\reference_data\\loading_proteins_features.xlsx")
p <- drawProteins::draw_canvas(features_df)
p <- drawProteins::draw_chains(p, features_df)
p
p <- drawProteins::draw_domains(p ,features_df)
p
p<- p + theme_bw(base_size = 20) + # white background
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(), legend.position = "none") +
  theme(axis.ticks = element_blank(),
        axis.text.y = element_blank()) +
  theme(panel.border = element_blank())
p
ggsave("drawProteinDomains.svg", p)

# names(df)[grep("Protein.Length", names(df))] <- "length"
# df$type <- toupper(df$type)
# df_features$accession <- toupper(df_features$accession)
# df$accession <- toupper(df$accession)
# accession_df <- unique(df_features[,c("accession", "entryName")])
# 
# merged_df <- merge(df, accession_df, by = 'accession', all.x = TRUE)
# 
# 
# index_to_exclude <- unlist(lapply(columns_to_exclude, function(x){
#   grep(x, names(merged_df))
# }))
# 
# df2 <- merged_df[grep("DOMAIN", merged_df$type),setdiff(1:ncol(merged_df), index_to_exclude)]
# df3<- cbind(df2,result_df)
# 
# df3$taxid <- df_features$taxid[1]
# 
# 
# df_features[grep("CHAIN", df_features$type),]
# names(df_features)
# names(df3)
# merge(df3, df_features[grep("CHAIN", df_features$type),], by = 'accession', all.x = TRUE)
# df_features[grep("CHAIN", df_features$type),]
# 
# 
# 
# columns_to_include <- c("entryName", "order")
# index_to_include <- unlist(lapply(columns_to_include, function(x){
#   grep(x, names(df_features))
# }))
# 
# 
# df_features[grep("CHAIN", df_features$type), index_to_include]
# 
# df3 <- merge(df3, df_features[grep("CHAIN", df_features$type), index_to_include], by = 'entryName', all.x = TRUE)
# test_df <- rbind(df3,df_features[grep("CHAIN", df_features$type),])
# 
# sort(test_df)
# test_df <- test_df[, names(df_features)]
# grep("AAA", test_df$description)
# grep("CHAIN", test_df$type)
# filtered_df <- unique(test_df[c(grep("AAA", test_df$description),grep("CHAIN", test_df$type)),])
# p <- drawProteins::draw_canvas(filtered_df)
# p <- drawProteins::draw_chains(p, filtered_df)
# p <- drawProteins::draw_domains(p ,filtered_df)
# p<-p + theme_bw(base_size = 20) + # white background
#   theme(panel.grid.minor=element_blank(),
#         panel.grid.major=element_blank(), legend.position = "none") +
#   theme(axis.ticks = element_blank(),
#         axis.text.y = element_blank()) +
#   theme(panel.border = element_blank())
# p
# 
# 
# df2 <- df_features[grep("DOMAIN|CHAIN", df_features$type),]
# p <- drawProteins::draw_canvas(df2)
# p <- drawProteins::draw_chains(p, df2)
# p <- drawProteins::draw_domains(p ,df2)
# p<-p + theme_bw(base_size = 20) + # white background
#   theme(panel.grid.minor=element_blank(),
#         panel.grid.major=element_blank(), legend.position = "none") +
#   theme(axis.ticks = element_blank(),
#         axis.text.y = element_blank()) +
#   theme(panel.border = element_blank())
# p

