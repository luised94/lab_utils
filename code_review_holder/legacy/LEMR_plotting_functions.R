library(tidyverse)
library(ggplot2)
library(plyr)
#Functions work on dataframes. Take a list of indexes to designate x and y variables, index to color data by, index to filter dataframe by, and index to facet data.
#Current functions only work for geom_point, _density,_box geometries of ggplot2. Filtering only tested for factor types using strings to filter by equality (==)
#to adjust size, plot windown is modified. May add extra input with default parameters to specify max size.
#Color index has to be specified always. Otherwise, color will be character(0) and will not work. Color by factors or boolean works well. 

plot_point_and_save <- function(df, x_index, y_index, facet_index = 0, filter_index = 0, filter_factor = "None", color_index, extra_String = "", file_type = "png"){
  df_columns = colnames(df)
  if(facet_index == 0 & filter_index == 0){
    for (i in 1:length(x_index)){
      for (j in 1:length(y_index)){
        if (df_columns[x_index[i]] != df_columns[y_index[j]]){
          plt <- df %>% ggplot(aes_string(x = df_columns[x_index[i]], y = df_columns[y_index[j]], color = df_columns[color_index])) + geom_point() 
          name_for_file <- paste(df_columns[x_index[i]], "vs", df_columns[y_index[j]], "col", df_columns[color_index], extra_String, ".png", sep = "_")
          ggsave(name_for_file, device = file_type, plot = plt)
        }
      }
    } 
  }else if (facet_index == 0 & filter_index != 0){
    for (i in 1:length(x_index)){
      for (j in 1:length(y_index)){
        if (df_columns[x_index[i]] != df_columns[y_index[j]]){
          plt <- df %>% filter(df[df_columns[filter_index]] == filter_factor) %>% ggplot(aes_string(x = df_columns[x_index[i]], y = df_columns[y_index[j]], color = df_columns[color_index])) + geom_point() 
          name_for_file <- paste(df_columns[x_index[i]], "vs", df_columns[y_index[j]], "filter", df_columns[filter_index], "col", df_columns[color_index], extra_String, ".png", sep = "_")
          ggsave(name_for_file, device = file_type, plot = plt)
        }
      }
    }
  }else if (facet_index != 0 & filter_index == 0){
    for (i in 1:length(x_index)){
      for (j in 1:length(y_index)){
        if (df_columns[x_index[i]] != df_columns[y_index[j]]){
          plt <- df %>% ggplot(aes_string(x = df_columns[x_index[i]], y = df_columns[y_index[j]], color = df_columns[color_index])) + geom_point() + 
           facet_wrap(as.formula(paste(df_columns[facet_index],"~", ".", sep = "")), nrow = 2)
          name_for_file <- paste(df_columns[x_index[i]], "vs", df_columns[y_index[j]], "by", df_columns[facet_index],"col", df_columns[color_index], extra_String, ".png", sep = "_")
          ggsave(name_for_file, device = file_type, plot = plt)
        }
      }
    }
  }else{
    for (i in 1:length(x_index)){
      for (j in 1:length(y_index)){
        if (df_columns[x_index[i]] != df_columns[y_index[j]]){
          plt <- df %>% filter(df[df_columns[filter_index]] == filter_factor) %>% ggplot(aes_string(x = df_columns[x_index[i]], y = df_columns[y_index[j]], color = df_columns[color_index])) + geom_point() + 
                facet_wrap(as.formula(paste(df_columns[facet_index],"~", ".", sep = "")), nrow = 2)
          name_for_file <- paste(df_columns[x_index[i]], "vs", df_columns[y_index[j]], "by", df_columns[facet_index], "filter", df_columns[filter_index], "col", df_columns[color_index], extra_String, ".png", sep = "_")
          ggsave(name_for_file, device = file_type, plot = plt)
        }
      }
    }
  } 
}
plot_density_and_save <- function(df, x_index, facet_index= 0, filter_index = 0, filter_factor = "None", color_index, extra_String="", file_type = "png"){
  df_columns = colnames(df)
  if(facet_index == 0 & filter_index == 0){
    for (i in 1:length(x_index)){
      plt <- df %>% ggplot(aes_string(x = df_columns[x_index[i]], color = df_columns[color_index])) + geom_density(alpha = .4)  
      name_for_file <- paste(df_columns[x_index[i]], "density", "col", df_columns[color_index], extra_String, ".png", sep = "_")
      ggsave(name_for_file, device = file_type, plot = plt)
    } 
  }else if (facet_index == 0 & filter_index != 0){
    for (i in 1:length(x_index)){
      plt <- df %>% filter(df[df_columns[filter_index]] == filter_factor) %>% ggplot(aes_string(x = df_columns[x_index[i]], color = df_columns[color_index])) + geom_density(alpha = .4) 
      name_for_file <- paste(df_columns[x_index[i]], "density", "filter", filter_factor, "col", df_columns[color_index], extra_String,".png", sep = "_")
      ggsave(name_for_file, device = file_type, plot = plt)
    }
  }else if (facet_index != 0 & filter_index == 0){
    for (i in 1:length(x_index)){
      plt <- df %>% ggplot(aes_string(x = df_columns[x_index[i]], color = df_columns[color_index])) + geom_density(alpha = .4) + 
        facet_wrap(as.formula(paste(df_columns[facet_index],"~", ".", sep = "")), nrow = 2)
      name_for_file <- paste(df_columns[x_index[i]], "density", "by", df_columns[facet_index], "col", df_columns[color_index], extra_String, ".png", sep = "_")
      ggsave(name_for_file, device = file_type, plot = plt)
    }
  }else{
    for (i in 1:length(x_index)){
    plt <- df %>% filter(df[df_columns[filter_index]] == filter_factor) %>% ggplot(aes_string(x = df_columns[x_index[i]], color = df_columns[color_index])) + geom_density(alpha = .4) + 
      facet_wrap(as.formula(paste(df_columns[facet_index],"~", ".", sep = "")), nrow = 2)
    name_for_file <- paste(df_columns[x_index[i]], "density", "by", df_columns[facet_index], "filter", filter_factor,"col", df_columns[color_index], extra_String, ".png", sep = "_")
    ggsave(name_for_file, device = file_type, plot = plt)
    }
  }
}  
plot_box_and_save <- function(df, x_index, y_index, facet_index = 0, filter_index = 0, filter_factor = "None", color_index, extra_String="", file_type = "png"){
  df_columns = colnames(df)
  if(facet_index == 0 & filter_index == 0){
    for (i in 1:length(x_index)){
      for (j in 1:length(y_index)){
        if (df_columns[x_index[i]] != df_columns[y_index[j]]){
          plt <- df %>% group_by(df_columns[x_index[i]])%>% ggplot(aes_string(x = df_columns[x_index[i]], y = df_columns[y_index[j]], color = df_columns[color_index])) +
            geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) 
          name_for_file <- paste("Box", df_columns[x_index[i]], "vs", df_columns[y_index[j]], "col", df_columns[color_index], extra_String, ".png", sep = "_")
          ggsave(name_for_file, device = file_type, plot = plt)
        }
      }
    }
  } else if (facet_index == 0 & filter_index != 0){
    for (i in 1:length(x_index)){
      
      for (j in 1:length(y_index)){
        
        if (df_columns[x_index[i]] != df_columns[y_index[j]]){
          
          plt <- df %>% group_by(df_columns[x_index[i]]) %>% filter(df[df_columns[filter_index]] == filter_factor) %>% ggplot(aes_string(x = df_columns[x_index[i]], y = df_columns[y_index[j]], color = df_columns[color_index])) +
            geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) 
          
          name_for_file <- paste("Box", df_columns[x_index[i]], "vs", df_columns[y_index[j]], "filter", df_columns[filter_index], "col", df_columns[color_index], extra_String, ".png", sep = "_")
          
          ggsave(name_for_file, device = file_type, plot = plt)
        }
      }
    }
  } else if (facet_index != 0 & filter_index == 0){
    for (i in 1:length(x_index)){
      
      for (j in 1:length(y_index)){
        
        if (df_columns[x_index[i]] != df_columns[y_index[j]]){
          
          plt <- df %>% group_by(df_columns[x_index[i]]) %>% ggplot(aes_string(x = df_columns[x_index[i]], y = df_columns[y_index[j]], color = df_columns[color_index])) +
            geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) +facet_wrap(as.formula(paste(df_columns[facet_index],"~", ".", sep = "")))
          
          name_for_file <- paste("Box", df_columns[x_index[i]], "vs", df_columns[y_index[j]], "by", df_columns[facet_index], "col", df_columns[color_index], extra_String, ".png", sep = "_")
          
          ggsave(name_for_file, device = file_type, plot = plt)
        }
      }
    }
  } else {
    for (i in 1:length(x_index)){
      
      for (j in 1:length(y_index)){
        
        if (df_columns[x_index[i]] != df_columns[y_index[j]]){
          
          plt <- df %>% group_by(df_columns[x_index[i]]) %>% filter(df[df_columns[filter_index]] == filter_factor) %>% ggplot(aes_string(x = df_columns[x_index[i]], y = df_columns[y_index[j]], color = df_columns[color_index])) +
            geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) +facet_wrap(as.formula(paste(df_columns[facet_index],"~", ".", sep = "")))
          
          name_for_file <- paste("Box", df_columns[x_index[i]], "vs", df_columns[y_index[j]], "by", df_columns[facet_index], "filter", df_columns[filter_index], "col", df_columns[color_index], extra_String, ".png", sep = "_")
          
          ggsave(name_for_file, device = file_type, plot = plt)
        }
      }
    }
  }
}

#Without color
plot_point_and_save_nocolor <- function(df, x_index, y_index, facet_index = 0, filter_index = 0, filter_factor = "None", extra_String = "", file_type = "png"){
  df_columns = colnames(df)
  if(facet_index == 0 & filter_index == 0){
    for (i in 1:length(x_index)){
      for (j in 1:length(y_index)){
        if (df_columns[x_index[i]] != df_columns[y_index[j]]){
          plt <- df %>% ggplot(aes_string(x = df_columns[x_index[i]], y = df_columns[y_index[j]])) + geom_point() 
          name_for_file <- paste(df_columns[x_index[i]], "vs", df_columns[y_index[j]], extra_String,".png", sep = "_")
          ggsave(name_for_file, device = file_type, plot = plt)
        }
      }
    } 
  }else if (facet_index == 0 & filter_index != 0){
    for (i in 1:length(x_index)){
      for (j in 1:length(y_index)){
        if (df_columns[x_index[i]] != df_columns[y_index[j]]){
          plt <- df %>% filter(df[df_columns[filter_index]] == filter_factor) %>% ggplot(aes_string(x = df_columns[x_index[i]], y = df_columns[y_index[j]])) + geom_point() 
          name_for_file <- paste(df_columns[x_index[i]], "vs", df_columns[y_index[j]], "filter", df_columns[filter_index], extra_String,".png", sep = "_")
          ggsave(name_for_file, device = file_type, plot = plt)
        }
      }
    }
  }else if (facet_index != 0 & filter_index == 0){
    for (i in 1:length(x_index)){
      for (j in 1:length(y_index)){
        if (df_columns[x_index[i]] != df_columns[y_index[j]]){
          plt <- df %>% ggplot(aes_string(x = df_columns[x_index[i]], y = df_columns[y_index[j]])) + geom_point() + 
            facet_wrap(as.formula(paste(df_columns[facet_index],"~", ".", sep = "")), nrow = 2)
          name_for_file <- paste(df_columns[x_index[i]], "vs", df_columns[y_index[j]], "by", df_columns[facet_index], extra_String,".png", sep = "_")
          ggsave(name_for_file, device = file_type, plot = plt)
        }
      }
    }
  }else{
    for (i in 1:length(x_index)){
      for (j in 1:length(y_index)){
        if (df_columns[x_index[i]] != df_columns[y_index[j]]){
          plt <- df %>% filter(df[df_columns[filter_index]] == filter_factor) %>% ggplot(aes_string(x = df_columns[x_index[i]], y = df_columns[y_index[j]])) + geom_point() + 
            facet_wrap(as.formula(paste(df_columns[facet_index],"~", ".", sep = "")), nrow = 2)
          name_for_file <- paste(df_columns[x_index[i]], "vs", df_columns[y_index[j]], "by", df_columns[facet_index], "filter", df_columns[filter_index], extra_String,".png", sep = "_")
          ggsave(name_for_file, device = file_type, plot = plt)
        }
      }
    }
  } 
}
plot_density_and_save_nocolor <- function(df, x_index, facet_index= 0, filter_index = 0, filter_factor = "None",  extra_String = "", file_type = "png"){
  df_columns = colnames(df)
  if(facet_index == 0 & filter_index == 0){
    for (i in 1:length(x_index)){
      plt <- df %>% ggplot(aes_string(x = df_columns[x_index[i]])) + geom_density(alpha = .4)  
      name_for_file <- paste(df_columns[x_index[i]], "density",extra_String,".png", sep = "_")
      ggsave(name_for_file, device = file_type, plot = plt)
    } 
  }else if (facet_index == 0 & filter_index != 0){
    for (i in 1:length(x_index)){
      plt <- df %>% filter(df[df_columns[filter_index]] == filter_factor) %>% ggplot(aes_string(x = df_columns[x_index[i]])) + geom_density(alpha = .4) 
      name_for_file <- paste(df_columns[x_index[i]], "density", "filter", filter_factor,extra_String,".png", sep = "_")
      ggsave(name_for_file, device = file_type, plot = plt)
    }
  }else if (facet_index != 0 & filter_index == 0){
    for (i in 1:length(x_index)){
      plt <- df %>% ggplot(aes_string(x = df_columns[x_index[i]])) + geom_density(alpha = .4) + 
        facet_wrap(as.formula(paste(df_columns[facet_index],"~", ".", sep = "")), nrow = 2)
      name_for_file <- paste(df_columns[x_index[i]], "density", "by", df_columns[facet_index], extra_String, ".png", sep = "_")
      ggsave(name_for_file, device = file_type, plot = plt)
    }
  }else{
    for (i in 1:length(x_index)){
      plt <- df %>% filter(df[df_columns[filter_index]] == filter_factor) %>% ggplot(aes_string(x = df_columns[x_index[i]])) + geom_density(alpha = .4) + 
        facet_wrap(as.formula(paste(df_columns[facet_index],"~", ".", sep = "")), nrow = 2)
      name_for_file <- paste(df_columns[x_index[i]], "density", "by", df_columns[facet_index], "filter", filter_factor, extra_String, ".png", sep = "_")
      ggsave(name_for_file, device = file_type, plot = plt)
    }
  }
}  
plot_box_and_save_nocolor <- function(df, x_index, y_index, facet_index = 0, filter_index = 0, filter_factor = "None", extra_String = "", file_type = "png"){
  df_columns = colnames(df)
  if(facet_index == 0 & filter_index == 0){
    for (i in 1:length(x_index)){
      for (j in 1:length(y_index)){
        if (df_columns[x_index[i]] != df_columns[y_index[j]]){
          plt <- df %>% group_by(df_columns[x_index[i]])%>% ggplot(aes_string(x = df_columns[x_index[i]], y = df_columns[y_index[j]])) +
            geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) 
          name_for_file <- paste("Box", df_columns[x_index[i]], "vs", df_columns[y_index[j]], extra_String,".png", sep = "_")
          ggsave(name_for_file, device = file_type, plot = plt)
        }
      }
    }
  } else if (facet_index == 0 & filter_index != 0){
    for (i in 1:length(x_index)){
      for (j in 1:length(y_index)){
        if (df_columns[x_index[i]] != df_columns[y_index[j]]){
          plt <- df %>% group_by(df_columns[x_index[i]]) %>% filter(df[df_columns[filter_index]] == filter_factor) %>% ggplot(aes_string(x = df_columns[x_index[i]], y = df_columns[y_index[j]])) +
            geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) 
          name_for_file <- paste("Box", df_columns[x_index[i]], "vs", df_columns[y_index[j]], "filter", df_columns[filter_index], extra_String, ".png", sep = "_")
          ggsave(name_for_file, device = file_type, plot = plt)
        }
      }
    }
  } else if (facet_index != 0 & filter_index == 0){
    for (i in 1:length(x_index)){
      for (j in 1:length(y_index)){
        if (df_columns[x_index[i]] != df_columns[y_index[j]]){
          plt <- df %>% group_by(df_columns[x_index[i]]) %>% ggplot(aes_string(x = df_columns[x_index[i]], y = df_columns[y_index[j]])) +
            geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) +facet_wrap(as.formula(paste(df_columns[facet_index],"~", ".", sep = "")))
          name_for_file <- paste("Box", df_columns[x_index[i]], "vs", df_columns[y_index[j]], "by", df_columns[facet_index], extra_String, ".png", sep = "_")
          ggsave(name_for_file, device = file_type, plot = plt)
        }
      }
    }
  } else {
    for (i in 1:length(x_index)){
      for (j in 1:length(y_index)){
        if (df_columns[x_index[i]] != df_columns[y_index[j]]){
          plt <- df %>% group_by(df_columns[x_index[i]]) %>% filter(df[df_columns[filter_index]] == filter_factor) %>% ggplot(aes_string(x = df_columns[x_index[i]], y = df_columns[y_index[j]])) +
            geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) + facet_wrap(as.formula(paste(df_columns[facet_index],"~", ".", sep = "")))
          name_for_file <- paste("Box", df_columns[x_index[i]], "vs", df_columns[y_index[j]], "by", df_columns[facet_index], "filter", df_columns[filter_index], extra_String, ".png", sep = "_")
          ggsave(name_for_file, device = file_type, plot = plt)
        }
      }
    }
  }
}

 



# plot_density_and_save <- function(df, x_index, facet_index, color_index, file_type = "png"){
#   df_columns = colnames(df)
#   if (facet_index == 0){
#     for (i in 1:length(x_index)){
#       plt <- Rhee_df %>% ggplot(aes_string(x = df_columns[x_index[i]], color = df_columns[color_index])) + geom_density(alpha = .4)  
#       name_for_file <- paste(df_columns[x_index[i]], "density", ".png", sep = "_")
#       ggsave(name_for_file, device = file_type, plot = plt)
#     } 
#   } else{
#     for (i in 1:length(x_index)){
#       plt <- Rhee_df %>% ggplot(aes_string(x = df_columns[x_index[i]], color = df_columns[color_index])) + geom_density(alpha = .4) + 
#         facet_wrap(as.formula(paste(df_columns[facet_index],"~", ".", sep = "")), nrow = 2)
#       name_for_file <- paste(df_columns[x_index[i]], "density", "by", df_columns[facet_index], ".png", sep = "_")
#       ggsave(name_for_file, device = file_type, plot = plt)
#     } 
#   }
# }
# plot_box_and_save <- function(df, x_index, y_index, facet_index, file_type = "png"){
#   df_columns = colnames(df)
#   if(facet_index == 0){
#     for (i in 1:length(x_index)){
#       
#       for (j in 1:length(y_index)){
#         
#         if (df_columns[x_index[i]] != df_columns[y_index[j]]){
#           
#           plt <- df %>% group_by(df_columns[x_index[i]])%>% ggplot(aes_string(x = df_columns[x_index[i]], y = df_columns[y_index[j]])) +
#             geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) 
#           
#           name_for_file <- paste("Box", df_columns[x_index[i]], "vs", df_columns[y_index[j]], ".png", sep = "_")
#           
#           ggsave(name_for_file, device = file_type, plot = plt)
#         }
#       }
#     }
#   } else{
#     for (i in 1:length(x_index)){
#       
#       for (j in 1:length(y_index)){
#         
#         if (df_columns[x_index[i]] != df_columns[y_index[j]]){
#           
#           plt <- df %>% group_by(df_columns[x_index[i]]) %>% ggplot(aes_string(x = df_columns[x_index[i]], y = df_columns[y_index[j]])) +
#             geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) +facet_wrap(as.formula(paste(df_columns[facet_index],"~", ".", sep = "")))
#           
#           name_for_file <- paste("Box", df_columns[x_index[i]], "vs", df_columns[y_index[j]], "by", df_columns[facet_index], ".png", sep = "_")
#           
#           ggsave(name_for_file, device = file_type, plot = plt)
#         }
#       }
#     }
#   }
# }