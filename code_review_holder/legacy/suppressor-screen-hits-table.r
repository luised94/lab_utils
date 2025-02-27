library(tidyverse)
library(gt)
hits_table <- data.frame(Gene=c("ORC1","ORC3","ORC4","ORC5","ORC6","TOA2"),
                         Mutation = c("E495K","P481L","P225S","E104K","E304K","G16E"))

simple_table <- hits_table %>% gt() %>% 
  tab_header(title = "ORC4-R267A Suppressor Screen Hits") %>%
  cols_label(Gene = md("**Gene**"), Mutation = md("**Amino-Acid Mutation**")) %>% cols_align(align ="center") %>% 
  cols_width(Gene ~ px(120), Mutation ~ px(120)) 

simple_table %>% gtsave(filename = "suppressor-screen-hits-simple-table.png")
