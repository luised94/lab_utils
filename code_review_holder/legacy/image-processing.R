library(tidyverse)
library(magick)
library(EBImage)

img <- image_read(list.files(file_folder, pattern = ".jpg", full.names = TRUE)[1])
image_browse(img)
text <- tesseract::ocr(img)
cat(text)
image_browse(img)
image_info(img)

img_processed <- img %>%
  image_crop(geometry = "3700x2900+500+400") %>%
  image_transparent(color = "white", fuzz=20)%>%
  image_quantize(colorspace = 'gray') %>%
  image_browse()
  # image_background("white", flatten = TRUE) %>%
  # image_trim() %>%
  # image_noise() %>%
  # image_enhance() %>%
  # image_normalize() %>%
  # image_contrast(sharpen = 1) %>%
  # image_deskew(threshold = 40) %>%
  # image_ocr()

#Using EBImage
img <- readImage(list.files(file_folder, pattern = ".jpg", full.names = TRUE)[1])
img %>% channel("gray") %>% watershed() %>% display()
