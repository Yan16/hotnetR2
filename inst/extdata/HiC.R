## code to prepare `HiC` dataset goes here

# Now that download_HiC is an exported function in the hotnetR package,
# you can use it directly once the package is loaded.

# Example usage to prepare the data:
# library(hotnetR)
# 
# hic_data <- download_HiC()
# HiC_PO <- hic_data$`P-O`
# HiC_PP <- hic_data$`P-P`
# 
# usethis::use_data(HiC_PO, overwrite = TRUE)
# usethis::use_data(HiC_PP, overwrite = TRUE)
