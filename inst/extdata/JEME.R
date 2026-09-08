## code to prepare `JEME` dataset goes here

# Now that download_JEME is an exported function in the hotnetR package,
# you can use it directly once the package is loaded.

# Example usage to prepare the data:
# library(hotnetR)
# 
# jeme_lasso <- download_JEME("lasso")
# JEME_lasso_combined <- jeme_lasso$combined
# 
# jeme_elasticnet <- download_JEME("elasticnet")
# JEME_elasticnet_combined <- jeme_elasticnet$combined
# 
# usethis::use_data(JEME_lasso_combined, overwrite = TRUE)
# usethis::use_data(JEME_elasticnet_combined, overwrite = TRUE)
