# Load the RDS file
mm_small_intestine_data <- readRDS("data/mm_small_intestine_data.RDS")

# Save the data in the correct format for R package
usethis::use_data(mm_small_intestine_data, overwrite = TRUE) 