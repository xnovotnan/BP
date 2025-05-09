library(tidyverse)
library(magrittr)
library(png)
options(scipen = 999)

# The get_fastqc_data function reads and returns the content of the fastqc_data.txt file 
# from the selected FASTQC folder, stopping with an error if the file is missing.
get_fastqc_data <- function(zip_file) {
  zip_contents <- unzip(zip_file, list = TRUE)
  fastqc_data_path <- zip_contents$Name[grepl("fastqc_data.txt", zip_contents$Name)]
  if (length(fastqc_data_path) == 0) {
    stop("File fastqc_data.txt not found in this folder.")
  }
  readLines(unz(zip_file, fastqc_data_path))
}

# The process_fastqc function extract relevant data from the FASTQC report 
# and creates a dictionary with these values. 
process_fastqc <- function(fastqc_folder, read) {
  incProgress(0.1, detail = "Reading file...")
  zip_file <- list.files(path = fastqc_folder, pattern = paste0("*", read, "_fastqc.zip"), full.names = TRUE)
  if (length(zip_file) == 0) {
    stop(paste0("Zip file containing read ", read," not found in selected folder."))
  }
  
  lines <- get_fastqc_data(zip_file)
  start_index <- grep("Filename", lines)
  end_index <- grep("Sequence length", lines)
  lines <- lines[start_index:end_index]
  results <- list()
  incProgress(0.4, detail = "Processing data...")
  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    key <- trimws(parts[1])
    value <- trimws(parts[2])
    results[[key]] <- value
  }
  results
}

# The find_fastqc_png function locates the .png file with name from the parameters   
# in the selected FASTQC read subfolder, stopping if the folder or image is missing.
find_fastqc_png <- function(fastqc_folder, read, png_name){
  zip_file <- list.files(path = fastqc_folder, pattern = paste0("*", read, "_fastqc.zip"), full.names = TRUE)
  if (length(zip_file) == 0) {
    stop("Zip file not found.")
  }
  
  temp_contents <- list.dirs(path = "Temp", full.names = FALSE, recursive = FALSE)
  unzipped <- FALSE
  for (folder in temp_contents) {
    if (grepl(folder, zip_file)){
      unzipped <- TRUE
      break
    }
  }
  # if file is not unzipped yet
  if (!unzipped){
    unzip(zip_file, exdir = "Temp")
  }
  directory <- getwd()
  folder_name <- file.path(getwd(), "Temp", substr(basename(zip_file), 1, nchar(basename(zip_file)) - 4), sep="/")
  print(getwd())
  file_path <- file.path(folder_name, paste("Images", png_name, sep="/"))
  
  if (!file.exists(file_path)) {
    stop("Image not found.")
  }
  file_path
}
