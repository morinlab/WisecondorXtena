#!/usr/bin/env Rscript

#R script created by Kristena Daley from Dr. Ryan Morins Lab on June 11, 2021 to determine the mean value that the bins are offset by, and this value will be input into wisecondorX with the offset_neut_peak argument to subsequently shift the segments to their proper place.
#The first model fitting here is used to determine the value to shift by, and the second model fitting with a specific cluster of 3 (G=3) is to calculate purity estimates to use as beta. 

library(mclust)
library(tidyverse)
library(glue)

args = commandArgs(trailingOnly = TRUE)

# Test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Input bins.bed file must be supplied (input file).n", call.=FALSE)
}  else if (length(args)==1) {
  # Default output file
  args[2] = "new_offset_mean_and_purity.txt"
}
in_bed <- read.delim(args[1])

# Get rid of all rows that have blacklisted values (called NaN)
in_bed_nonull <- na.omit(in_bed, cols_only(ratio))


##
## CALCULATE THE OFFSET NEUTRAL PEAK VALUE/MEAN USING G=1
##

# Fit mixture model to the logratio portion of the cleaned up input bed file
model_fit = mclust::Mclust(in_bed_nonull$ratio, G=1, warn = TRUE)

# This makes a beautiful density plot showing the mean of each peak and putting a line through them
# Looks into the parameters section of the model, and find the mean of the 1 cluster that is automatically calculated, 
model_fit$parameters$mean
# modeldf takes the mean from above and puts the  
modeldf=data.frame(mean=model_fit$parameters$mean)
offset_mean=modeldf$mean

# Save the shifted mean peak output as a value to export to a text file
print_offset_mean <- glue("offset_mean = ", offset_mean)


##
## CALCULATE THE PURITY USING G=3
##

# Fit mixture model to the logratio portion of the cleaned up input bed file, specifying 3 clusters to calculate purity (tumour fraction)
model_fit = mclust::Mclust(in_bed_nonull$ratio, G=3, warn = TRUE)

# Take the Mean and Pro of the mixture model, since the highest pro correlates to the highest peak
peak_info <- data.frame(model_fit$parameters$mean, model_fit$parameters$pro)
# Renames the columns
names(peak_info) <- c("Mean", "Intensity")
# Orders the means from lowest to highest
peak_info <- peak_info[order(peak_info$Mean),]
# Calculates the average log ratio of every data point in the bins file
avg_log2 = mean(model_fit$data)

# %>% means pipe, mutate is used to add new variables and preserve existing ones
# This takes the peak_info table that has the means of each classification and matched to their "intensity" or "pro" number form the mixture model, and calculates the "mean_offset" which is the mean of the classification minus the average log2ratio of all bins. 
peak_info <- peak_info %>%
  mutate(mean_offset = abs(Mean - avg_log2))

# Takes the offset mean from the neutral/center peak table
center_peak_offset = peak_info[2,3]

# Adds how offset the neutral peak is from log 0 to each log ratio from wisecondorX and adds that as a new column in the bins.bed dataframe
in_bed_nonull = in_bed_nonull %>%
  mutate(mean_offset = (ratio + center_peak_offset))

# Adds another column called "peak" with the matching classifications from the mixture model
in_bed_nonull$peak = model_fit$classification

# Creates a new file with the original columns from the input_bed_nonull, 
# filters everything out that is in the neutral peak/classification, (keep values classified as 1 or 3)
# case_when is a kind of ifelse statement, saying if the peak is classified as 1 its a loss (this basically calculates tumour content), and 3 is a gain
# divides by the ideal log2(read depth ratio) if purity was at 100%
# case_when(condition ~ output_value) (case_when = name of the function, condition = evaluates as true, output_value = the value to output if condition is true)
# TRUE ~ 1000 is a catch-all
bed_purity = in_bed_nonull %>%
  filter(peak != 2) %>%
  mutate(purity_estimate = case_when(
    peak == 1 ~ mean_offset/-1,
    peak == 3 ~ mean_offset/0.58,
    TRUE ~ 1000
  ))

# the finished purity estimation takes the average of all purity estimates (from all associated/relevant bins)
purity = mean(bed_purity$purity_estimate)

rerun <- 
  if (purity >= 1 & !is.na(purity)){
  print ("TRUE")
} else{
  print ("FALSE")
}

# Save the purity estimate output as a value to export to a text file
print_purity <- glue("improved_purity_estimate = ", purity)
print_rerun <- glue("rerun = ", rerun)

# Output the number that the neutral peak is offset by (first), and the purity estimate (second) into a text file.
capture.output(print_offset_mean, print_purity, print_rerun, file =args[2])
