# Load necessary library
library(MASS)
library(readr)
library(boot)

#import the measured and predicted microbe abundances
micro_data <- read.csv("microbes_franzosa_IBD_2019_cr_healthy.csv",header = TRUE,row.names = 1)

microbes = row.names(micro_data)
samples = colnames(micro_data)

# Define a function to calculate the 2.5th and 97.5th percentile
percentiles_2.5 <- function(data, indices) {
  return(quantile(data[indices],0.025))
}

# Define a function to calculate the 2.5th and 97.5th percentile
percentiles_97.5 <- function(data, indices) {
  return(quantile(data[indices],0.975))
}

ci = data.frame(LL = rep(0,length(microbes)),UL = rep(0,length(microbes)),row.names = microbes)

for(mic in microbes){
  mic.abundace = as.numeric(micro_data[mic,])
  
  # Perform bootstrapping
  set.seed(123)  # For reproducibility
  bootstrap_results2.5 <- boot(mic.abundace, percentiles_2.5, R = 1000)
  bootstrap_results97.5 <- boot(mic.abundace, percentiles_97.5, R = 1000)
  
  # Original statistic
  original_stat2.5 <- bootstrap_results2.5$t0
  original_stat97.5 <- bootstrap_results97.5$t0
  
  ci[mic,"LL"] = original_stat2.5
  ci[mic,"UL"] = original_stat97.5
}

write.csv(ci,file = "Confidence_interval_franzosa.csv")


###Other formulas
# Perform bootstrapping
set.seed(123)  # For reproducibility
bootstrap_results2.5 <- boot(data, percentiles_2.5, R = 1000)
bootstrap_results97.5 <- boot(data, percentiles_97.5, R = 1000)

# Calculate the 95% confidence interval
ci2.5 <- boot.ci(bootstrap_results2.5, type = "norm")
ci97.5 <- boot.ci(bootstrap_results97.5, type = "norm")

# Original statistic
original_stat2.5 <- bootstrap_results2.5$t0
original_stat97.5 <- bootstrap_results97.5$t0

# Bias
bias2.5 <- mean(bootstrap_results2.5$t) - bootstrap_results2.5$t0
bias97.5 <- mean(bootstrap_results97.5$t) - bootstrap_results97.5$t0

# Standard Error
std_error2.5 <- sd(bootstrap_results2.5$t)
std_error97.5 <- sd(bootstrap_results97.5$t)



