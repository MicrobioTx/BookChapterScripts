# Loading all required libraries
library(melonnpan)
library(readxl)

# Changes were made in the original scripts. We use the source function to call these modified functions to the global environemnt
source("C:/Users/USER/Desktop/Melonnpan/melonnpan/R/melonnpan_train.R")
source("C:/Users/USER/Desktop/Melonnpan/melonnpan/R/melonnpan_predict.R")

# Creating an output file
test1 = "C:/Users/USER/Desktop/Melonnpan/Franzosa"

# acquiring data. Data from all three datasets have been merged
fz.microb.data = read.csv("microbes_franzosa_IBD_2019_cr.csv",row.names = 1)
fz.metab.data = read.csv("metabolites_franzosa_IBD_2019.csv",row.names = 1)
#totalmicrobes = nrow(read.csv("New_RAW-Merged-GenusData.csv"))

# changing the row names to microbe/metabolite name
# Make the first column as row names
#rownames(mtx.metab.data) <- mtx.metab.data$Metabolite
#row.names(mtx.microb.data) <- mtx.microb.data$Genus

# Remove the first column after setting row names
#mtx.microb.data = mtx.microb.data[,-1]
#mtx.metab.data = mtx.metab.data[,-1]

# pre-processing of data
#NA microbial abundances are set to 0
fz.microb.data[is.na(fz.microb.data)] = 0 

#Function to replace NA metabolite concentrations with median value accross all individuals
replace_na_with_median <- function(row) {
  row[is.na(row)] <- median(row, na.rm = TRUE)
  return(row)}

#Apply the function to rows of metabolite dataset
fz.metab.data <- t(apply(fz.metab.data, 1, replace_na_with_median))

#Total Sum Scaling (TSS) normalization. 
column.sum.met <- colSums(fz.metab.data)
fz.metab.data = as.data.frame(sweep(fz.metab.data, 2, column.sum.met, `/`))

#column.sum.mic <- colSums(mtx.microb.data)
#mtx.microb.data = as.data.frame(sweep(mtx.microb.data, 2, column.sum.mic, `/`))

#We set the threshold as 1e-4. Any value below this, either measured or predicted is set to 0.
#This is done for microbe data, for accuracy calculations
#fz.microb.data[fz.microb.data<1e-4] = 0

#Transposing both dataframes
fz.metab.data = as.data.frame(t(fz.metab.data))
fz.microb.data = as.data.frame(t(fz.microb.data))

# Determining the threshold for removing low abundance microbes

# Function to calculate binomial coefficient (nCr)
nCr <- function(n, r) {
  factorial(n) / (factorial(r) * factorial(n - r))
}

# Function to calculate probability of the train data in each fold of CV having only 0
#Calculations have been done based on calculation that 146 points will be part of the train set
prob = function(row){
  total_points <- length(row)
  zero_points <- length(row[row==0])
  nonzero_points <- total_points - zero_points
  selected_points <- 146
  
  # Check if there are at least 50 zero points
  if (zero_points < selected_points) {
    return(0)
  }
  
  
  probs = nCr(zero_points, selected_points) * nCr(nonzero_points, 0) / nCr(total_points, selected_points)
  return(probs)
}


# Calculate the probability
probability.micro <- sapply(fz.microb.data, prob)


#If we subset only those microbes with less than 0.05 probability for having all 0 values as output in CV
fz.microb.data <- fz.microb.data[,probability.micro < 0.005]

# transforming datasets. rntransformation for microbe abundance and ArcSin transformation for metabolite abundances 
#mtx.microb.data = as.data.frame(apply(mtx.microb.data,2,rntransform))
#row.names(mtx.microb.data) = row.names(mtx.metab.data)
#mtx.metab.data = as.data.frame(apply(mtx.metab.data,2,ArcSin))

# creating test and train
set.seed(123)  # For reproducibility
Test = sample(nrow(fz.metab.data), size = 1)
test.metab = fz.metab.data[Test,]
test.metag = fz.microb.data[Test,]
train.metab = fz.metab.data[-Test,]
train.metag = fz.microb.data[-Test,]

# This indentifies any column/microbe that has less than 3 unique values. Ideally it should be 2. But CV will again split the training set.
test.uniqueness = train.metag[,sapply(train.metag, function(x) length(unique(x)))< 3]
if(length(test.uniqueness)!=0){print('Dont train model, Some microbes in the training set have constant values')}

# training the model. Data is already transformed
output = melonnpan.train(metab = train.metag,metag = train.metab,test1,no.transform.metab = TRUE,no.transform.metag = TRUE,discard.poor.predictions = TRUE)


# Read the weights text file into R
weights_table = output

# Transpose the data frame
weights_table_transposed <- t(weights_table)

# Convert the transposed matrix back to a data frame if necessary
weights <- as.data.frame(weights_table_transposed)

# Predict the microbes from the test metabolites
pred = melonnpan.predict(test.metab,output = test1, weight.matrix = weights,train.metag = train.metab,no.transform.metab = TRUE,no.transform.metag = TRUE)

# Retrieving the predicted microbe data for the test sample
predicted.metag = pred

#We set the threshold as 1e-4. Any value below this, either measured or predicted is set to 0.
#This is done for microbe data, for accuracy calculations
#predicted.metag[predicted.metag<1e-4] = 0

write.csv(predicted.metag, file = file.path(test1, 'Predicted_Microb.csv'), row.names = TRUE)
write.csv(test.metag, file = file.path(test1, 'Test_Microb.csv'), row.names = TRUE)


# accuracy
# We combine the entire data to get confidence intervals
data = as.data.frame(rbind(train.metag,test.metag))

# computing confidence intervals

# Function to calculate 2.5th and 97.5th percentile of the data
percentile_ci <- function(data) {
  
  nonzero.data = data[data!=0]
  
  percentiles <- quantile(nonzero.data, probs = c(0.025, 0.975))
  
  ci_lower <- as.numeric(percentiles[1])
  ci_upper <- as.numeric(percentiles[2])
  
  # Return the confidence interval as a vector
  return(c(ci_lower,ci_upper))
}

# Apply the function to each column
ci.results <- apply(data,2, percentile_ci)

# Convert results to a data frame
ci.metag <- as.data.frame(t(ci.results))
names(ci.metag) <- c("Lower Bound", "Upper Bound")
ci.metag = t(ci.metag)

write.csv(ci.metag, file = file.path(test1, 'Confidence_Intervals.csv'), row.names = TRUE)

# accuracy
predicted.genera = intersect(colnames(ci.metag), colnames(test.metag))

# Create a vector of genera
genera_idx <- colnames(predicted.metag)

#number of patients
num_genera = length(predicted.metag[1,])

# Create an empty dataframe to collect TN,TP,FN and FP
con.table <- data.frame(TP=integer(num_genera), 
                        FP=integer(num_genera), 
                        TN=integer(num_genera), 
                        FN=integer(num_genera),
                        Accuracy=integer(num_genera),
                        stringsAsFactors=FALSE,
                        row.names = genera_idx)

for(genus in row.names(con.table)){
  for(patient in rownames(predicted.metag)){
    measured = test.metag[patient,genus]
    prediction = predicted.metag[patient,genus]
    low.int = ci.metag["Lower Bound",genus]
    up.int = ci.metag["Upper Bound",genus]
    
    if(measured == 0){
      ifelse(prediction==0,con.table[genus,"TN"] <- con.table[genus,"TN"] + 1,con.table[genus,"FP"] <- con.table[genus,"FP"] + 1)
    }
    
    if(measured != 0){
      ifelse(low.int<=prediction & prediction<=up.int,con.table[genus,"TP"] <- con.table[genus,"TP"] + 1,con.table[genus,"FN"] <- con.table[genus,"FN"] + 1)
    }
  }
  con.table[genus,"Accuracy"] = (con.table[genus,"TP"]+con.table[genus,"TN"])/(con.table[genus,"TP"]+con.table[genus,"TN"]+con.table[genus,"FP"]+con.table[genus,"FN"])
}

write.csv(con.table, file = file.path(test1, 'Confusion_Matrix.csv'), row.names = TRUE)

tp = sum(con.table[,"TP"])
tn = sum(con.table[,"TN"])
fp = sum(con.table[,"FP"])
fn = sum(con.table[,"FN"])

# Total accuracy for only modeled microbes
Accuracy.total.modeled = (tp+tn)/(tp+tn+fp+fn)
#Accuracy.total = (tp+tn)/(tp+tn+fp+fn+nrow(test.metag)*(totalmicrobes-nrow(con.table)))

print(paste0('Prediction Accuracy for only modelled microbes is ',round(Accuracy.total.modeled,2)))

# Number of microbes predicted with more than 90 percentage accuracy
high.acc = con.table[,"Accuracy"]>=0.90
print(paste0('Total Microbes predicted with more then 90% accuracy is ',sum(high.acc)))
print(con.table[high.acc,])

Sensitivity = tp/(tp+fn)
Specificity = tn/(tn+fp)
print(paste0('Prediction Sensitivity of modelled microbes is ',round(Sensitivity,2)))
print(paste0('Prediction Specificity of modelled microbes is ',round(Specificity,2)))

print(paste0('Prediction Accuracy considering all microbes is ',round(Accuracy.total,2)))

# Comparing with list of clinically relevant microbes
#library(tibble)
#library(dplyr)
#Clinically.relevant.microbes = read_excel("Working File Genera Report.xlsx")
#Clinically.relevant.microbes = as.character(Clinically.relevant.microbes %>% pull(`Genus Name`))
#Clinically.relevant.microbes = unique(Clinically.relevant.microbes)

#Finding the clinically relevant microbes for which model was generated
#con.table.clinical.relevant = con.table[intersect(row.names(con.table),Clinically.relevant.microbes),]

# Printing the clinically relevant microbes that have not been modelled
#setdiff(Clinically.relevant.microbes,row.names(con.table.clinical.relevant))
#not.modelled.cr = setdiff(Clinically.relevant.microbes,row.names(con.table.clinical.relevant))

# Accuracy calculations
#tp.cr = sum(con.table.clinical.relevant[,"TP"])
#tn.cr = sum(con.table.clinical.relevant[,"TN"])
#fp.cr = sum(con.table.clinical.relevant[,"FP"])
#fn.cr = sum(con.table.clinical.relevant[,"FN"])


#Accuracy.total.modeled.cr = (tp.cr+tn.cr)/(tp.cr+tn.cr+fp.cr+fn.cr)
#Accuracy.total.cr = (tp.cr+tn.cr)/(tp.cr+tn.cr+fp.cr+fn.cr+length(not.modelled.cr))

# Number of microbes predicted with more than 90 accuracy
#high.acc.cr = con.table.clinical.relevant[,"Accuracy"]>=0.90
#print(paste0('Prediction Accuracy for modelled clinically relevant microbes is ',round(Accuracy.total.modeled.cr,2)))

#print(paste0('Total Clinically relevant Microbes predicted with more then 90% accuracy is ',sum(high.acc.cr)))
#print(con.table.clinical.relevant[high.acc.cr,])


#Sensitivity.cr = tp.cr/(tp.cr+fn.cr)
#Specificity.cr = tn.cr/(tn.cr+fp.cr)
#print(paste0('Prediction Sensitivity of modelled clinically relevant microbes is ',round(Sensitivity,2)))
#print(paste0('Prediction Specificity of modelled clinically relevant microbes is ',round(Specificity,2)))

#print(paste0('Prediction Accuracy considering all clinical relevant microbes is ',round(Accuracy.total.cr,2)))
