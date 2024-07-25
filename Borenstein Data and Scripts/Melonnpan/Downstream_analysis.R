#Finding significant microbes
library(readxl)

healthy.microbes = read.csv("microbes_franzosa_IBD_2019_cr_healthy.csv",header = TRUE,row.names = 1)
diseased.microbes = read.csv("microbes_franzosa_IBD_2019_cr.csv",header = TRUE,row.names = 1)

microbes = row.names(healthy.microbes)
#compared to control
significant.microbes = data.frame(Median = rep('empty',length(microbes)),row.names = microbes)
medians.healthy = data.frame(Median = rep('empty',length(microbes)),row.names = microbes)

for(mic in microbes){
  
  control <- as.numeric(healthy.microbes[mic,])  
  diseased <- as.numeric(diseased.microbes[mic,])  
  
  # Combine data into a data frame
  data <- data.frame(
    group = factor(c(rep("Control", length(control)), rep("Diseased", length(diseased)))),
    value = c(control, diseased)
  )
  
  # Perform the Mann-Whitney U test
  result <- wilcox.test(value ~ group, data = data, exact = FALSE)
  medians.healthy[mic,'Median'] = median(data$value[data$group == "Control"])
  
  # Check if the test is significant
  if (result$p.value < 0.05) {

    # Determine the direction of the difference
    median_control <- median(data$value[data$group == "Control"])
    median_diseased <- median(data$value[data$group == "Diseased"])
    
    if (median_diseased > median_control) {
      significant.microbes[mic,'Median'] = 'Higher'
    } else {
      significant.microbes[mic,'Median'] = 'Lower'
    }
  } else {
    significant.microbes[mic,'Median'] = 'Same'
  }
}

write.csv(significant.microbes,file = "Significant_microbes_franzosa.csv")

#########################################################################
#Comparing with literature
#########################################################################

sig.microbes = read.csv("Significant_microbes_franzosa_lit.csv",header = TRUE,row.names = 1)
sig.microbes.excess = character()
sig.microbes.deficient = character()

for (mic in row.names(sig.microbes)) {
  
  sig.microbes.analysis = significant.microbes[mic,'Median']
  sig.microbes.literature = sig.microbes[mic,'Literature']
  is.same = grepl(sig.microbes.analysis, sig.microbes.literature)
  if(is.same){
    if(sig.microbes.analysis=="Higher"){
      sig.microbes.excess<-c(sig.microbes.excess,mic)
    }
    
    if(sig.microbes.analysis=="Lower"){
      sig.microbes.deficient<-c(sig.microbes.deficient,mic)
    }
  }
}

#########################################################################
#Finding microbial deficiencies in the predictions
#########################################################################

ci <- read.csv("Confidence_interval_franzosa.csv",header = TRUE,row.names = 1)
test.microbes = read.csv("Test_Microb_franzosa.csv",header = TRUE,row.names = 1)
pred.microbes = read.csv("Predicted_Microb_franzosa.csv",header = TRUE,row.names = 1)

microbes = row.names(ci)
more.microbes = character()
less.microbes = character()

for(mic in microbes){
  
  mic.abundace = pred.microbes[,mic]
  
  lower.limit = ci[mic,"LL"] 
  upper.limit = ci[mic,"UL"]
  
  # Comparing the abundance with ci
  ifelse(mic.abundace > upper.limit,more.microbes <- c(more.microbes,mic),ifelse(mic.abundace < lower.limit,less.microbes <- c(less.microbes,mic),NA))
}



#########################################################################
#Filtering the results with literature
#########################################################################

more.microbes.proven = more.microbes[more.microbes %in% sig.microbes.excess]
less.microbes.proven = less.microbes[less.microbes %in% sig.microbes.deficient]

#########################################################################
#Ranking them based on their normalized (with median) abundances
#########################################################################

pred.microbes.sig = pred.microbes

for (mic in colnames(pred.microbes)) {
  pred.microbes.sig[1,mic] = pred.microbes[1,mic]/as.numeric(medians.healthy[mic,'Median'])
}

pred.microbes.sig.less = pred.microbes.sig[colnames(pred.microbes)%in%less.microbes.proven]
pred.microbes.sig.more = pred.microbes.sig[colnames(pred.microbes)%in%more.microbes.proven]

order.more = rank(pred.microbes.sig.more)
order.less = rank(pred.microbes.sig.less)
pred.microbes.sig.more.ordered = pred.microbes.sig.more[,order(order.more)]
pred.microbes.sig.less.ordered = pred.microbes.sig.less[,order(order.less)]

#########################################################################
#Accuracy calculation. 3 classes
#########################################################################

library(caret)
microbes = row.names(ci)
pred.microbes.classes = pred.microbes
test.microbes.classes = test.microbes

for(mic in microbes){
  
  mic.abundace.test = test.microbes[,mic]
  mic.abundace.pred = pred.microbes[,mic]
  
  lower.limit = ci[mic,"LL"] 
  upper.limit = ci[mic,"UL"]
  
  # Comparing the abundance with ci
  ifelse(mic.abundace.test > upper.limit,test.microbes.classes[,mic] <- 'High',ifelse(mic.abundace.test < lower.limit,test.microbes.classes[,mic] <- 'Low',test.microbes.classes[,mic] <- 'In Range'))
  ifelse(mic.abundace.pred > upper.limit,pred.microbes.classes[,mic] <- 'High',ifelse(mic.abundace.pred < lower.limit,pred.microbes.classes[,mic] <- 'Low',pred.microbes.classes[,mic] <- 'In Range'))
}
pred.microbes.classes = as.data.frame(t(pred.microbes.classes))
test.microbes.classes = as.data.frame(t(test.microbes.classes))

test.patient = colnames(pred.microbes.classes)[1]

pred.microbes.classes[,1] <- as.factor(pred.microbes.classes[,1])
test.microbes.classes[,1] <- as.factor(test.microbes.classes[,1])

all_levels <- union(levels(pred.microbes.classes[,1]), levels(test.microbes.classes[,1]))
pred.microbes.classes[,1] <- factor(pred.microbes.classes[,1], levels = all_levels)
test.microbes.classes[,1] <- factor(test.microbes.classes[,1], levels = all_levels)

conf.matrix <- confusionMatrix(pred.microbes.classes[,1], test.microbes.classes[,1])

sink("Accuracy results.txt")
print(conf.matrix)
sink()

#########################################################################
#Probiotics
#########################################################################
food.recs = read_excel("food_version2_withImpact_CleanedAndGrouped.xlsx")
pred.microbes.sig.less.ordered.list = colnames(pred.microbes.sig.less.ordered)
pred.microbes.sig.more.ordered.list = colnames(pred.microbes.sig.more.ordered)

#subsetting only microbes into harmful and helpful
pred.microbes.sig.less.helpful = unique(food.recs[food.recs$Genus %in% pred.microbes.sig.less.ordered.list & food.recs$Impact=='Helpful',]$Genus)
pred.microbes.sig.less.harmful = unique(food.recs[food.recs$Genus %in% pred.microbes.sig.less.ordered.list & food.recs$Impact=='Harmful',]$Genus)
pred.microbes.sig.more.helpful = unique(food.recs[food.recs$Genus %in% pred.microbes.sig.more.ordered.list & food.recs$Impact=='Helpful',]$Genus)
pred.microbes.sig.more.harmful = unique(food.recs[food.recs$Genus %in% pred.microbes.sig.more.ordered.list & food.recs$Impact=='Harmful',]$Genus)

sink('Microbial Irregularities.txt')
cat('Microbial Abundaces Abnormally Low :\n')
print(pred.microbes.sig.less.ordered)
cat('\nMicrobial Abundaces Abnormally High :\n')
print(pred.microbes.sig.more.ordered)
cat('\nBeneficial Microbes :\n')
print(pred.microbes.sig.less.helpful)
print(pred.microbes.sig.more.helpful)
cat('\nHarmful Microbes :\n')
print(pred.microbes.sig.less.harmful)
print(pred.microbes.sig.more.harmful)
sink()

sink("Probiotic Recommendation.txt")
cat('You are deficient in the following microbial genera :\n')
ifelse(length(pred.microbes.sig.less.helpful)==0, 'NA',print(colnames(pred.microbes.sig.less.helpful)))
sink()

#########################################################################
#Food Recommendations
#########################################################################

#Foods to intake. Increases growth of beneficial low level microbes and Decreases growth of bad abundant bacteria
#Foods to limit. Decreases growth of beneficial low level microbes and Increases growth of bad abundant bacteria

sink("Food Recommendation.txt")
cat('Add these foods to your diet :\n')
unique(food.recs[food.recs$Genus %in% pred.microbes.sig.less.helpful & food.recs$Effect=='Increase',]$`Examples to be included in report`)
unique(food.recs[food.recs$Genus %in% pred.microbes.sig.more.harmful & food.recs$Effect=='Decrease',]$`Examples to be included in report`)

cat('\nLimit these foods from your diet :\n')
unique(food.recs[food.recs$Genus %in% pred.microbes.sig.less.helpful & food.recs$Effect=='Decrease',]$`Examples to be included in report`)
unique(food.recs[food.recs$Genus %in% pred.microbes.sig.more.harmful & food.recs$Effect=='Increase',]$`Examples to be included in report`)
sink()

#########################################################################
#Phyla level
#########################################################################
medians.subset = medians.healthy
medians.subset$rownames = row.names(medians.subset)
medians.subset = medians.subset[medians.subset$rownames %in% sig.microbes.excess|medians.subset$rownames %in% sig.microbes.deficient ,]


pred.microbes.subset = pred.microbes
pred.microbes.subset = as.data.frame(t(pred.microbes.subset))
pred.microbes.subset$rownames = row.names(pred.microbes.subset)
pred.microbes.subset = pred.microbes.subset[pred.microbes.subset$rownames %in% sig.microbes.excess|pred.microbes.subset$rownames %in% sig.microbes.deficient ,]

