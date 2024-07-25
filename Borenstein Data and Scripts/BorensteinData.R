#Loading Libraries
library(ggplot2)
library(viridis)
library(broom)
library(dplyr)
library(gt)
library(formula.tools)
library(logger)
library(future.apply)
library(meta)
library(kableExtra)
library(tibble)

# Notebook settings
future::plan("multisession", workers = 4)
options(scipen = 999)

# Load utility scripts
source("scripts/data_organization/utils.R")
source("scripts/data_analysis/hmdb_utils.R")

# Load all data available in the curated gut microbiome-metabolome data resource
all.data <- load.all.datasets(parent.folder = "processed_data")
for(i in 1:length(all.data)) assign(names(all.data)[i], all.data[[i]])
rm(all.data)

#Loading metabolites
metab = mtb[["FRANZOSA_IBD_2019"]]
metabmap = mtb.map[["FRANZOSA_IBD_2019"]]
row.names(metab) = metab$Sample
metab = as.data.frame(t(metab[,-1]))

# Convert row names to a column
metab <- rownames_to_column(metab, var = "Metabolite")

# Filter out all metabolites with no name or low confidence annotations
metabmap_filtered = metabmap[metabmap$High.Confidence.Annotation & !is.na(metabmap$Compound.Name),]

# Replace the current name of each metabolite with that from metab.map
matched_indices <- match(metab$Metabolite, metabmap_filtered$Compound)
metab$Metabolite <- metabmap$Compound.Name[matched_indices]
metab_filtered = metab[!is.na(metab$Metabolite),]

# Group by Metabolite name and calculate the mean for each sample
metab_abundance <- metab_filtered %>% group_by(Metabolite) %>% summarise(across(everything(),mean))
rownamesmetab_abundance
  
#Loading microbes
microb = genera[["FRANZOSA_IBD_2019"]]
row.names(microb) = microb$Sample
microb = as.data.frame(t(microb[,-1]))

# Extract genus names
extract_genus_name =  function(row.name){
  #Split at semicolons
  taxas = strsplit(row.name,";")[[1]]
  #Find the genus
  genus = taxas[grep("^g_",taxas)]
  #Remove the prefix
  genus.name = sub("g__","",genus)
  
  return(genus.name)
}

# Replacing the nomenclature with only genus names
genus.names = sapply(row.names(microb),extract_genus_name)
row.names(microb) = genus.names

# Removing control patients
meta.data = metadata[["FRANZOSA_IBD_2019"]]
control.patients = meta.data$Study.Group=='Control'
microb_healthy = microb[,control.patients]

cr <- read.csv("C:/Users/USER/Desktop/Borenstein/Working File Genera Report.csv",header = TRUE)
microb_healthy_cr = microb_healthy[row.names(microb_healthy)%in%cr$Genus.Name,]
write.csv(microb_healthy_cr,file = 'microbes_franzosa_IBD_2019_cr_healthy.csv',row.names = TRUE)


