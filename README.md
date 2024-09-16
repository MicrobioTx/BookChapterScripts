Folder 1. Data and Scripts

Both the sub folders - MiMeNet and Melonnpan have the same data and the same scripts. Only the file Predicted_Microb_franzosa.csv is changed in each folder. This file is the predicted microbe abundance by the respective model (MiMeNet/ Melonnpan).

The scripts -

1. Data.R
    Script to retrieve the franzosa dataset and metadata from the Borenstein Data repository

2. Melonnpan/Downstream_analysis.R

All downstream statistically analysis on the results

3. Melonnpan/Non-parametric.R

Script to find reference intervals for microbe abundances from healthy populations

Folder 2. Melonnpan Scripts

All scripts within the melonnpan folder are unchanged from the original Melonnpan scripts. The main script however has been reversed and curated towards the franzosa dataset.

The scripts-

1. Melonnpan-Fransoza

Reversed Melonnpan script to predict microbes from metabolites on the franzosa dataset

Folder 3. MiMeNet Scripts

The scripts -

1. MiMeNet_train_Borenstein.py

Reversed Mimenet script to predict microbes from metabolites specifically for the test patient 'Validation.UMCG8089650' in the franzosa dataset.
