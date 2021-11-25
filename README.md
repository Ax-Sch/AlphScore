# AlphScore
Improving pathogenicity prediction of missense variants by using AlphaFold-derived features

This Repository is structured into two parts: 
- extract features
Extract features from Alphafold2 derived structures ( https://alphafold.ebi.ac.uk/ ) using the tools DSSP, FEATURE, protinter and biop... . Finally these Features are combined with data from dbNSFP 4.2.
![alt text](https://github.com/Ax-Sch/AlphScore/blob/main/Overview.png?raw=true)


- Analysis
Here the extracted and merged features are analysed and used for prediction of pathogenicity. This part is written in R.
