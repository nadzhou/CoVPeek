# CoVPeek - proteome divergence visualzation tool
(Yes, I'm looking for contributions. I think I'm un to something important and would 
appreciate others' input on this.)

This tool visualizes divergence divergence amongstt different global strains of SARS CoV2. 


## Motivation 
While genome sequence visualization is done by Nextstrain, we do not yet have a proteome level variation among different viral sequences. This tool aspires to find variation per amino acid position. This is done by comparing the known Uniprot protein sequence with the available data from GISAID. 

## Goals 
Thiis tool will do the following: 
  1. Visualize the divergence on a bar chart. The bar chart 
  will contain by country variation. 
  1. Write a PyMol script to see the changes that have been 
  documented with respect to the canonical Uniprot sequence. 


# Update 
## Incorporating Jupyter Notebook 
Showing the code in real time is important, so I have added some notebooks to go through the code step by step. Here, I show that the spike protein does indeed show divergence at the 613th position (starting from 0). 

![image](https://user-images.githubusercontent.com/25282805/82142814-8db1b880-9858-11ea-92f1-b433f38b8dd0.png)

This graph clearly shows there is a significant mutation rate at the 613th position, while the others either contain an amino acid or a gap '-' while this specific position contains, as is evident in the notebook, there is occurence of all amino acids in that position. 
