# Analysis of tissue samples stained with markers to delineate EMT tumor heterogeneity and phenotypes

Full methods and results can be found in the corresponding:

- methods paper: https://doi.org/10.1126/sciadv.abj8002
- results paper: https://doi.org/10.1126/sciadv.abj8002

## Brief overview of staining method 

Tyramide Signal Amplification (TSA) is a multiround, multiplexed method of staining FFPE tissue with high specificity for multiple overlapping markers. It uses enzyme-mediated catalysis of HRP-conjugated secondary antibodies to deposit a labeled tyramide or OPAL fluorophore at the site of primary antibody detection, which remains in place after primary and secondary antibodies are stripped away by a microwaving step. The order and concentration of primary and secondary antibodies, as well as pairing with each OPAL fluorophore must be carefully considered based on overall intensity of protein expression and detectability of epitopes. 

![](https://ars.els-cdn.com/content/image/1-s2.0-S0091679X22000796-f08-01-9780323900188.jpg)

1. Human tumor samples are fixed and prepared before being stained for six critical markers using a multiplexed Tyramide Signal Amplification method described in methods paper above
* EMT markers
  + Transcription Factors: Snail & ZEB1
  + Intermediate Filaments: K8, K14, VIM
  + Cell Adherens Junction: E-Cadherin
2.  Stained tissue slides are imaged and Regions of Interest (ROIs) are chosen to represent and encompass the totality of the tumor
3.  Imaging analysis conducted using Perkin Elmer software to segment cells as well as assign and train algorithm to detect cellular phenotypes based on co-expression and localization of markers
4. Cellular locations, phenotypes, and marker expression exported for futher analysis

# Analysis 

`TMA_TSA_Analysis_annotated.R` 
- annotates inform file outputs with sample IDs
- Formatted for use with ROIs from a tumor slice or from Tumor Microarray punches
- Checks mean counts to confirm good overall staining. Can be done BEFORE phenotyping
- Compiles Total normalized Counts from the entire cell for each sample. This can be altered if there is highly specific cell localization i.e. choosing "Nucleus" or "Cellular Membrane" for specific markers
- Removes cells that have not been phenotyped and tidies up data and normalize markers to a percentile 0-1 for violin plots

_Note: TSA_Analysis_annotated.R contains an earlier version of the code but more or less does the same things_

### Generation of heterogeneity scores can be found in the following github location 
https://github.com/drmsbrown/cell-heterogeneity-emtscore

`TMA_TSA_EMT&Het_annotated.R` provides code to integrate EMT and heterogeneity scores from Machine learning algorithm and provides plots 

![](https://ars.els-cdn.com/content/image/1-s2.0-S0091679X22000796-f08-02-9780323900188.jpg)
