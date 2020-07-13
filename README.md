# TS-COL\_0137-Vathiotis\_Rimm

### Description

Scripts related to IO 360 and GeoMx DSP analysis of melanoma samples in collaboration with Yale through the Rimm lab lead by Ioannis Vathiotis. The files contained in this document serve as reference during and after review of the related manuscript: *Models that Combine Transcriptomic with Spatial Protein Information Exceed the Predictive Value for Either Single Modality*

**Authors:** Zhi Yang, Jason Reeves

**Author Affiliations:** NanoString Technologies Inc

**Files contained:**

- 01\_TSCOL\_0137\_Initialize.R: script to read in and process univariate analysis of the GeoMx and IO 360 data. Outputs figures and a working environment for downstream multivariate analysis via regularization.
- 02\_TSCOL\_0137\_Building\_Models.R: script to build and test multivariate regularization models based on GeoMx and IO 360 data.
- data\/graphing\_adjusted\_log2\_expression\_data\_WITH DSP.xlsx: combined GeoMx and IO 360 expression results provided to NanoString for analysis.
- data\/TargetsNotInProtein.csv: lookup table for targets which have different gene symbols than protein target IDs
- output\/permutations/\*.rdata: files related to permutation testing based run during analysis for ease of access. Files can be generated on the fly, but these are provided to reduce run time to test specific cases if necessary.

### Software Requirements
The code has been test on the following systems:
Linux: x86_64-pc-linux-gnu (64-bit)
Windows: x86_64-w64-mingw32/x64 (64-bit)

Users need to install `R` with version of 3.6.3 or higher and the following R packages before running the script.   

```
install.packages(c('readxl', 'ggplot2', 'ggrepel', 'pheatmap', 'survival', 
'glmnet', 'DT', 'pROC', 'survminer', 'patchwork', 'glmnetUtils', 
'OptimalCutpoints'', 'caret', 'kernelboot'))
```

However, for the `rsq` package, please make sure you download the `1.1` version from the website here https://cran.r-project.org/src/contrib/Archive/rsq/ followed by installing it using `install.packages("../rsq_1.1.tar.gz", repos = NULL, type = "source")`. 

Below lists the version of all R packages mentioned above.   
```
kernelboot            caret          lattice OptimalCutpoints      glmnetUtils        patchwork 
   "0.1.7"         "6.0-86"        "0.20-41"          "1.1-4"          "1.1.5"          "1.0.0" 
 survminer           ggpubr         magrittr             pROC               DT           glmnet 
   "0.4.6"          "0.2.5"            "1.5"         "1.16.2"           "0.13"          "3.0-2" 
    Matrix         survival         pheatmap          ggrepel          ggplot2           readxl 
  "1.2-18"          "3.1-8"         "1.0.12"          "0.8.2"          "3.3.0"          "1.3.1" 
       rsq 
     "1.1" 
```
 

**Contact us:**\
NanoString Technologies, Inc.\
530 Fairview Avenue N\
Seattle, WA 98109\
Tel: (888) 358-6266\
zyang@nanostring.com

### License
Copyright (C) 2020, NanoString Technologies, Inc.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see https://www.gnu.org/licenses/.
