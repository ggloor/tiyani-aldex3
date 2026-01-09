## This is intended to be the code repository for Pranav's Biochem 3383 project

 The idea is to analyze the Tiyani data with age as a continuous variable with ALDEx3 linear model and to include scale uncertainty

Task 1 is to get ALDEx2 and ALDEx3 working on the selex dataset

- download ALDEx2 from Bioconductor

- download ALDEx3 from https://github.com/jsilve24/ALDEx3


 ```{r}   
  if (!require("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
devtools::install_github("jsilve24/ALDEx3")
```
## - Pranav's edits

- Downloaded R and R-studio
- Got a text editor setup and working
- Finished the terminal setup and functional

## Currently working on how to download ADLEX2 through Bioconductor

# this is command from the the Bioconductor - Will run it in the terminal
``` {r}
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
    }
  BiocManager::install("ALDEx2")

  ```
## Finished downloading the ALDEx2 and BiocManager

Testing around and was following a video to get some feel to use R studio and how to use a data and set and other concepts in one of the paper

-----------------

## Now Moving on the downloading packages such as "tidyverse" and other packages
  # ggplot2, ggpattern, cowplot

Now that I have ADLEx2 and other packages downloaded and I will move on to the test analysis of sim_seq_data from a Michelle Nixon

    "https://github.com/michellepistner/DAWG_workshop/blob/main/scripts/1.0_sim_analysis.R"

 - I will provide the code and my understanding of the code below:

 - ``` {r}
## Analysis code for sim_seq_data

# importing the libraries

    library(ALDEx2)
    library(tidyverse)
    library(ggplot2)
    library(ggpattern)
    library(cowplot)

# To ensure reproducibility
    set.seed(12345)

# Reading in Data  
  rdat <- read.cvs(file.path("/Users/pranavdivvela/Desktop/3383/0_git/Data_files/TEST_DATA/sim_seq_dat.csv"))

  ##Reading in the simulated flow data for bulding the scale model
flow_data <- read.csv(file.path("/Users/pranavdivvela/Desktop/3383/0_git/Data_files/TEST_DATA/sim_seq_dat.csv"))

## Inspecting elements

  ## "Y" represents the OTU table
Y <- t(rdat[,-1])

## Vector denoting whether samples was in pre- or post- antibiotic condition.
conds <- as.character(rdat[,1])

## Fitting and the analyzing the orginal ALDEx2 model
mod.base <- aldex(Y,conds) ## gamma=NULL meaning there is no level of unceritainty in the analysis
mod.base %>% filter(we.eBH < 0.05)

## Recreating ALDEx2
mod.clr <- aldex(Y,conds, gamma = 1e-3)
## even though the function was built to add uncertainty the value we assigned it so small that it doesn't affect anything
 mod.clr %>% filter (we.eBH < 0.05)

## Checking for concordance in effect sizes
plot(mod.base$effect, mod.clr$effect, xlab = "Original ALDEx2
  Effect Size", ylab = "CLR Scale Model Effect Size")
abline(a=0,b=1, col = "red", lty = "dashed")
