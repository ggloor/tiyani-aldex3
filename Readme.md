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

## Now Moving on the downloading packages such as "tidyverse" and other packages
  # ggplot2, ggpattern, cowplot

Now that I have ADLEx2 and other packages downloaded and I will move on to the test analysis of sim_seq_data from a Michelle Nixon

    "https://github.com/michellepistner/DAWG_workshop/blob/main/scripts/1.0_sim_analysis.R"

 - I will provide the code and my understanding of the code below:

 - 
