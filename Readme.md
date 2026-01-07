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