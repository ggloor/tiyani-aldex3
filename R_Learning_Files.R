#basic math calculation
5+5

#assigning values to variables
a<- 11
b<-9

#Commands for the variables
sum(a,b)

# VECTOR: an ordered collection of values of the SAME TYPE
  # Create with c(): ages <- c(25,35,45,55)
ages <-c(25,35,45,55)
  # Pull data by position using brackets:
  #   ages[1] = first value (25)
  #   ages[3] = third value (67)
  #   ages[2:4] = values 2 through 4 (35,45,55)
  #   ages[c(1,3)] = first and third values (25, 45)


names <- c("Pranav","Heer","Vibs")
age <- c(25,35,45)
gender <- c("M","F","F")

# ACCESSING DATAFRAME DATA
# df[2] = column 2 as a dataframe
# df[[2]] = column 2 as a vector
# df$column_name = column as a vector
# df[,2] = column 2 as a vector
# df[2,] = row 2
# df[2,3] = cell at row 2, column 3
friends <- data.frame(names,age,gender)

library(tidyverse)

friends %>% 
  select(names,age) %>%
  filter(age>30) %>%
  arrange(age)
  
