
library(dplyr)
library(readr)
library(stringr)

 matrix_all <- read_table2("matrix.txt", col_names = F)

## DEF function for selecting kmers of proper species
replace_correct3  <- function(df,i){
     
  df[[i]]  <- if_else(mapply(grepl, paste0("^",df$X5,":"), df[[i]]), # creates TRUE/FALSE column saying whether the final species matches to the species of the k-mer
                       str_replace(df[[i]], paste0(df$X5,":"),""), # if it matches, we will keep the number of alighed k-mers
                      NA_character_) # if it doesn't match, we will change it to NA
 return(df)
}

## Run matrix transformation 

for(i in c(7:ncol(matrix_all))){ # loop over all columns that include k-mer alignment information
    matrix_all  <- replace_correct3(matrix_all,i) 
}

matrix_all  <- matrix_all %>% mutate_at(vars(-X1,-X2,-X3,-X4,-X5,-X6), as.numeric) 
matrix_all  <- matrix_all %>% 
    mutate("sum" = rowSums((matrix_all[,7:ncol(matrix_all)]), na.rm = TRUE))  %>% 
    mutate(pct_kmers = sum/(X6-34))  %>% select(barcode = X2, species = X5, pct_kmers)

write.table(matrix_all, "matrix_all_with_scores.txt", row.names = F)