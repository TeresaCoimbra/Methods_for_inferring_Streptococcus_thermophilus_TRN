####### Script used to obtain the final dataset ######
# Teresa Coimbra

## Final Dataset includes: 
# normalized_df.csv - containing only the normalized readcounts after SARTools
# fixed_metadata.csv - containing the metadata regarding all the datasets

# Dataset after SARTools:
setwd("")                                                                # insert path
data <- read.delim("ERR4C1vsSRR5CD2.complete.txt", row.names=1)
metadata <- read.csv("metadata.csv",row.names=1, sep=";")
dim(data)

#A new dataframe was created containing only the columns with the normalized readcounts:
data_norm= data[ ,grepl('norm', names(data))]

# Dimensions of the new dataset (confirms that all the important information is retained):
dim(data_norm)


#### DATA:
# compare if it would be the same getting the columns from a different 'complete' file:
data2 <- read.delim("SRR5CD1vsSRR3CHCC5086.complete.txt", row.names=1)
data_norm2= data2[ ,grepl('norm', names(data2))]
all.equal(data_norm2,data_norm)

colnames(data_norm) <- gsub('norm.','',colnames(data_norm))  # remove the 'norm' from the column names
which(colnames(data_norm) == "SRR390316_1")                  # which column is "SRR390316_1"
which(colnames(data_norm) == "SRR390317_1")
colnames(data_norm)[16] <- 'SRR390317'                       # change the names
colnames(data_norm)[17] <- 'SRR390316'
df = data_norm[,order(colnames(data_norm))]                  # obtain a new df which has the ordered columns   
head(df)

write.csv(df,'normalized_df.csv')                           # create csv file for normalized readcounts

#### METADATA: 
metadata = metadata[order(rownames(metadata)),]              # ordering the row names 
#head(metadata)
#dim(metadata)
# metadata %>%
#   mutate(exp_or_stat = case_when(
#     growth_phase == 'Late-exponential' ~ "Exponential",
#     growth_phase == 'Exponential' ~ "Exponential",
#     growth_phase == 'Stationary ' ~ "Stationary"
#   ))

write.csv(metadata,'fixed_metadata.csv')                    # create csv file for the metadata
