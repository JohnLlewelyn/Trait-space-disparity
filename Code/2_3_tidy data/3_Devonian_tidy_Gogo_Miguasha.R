#combine RDSs (after adding names) and check/tidy data
library(ggplot2)
library(dplyr)
library(readxl)

#Adjust file path at"###"

# Set working directory to  folder containing the RDS files
setwd("###/")

#new_env1 <- new.env()
DF <- readRDS("data/Devonian_fish_traits_RDS/Mig&GoImp_2024.RDS")

#subset to communities/sites to keep
kp <- c("Miguasha","Gogo")
mg <- DF
mg <- mg[mg$site%in%kp,]

#fix names so they match those in modern fish data set, add missing columns
names(mg)[names(mg)=="taxon"] <- "Species"
names(mg)[names(mg)=="total length"] <- "TL"
names(mg)[names(mg)=="standard length"] <- "SL"
names(mg)[names(mg)=="head length"] <- "HL"
names(mg)[names(mg)=="pre-orbital length (snout length)"] <- "POL"
names(mg)[names(mg)=="body depth (body height)"] <- "BD"
names(mg)[names(mg)=="eye size (diameter)"] <- "ED"
names(mg)[names(mg)=="body shape I (saggital plane)"] <- "BodyShapeI"
names(mg)[names(mg)=="body shape II (transverse plane)"] <- "BodyShapeII"
names(mg)[names(mg)=="mouth position"] <- "PosofMouth"
names(mg)[names(mg)=="eye position"] <- "eye.position"
names(mg)[names(mg)=="spiracular"] <- "spiracle"
names(mg)[names(mg)=="caudal fin shape (koaw.org)"] <- "caudal.fin.shape"
names(mg)[names(mg) == "site"] <- "community" 

#add some columns
mg$habitat = ifelse(mg$community=="Miguasha","estuary", ifelse(mg$community=="Gogo","reef","don't know"))
mg$degsFromEquator = "?"

#cut to columns to keep
trait_names <- c("Species","BodyShapeI","BodyShapeII","TL","HL","ED","POL","BD","PosofMouth","mandible","eye.position","spiracle","caudal.fin.shape","community","habitat","degsFromEquator")
mg <- mg[,trait_names]

#change numeric columns to numeric
mg[] <- lapply(mg, function(x) {
  if(all(grepl("^-?\\d+\\.?\\d*$", x))) return(as.numeric(x))
  else return(x)
})

#change Devonian morphometrics (apart from total length [TL]) to be fraction of TL
# Identify the numeric columns (excluding "TL")
numeric_columns <- sapply(mg, is.numeric)
numeric_columns <- numeric_columns & !colnames(mg) == "TL"

# Divide the numeric columns by "TL"
mg[, numeric_columns] <- mg[, numeric_columns] / mg$TL

#get all the numeric columns again
numeric_columns <- sapply(mg, is.numeric)
#visually check for outliers
# Create scatter plots for each combination of numeric columns
for (i in which(numeric_columns)) {
  for (j in which(numeric_columns)) {
    if (i != j) {
      plot(mg[, i], mg[, j], 
           xlab = colnames(mg)[i], 
           ylab = colnames(mg)[j],
           #method = "overplot",
           main = paste("Scatter plot of:", colnames(mg)[i], "vs.", colnames(mg)[j]))
    }
  }
}

#change TL to cms
mg$TL <- mg$TL/10

#tables to check categorical variables, and tidy where required
table(mg$BodyShapeI) #looks good
table(mg$BodyShapeII) #looks good
table(mg$PosofMouth) #good
table(mg$mandible) #good
table(mg$eye.position) #good
table(mg$spiracle) #good
table(mg$caudal.fin.shape) #good

rm(list= ls()[! (ls() %in% c('mg'))])
#saveRDS(mg,"##/data/Devonian_fish_traits_tidy_RDS/Devonian_traits_Miguasha&Gogo_tidy.rds")