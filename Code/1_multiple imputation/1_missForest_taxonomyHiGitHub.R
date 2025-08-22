#add supertree information and run missForest using taxonomy/supertree and traits to impute missing data
library(readxl)
library(dplyr)
library(tidyr)
library(missForest) 

#Adjust file path at"filepath/"
filepath <- "filepath/"
setwd(filepath)

#get taxonomy
tax <- read_excel(paste(filepath,"/data/taxonomy/taxonomy.xlsx",sep=""))
tax$site <- NULL
#get traits
traits <- read_excel(paste(filepath,"/data/Devonian fish traits/Devonian_traits_and_taxonomy-with Canowindra_v1.xlsx",sep=""))

#remove traits not needed for analyses
traits <- traits[, setdiff(names(traits),c("similar-shaped extant species","specimen number","lit/Reference","different max. total length", 
                                           "pre-caudal length if different from standard length", "body shape","lungs", 
                                           "caudal_fin_shape_known_inferred", "caudal height", "unpaired fins", "mouth gape length", 
                                           "stomach contents", "NOTES", "genus", "family", "order", "superorder/subclass"))]
#add higher taxonomic groupings
tax$supertree1 <- ifelse(tax$class=="Sarcopterygii"|tax$class=="Actinopterygii", "Sarcopt-Actinopt",
                         ifelse(tax$class=="Acanthodii"|tax$class=="Chondrichthyes", "Acanth-Chondri",tax$class))

tax$supertree2 <- ifelse(tax$class=="Placodermi"|tax$supertree1=="Sarcopt-Actinopt"|tax$supertree1=="Acanth-Chondri", 
                         "pelvic_fin_grp", tax$class)

tax$supertree3 <- ifelse(tax$class=="Osteostraci"|tax$supertree2=="pelvic_fin_grp", 
                         "osteo_and_higher", "basal_fish")

#combine with trait data
traits <- merge(traits, tax, by = "taxon", all.x = TRUE)

#cut to traits needed for trait space
keep <- c("taxon","site", "total length", "standard length", "head length", "pre-orbital length (snout length)",
          "body depth (body height)", "eye size (diameter)", "body shape I (saggital plane)",
          "body shape II (transverse plane)", "mandible", "mouth position","eye position", "spiracular",
          "caudal fin shape (koaw.org)", "group", "supertree3", "supertree2", "supertree1", "class","superorder/subclass","order","family")

traits <- traits[, keep[keep %in% names(traits)]]

#number of NAs by species
missing <- data.frame(taxon = traits$taxon, NAs = as.numeric(rowSums(is.na(traits))))

#fix BodyShapeI - fusiform had been split into two groups
traits$`body shape I (saggital plane)` <- gsub("[^a-zA-Z0-9 / -]", "", traits$`body shape I (saggital plane)`)
traits$`body shape I (saggital plane)`[traits$`body shape I (saggital plane)`=="fusiform / normal"] <-  "fusiform/normal"
traits$`body shape I (saggital plane)`[traits$`body shape I (saggital plane)`=="short and / or deep"] <-  "short and/or deep"
#fix mandible data
traits$mandible[traits$taxon=="Cainocara enigma"] <- "present"
traits$mandible[traits$taxon=="Kapitany sarcopt"] <- "present"

#make sure numeric columns are numeric, replace "-" with NA, and remove c. from numeric columns
nums <- c("total length", "standard length", "head length", "pre-orbital length (snout length)", 
          "body depth (body height)", "eye size (diameter)")
traits[] <- lapply(traits, function(x) ifelse(x == "-", NA, x))
traits[] <- lapply(traits, function(x) gsub("c\\.", "", x))
traits[nums] <- lapply(traits[nums], as.numeric)
traits[!names(traits)%in%nums] <- lapply(traits[!names(traits)%in%nums], as.factor)

#######Calculate proportion of NAs for each community####################
#identify the column range
start_col <- which(names(traits) == "total length")
end_col <- which(names(traits) == "caudal fin shape (koaw.org)")

#get the relevant column names (excluding "standard length")
trait_cols <- names(traits)[start_col:end_col]
trait_cols <- trait_cols[trait_cols != "standard length"]

#calculate NA proportion per community
na_proportion <- traits %>%
  group_by(site) %>%
  summarise(across(all_of(trait_cols), ~mean(is.na(.)), .names = "na_prop_{.col}")) %>%
  ungroup()

rowMeans(data.frame(na_proportion[,2:13]))

#########################################################################

#how many are missing per column
na_perc <- function(x) {
  mis <- table(is.na(x))[2]/length(x)*100
  mis <- ifelse(is.na(mis), 0, mis)
  return(mis)}
mis <- data.frame(trait=names(traits),missingPerc=unlist(lapply(traits,na_perc)))

#how many are missing per family
missing_family <- function(df, grp) {
  df %>% 
    gather(key = "trait", value = "value", -!!grp) %>% 
    group_by(!!grp) %>% 
    summarise(proportion_NA = mean(is.na(value)))
}

#how many are missing per species
missSp <-  data.frame(species=as.character(traits$taxon),prop_miss=rowMeans(is.na(traits[,])))

#list of species with few traits/taxonomic info
remSp <- missSp$species[missSp$prop_miss>0.52]

#remove standard length
#traits$`standard length` <-NULL

# Save traits for Devonian traits supp mat
write.csv(traits, "filepath/Table S1 Devonian species traits.csv")

library(flextable)
library(officer)
library(dplyr)

# Create a Word document
doc <- read_docx() %>%
  body_add_flextable(flextable(traits)) %>%
  body_add_par("", style = "Normal")

# Save as Word
print(doc, target = "filepath/Table S1 Devonian species traits.docx")

##now impute traits using missforest#######
set.seed(11)
imputed_data <- missForest(traits[,c(3:22)])
imputed_data$OOBerror
# NRMSE 0.30
# PFC 0.11

#optimise with number of taxonomic levels included, order columns so broadest taxonomy first
#traitsT <- traits[,c(1:14,22:21,20:17)]

#hyper grid to test combinations
hyper_grid <- expand.grid(
  colMax = c(16:22),
  mtry_frac = seq(0.4, 0.9, 0.05),
  maxiter = c(10, 20),
  ntrees = c(300)  #could just do 300 trees to be quicker
)
hyper_grid$mtry = round((hyper_grid$colMax-3)*hyper_grid$mtry_frac)

set.seed(15)
for(i in seq_len(nrow(hyper_grid))) {
  # fit model for ith hyperparameter combination
  imputed_data <- missForest(traits[,3:hyper_grid$colMax[i]],
                             mtry = hyper_grid$mtry[i],
                             ntree = hyper_grid$ntree[i],
                             maxiter = hyper_grid$maxiter[i])
  print(i)
  print(c(hyper_grid$colMax[i],hyper_grid$mtry[i],  hyper_grid$ntree[i],hyper_grid$maxiter[i], imputed_data$OOBerror))
  # export OOB error
  hyper_grid$NRMSE[i] <- imputed_data$OOBerror[1]
  hyper_grid$PFC[i] <- imputed_data$OOBerror[2]
}
#scale and combine performance measures
hyper_grid$NRMSE_scale <- scale(hyper_grid$NRMSE)
hyper_grid$PFC_scale <- scale(hyper_grid$PFC)
hyper_grid$comb_perf <- (hyper_grid$NRMSE_scale+hyper_grid$PFC_scale)/2
hyper_grid[which(hyper_grid$comb_perf==(min(hyper_grid$comb_perf))),c(1,3,4)]

#best performance overall
set.seed(2)
bestM1 <- missForest(traits[,3:20],mtry = 15, ntree = 300, maxiter = 20)
set.seed(2)
bestM <- missForest(traits[,3:20],mtry = 15, ntree = 300, maxiter = 20, variablewise = TRUE)
# 
bestM1$OOBerror
# NRMSE 0.14 PFC 0.09 # perform well, but some of this is due to performance on taxonomy columns

# Calculate normalized performance for each column
# Columns of interest: 3 to 15
cols_to_check <- 3:15

# Extract the MSEs from missForest output
mse_vec <- bestM$OOBerror[1:length(cols_to_check)]
rmse_vec <- sqrt(mse_vec)

# Initialise result storage
nrmse_mean_vec <- numeric(length(cols_to_check))
names(nrmse_mean_vec) <- colnames(traits)[cols_to_check]

# Loop through columns and calculate NRMSE (mean)
for (i in seq_along(cols_to_check)) {
  col_data <- traits[[cols_to_check[i]]]
  mean_val <- mean(col_data, na.rm = TRUE)
  nrmse_mean_vec[i] <- rmse_vec[i] / mean_val
}

# View results
round(nrmse_mean_vec, 4)
tb <- data.frame(cbind(names(traits[,3:20]),bestM$OOBerror))
tb <- tb[grepl("PFC",row.names(tb)),]
tb$X2 <- as.numeric(tb$X2)
tb <- tb[tb$X2>0,]
mean(tb$X2)

# mean performance:NRMSE 0.29; PFC: 0.17

# extract names of columns used
names_used <- names(traits[,3:20])

# names of columns where perf is high (< 0.25 for NRMSE, < 0.3 PCF)
# Remove NA values before filtering and extracting names
filtered <- round(nrmse_mean_vec, 4)
filtered <- filtered[!is.na(filtered) & filtered < 0.3]
names_cut <- names(filtered)
cat_nms <- cbind(names(traits[,3:20]),bestM$OOBerror)[,1][cbind(names(traits[,3:20]),bestM$OOBerror)[,2]<0.25][1:6]
names_cut <- c(names_cut, cat_nms)


#assign data from the best model
data <- bestM$ximp
data$taxon <- traits$taxon

#add site information back in
site <- traits[,names(traits)%in%c("taxon","site")]
data <- merge(data,site,by="taxon",all.x=TRUE)

# remove higher taxonomy columns
data <- data[, !names(data)%in%c("group", "supertree3", "supertree2" , "supertree1", "class","superorder/subclass" )]

# cut to more accurately imputed traits
data_cut <- data[,names(data)%in%names_cut]

#tidy
rm(list= ls()[! (ls() %in% c('data',"data_cut","filepath"))])
saveRDS(data, paste(filepath,"/output/Mig&GoImp_2025.RDS", sep=""))
saveRDS(data_cut, paste(filepath,"/output/DevCut_2025.RDS", sep=""))
