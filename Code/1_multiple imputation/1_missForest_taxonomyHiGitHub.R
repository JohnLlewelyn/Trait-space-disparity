#add supertree information and run missForest using taxonomy/supertree and traits
library(readxl)
library(dplyr)
library(tidyr)
library(missForest)

#Adjust file path at"###"

#get taxonomy
tax <- read_excel("###/data/taxonomy/taxonomy.xlsx")
tax$site <- NULL

#get traits
traits <- read_excel("###/Gogo&Miguasha_fish_traits.xlsx")

#remove pre-caudal length (too many missing in Devonian and modern data sets)
traits$`pre-caudal length if different from standard length`<-NULL

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
          "caudal fin shape (koaw.org)", "genus", "family", "order", "superorder/subclass", "class", "supertree1",
          "supertree2", "supertree3")

traits <- traits[,names(traits)%in%keep]

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

traits$`standard length` <-NULL

##now impute traits using missforest#######
set.seed(11)
imputed_data <- missForest(traits[,c(3:14,17:22)],mtry=17)
imputed_data$OOBerror

#optimise with number of taxonomic levels included, order columns so broadest taxonomy first
traitsT <- traits[,c(1:14,22:21,20:17)]

#hyper grid to test combinations
hyper_grid <- expand.grid(
  colMax = c(14:20),
  mtry_frac = seq(0.4, 0.9, 0.05),
  ntrees = c(100,200,300)
)
hyper_grid$mtry = round((hyper_grid$colMax-3)*hyper_grid$mtry_frac)

set.seed(15)
for(i in seq_len(nrow(hyper_grid))) {
  # fit model for ith hyperparameter combination
  imputed_data <- missForest(traitsT[,3:hyper_grid$colMax[i]],mtry = hyper_grid$mtry[i],ntree = hyper_grid$ntree[i])
  print(c(hyper_grid$colMax[i],hyper_grid$mtry[i],  hyper_grid$ntree[i], imputed_data$OOBerror))
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
set.seed(4)
bestM <- missForest(traitsT[,3:18],mtry = 10, ntree = 100)
bestM$OOBerror
#NRMSE 0.29 PFC 0.11
names_used <- names(traitsT[,3:18])

#assign data from the best model
data <- bestM$ximp
data$taxon <- traits$taxon

#add site information back in
site <- traits[,names(traits)%in%c("taxon","site")]
data <- merge(data,site,by="taxon",all.x=TRUE)

#tidy
rm(list= ls()[! (ls() %in% c('data'))])
saveRDS(data, "###/data/Devonian_fish_traits_RDS/Mig&GoImp_2024.RDS")
