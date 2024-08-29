#combine RDSs (after adding names) and check/tidy data

library(ggplot2)
library(dplyr)
library(rfishbase)

#Adjust file path at"###"

# Set working directory to  folder containing the RDS files for modern fish
setwd("###/data/modern_fish_traits_RDSs")

# List all the .rds files in the folder
rds_files <- list.files(pattern = "\\.[rR][dD][sS]$")

# Initialize an empty list to store the data
data_list <- list()

# Loop through each .rds file, read it, and add a column with the file name
for (file in rds_files) {
  # Read the .rds file
  data <- readRDS(file)
  
  # Extract the file name (without extension)
  file_name <- tools::file_path_sans_ext(file)
  
  # Add a column with the file name
  data$file_name <- file_name
  
  # Add the data to the list
  data_list[[file_name]] <- data
}

#use function to keep (and order) the columns needed
trimmed <- lapply(data_list,function(x) {x <- x[,c("Species","BodyShapeI","BodyShapeII","TL","SL","HL","ED",
                       "POL","BD","PosofMouth","mandible","eye.position","spiracle",
                       "caudal.fin.shape","file_name")]})

#combine into one df
traits <- as.data.frame(do.call("rbind",trimmed))
dim(traits) 
length(unique(traits$Species)) #some fish are found in more than one assemblage

#check they have the same values (unlikely to have exactly the same morphometrics if measured manually)
dups <- traits$Species[duplicated(traits$Species)]
dups <- traits[traits$Species%in%dups,] 

#most have identical records, filter to see which don't
#need to round the morphometrics first
traits$SL <- round(traits$SL,6)
traits$HL <- round(traits$HL,6)
traits$ED <- round(traits$ED,6)
traits$POL <- round(traits$POL,6)
traits$BD <- round(traits$BD,6)
traits$TL <- round(traits$TL,6)

dups_cut <- dups[,names(dups)!="file_name"]
dups_cut <- distinct(dups_cut) #got rid of a bunch of duplicates, keeping only one row for these species. Need to cut back to only those that are duplicated now.
dups_cut <- dups_cut[duplicated(dups_cut$Specie) | duplicated(dups_cut$Species, fromLast = TRUE), ]

#fix these so there's only one set of traits for each species - make changes directly in "traits"
traits$PosofMouth[traits$Species=="Eucinostomus argenteus"] <- "terminal"
traits$caudal.fin.shape[traits$Species=="Hemiramphus brasiliensis"] <- "homocercal"
traits$BodyShapeII[traits$Species=="Hyporhamphus unifasciatus"] <- "compressed"
traits$POL[traits$Species=="Lutjanus jocu"] <- 14.19868
traits$TL[traits$Species=="Lutjanus synagris"] <- 60
traits$TL[traits$Species=="Mycteroperca bonaci"] <- 150
traits$PosofMouth[traits$Species=="Oligoplites saurus"] <- "superior"
traits$TL[traits$Species=="Synodus foetens"] <- 53.8

#now check again
#check they have the same values (unlikely to have exactly the same morphometrics if I measured them manually)
dups <- traits$Species[duplicated(traits$Species)]
dups <- traits[traits$Species%in%dups,]
dups_cut <- dups[,names(dups)!="file_name"]
dups_cut <- distinct(dups_cut) #got rid of a bunch of duplicates, keeping only one row for these species. Need to cut back to only those that are duplicated now
dups_cut <- dups_cut[duplicated(dups_cut$Specie) | duplicated(dups_cut$Species, fromLast = TRUE), ] #all good now

#which fish had TL added manually (not from fishbase) 
fsp <- species(traits$Species, fields = c("Species", "Length","LTypeMaxM"))
#3 rays that had dsic width rather than length: Aetobatus narinari, Hypanus americanus, Hypanus guttatus  
#1 fish NG - Not Given: Etropus longimanus. 
#2 fish NAs: Acanthurus bahianus, Hemisorubim platyrhynchos
#the TLs for these species have been double checked and references added where needed.

#visually check for outliers in each column - absolute measures and relative lengths
#numeric_columns <- traits[, sapply(traits, is.numeric)]
#pairs(numeric_columns)
plot(traits$TL,traits$TL)

#two TLs are way out
TL <- head(traits[order(traits$TL, decreasing = TRUE), ], n = 2) #tiger shark is 750 on fishbase, but Aetobatus narinari seems too long
#Aetobatus narinari is a  large stingray - some sources say it can reach 8.8 meters but this might be an over estimation. 
#Using information from Florida Museum, we now assume a maximum total length of 5 meters and adjust the other measures.
traits$SL = ifelse(traits$Species == "Aetobatus narinari"& traits$TL==880, traits$SL / 880 * 500, traits$SL)
traits$HL = ifelse(traits$Species == "Aetobatus narinari"& traits$TL==880, traits$HL / 880 * 500, traits$HL)
traits$ED = ifelse(traits$Species == "Aetobatus narinari"& traits$TL==880, traits$ED / 880 * 500, traits$ED)
traits$POL = ifelse(traits$Species == "Aetobatus narinari"& traits$TL==880, traits$POL / 880 * 500, traits$POL)
traits$BD = ifelse(traits$Species == "Aetobatus narinari"& traits$TL==880, traits$BD / 880 * 500, traits$BD)
traits$TL = ifelse(traits$Species == "Aetobatus narinari"& traits$TL==880, traits$TL / 880 * 500, traits$TL)

#change morphometrics to proportion of TL 
# Identify the numeric columns (excluding "TL")
numeric_columns <- sapply(traits, is.numeric)
numeric_columns <- numeric_columns & !colnames(traits) == "TL"
# Divide the numeric columns by "TL"
traits[, numeric_columns] <- traits[, numeric_columns] / traits$TL

#visually check for outliers
# Create scatter plots for each combination of numeric columns
for (i in which(numeric_columns)) {
  for (j in which(numeric_columns)) {
    if (i != j) {
      plot(traits[, i], traits[, j], 
           xlab = colnames(traits)[i], 
           ylab = colnames(traits)[j],
           main = paste("Scatter plot of", colnames(traits)[i], "vs.", colnames(traits)[j]))
    }
  }
}


#three SLs < 0.7 
SL <- tail(traits[order(traits$SL, decreasing = TRUE), ], n = 3) #Gymnogeophagus balzanii checks out; Gobionellus oceanicus SL from fishbase looks a bit short, should ~0.73; Ctenogobius smaragdus checks out -> big tail
traits$SL[traits$Species=="Gobionellus oceanicus"] <- 0.73
#one ED > 0.13
ED <- head(traits[order(traits$ED, decreasing = TRUE), ], n = 1) #seems a bit large, but measurement comes from fishabse. Check two more images (0.13 and 0.1), so take average
traits$ED[traits$Species=="Myripristis jacobus"] <- 0.12
#two BDs > 0.6
BD <- head(traits[order(traits$BD, decreasing = TRUE), ], n = 2) #Chaetodon ocellatus from fishbase, looks about right, same for Pomacanthus paru
#two head lengths < 0.1
HL <- tail(traits[order(traits$HL, decreasing = TRUE), ], n = 2) #checks out, Myrichthys have small heads and long bodies
#five POL
POL <- traits[traits$POL>0.21,] #all check out

#fix some tail traits
traits$caudal.fin.shape[traits$Species=="Eigenmannia trilineata"] <- "reduced or absent"
traits$caudal.fin.shape[traits$Species=="Rhamphichthys hahni"] <- "reduced or absent"
traits$caudal.fin.shape[traits$Species=="Aetobatus narinari"] <- "whip-like"
traits$caudal.fin.shape[traits$Species=="Hypanus americanus"] <- "whip-like"

#check for missing data
# Check for NAs in each column
has_nas <- sapply(traits, function(x) any(is.na(x))) #no NAs
# Check for empty cells in each column
has_empty_cells <- sapply(traits, function(x) any(x == "")) #position of mouth
table(traits$PosofMouth) #one appears to be empty -> Melichthys niger
traits$PosofMouth[traits$Species=="Melichthys niger"] <- "terminal"

#tables to check categorical variables, and tidy where required
table(traits$BodyShapeI) #some need fixing, 5 classed as "other" or "other (see remarks)", make new category for flatfish in bodyshapeII
traits$BodyShapeI[traits$BodyShapeI=="Elongated"] <- "elongated"
traits$BodyShapeI[traits$BodyShapeI=="other (see remarks)"] <- "other"
flatfish <- c("Etropus longimanus","Paralichthys brasiliensis","Citharichthys spilopterus","Platichthys flesus",
                "Paralichthys microps","Pleuronectes platessa","Etropus crossotus","Paralichthys adspersus",
                "Bothus lunatus","Achirus declivis","Bothus ocellatus","Achirus lineatus")
traits$BodyShapeII[traits$Species%in%flatfish] <- "compressed and lies on side"
traits$BodyShapeI[traits$Species=="Etropus longimanus"] <- "fusiform / normal"
traits$BodyShapeI[traits$Species=="Etropus crossotus"] <- "short and / or deep"
#fix bat fish, strange body shape
traits$BodyShapeI[traits$Species=="Ogcocephalus nasutus"] <- "other"
traits$BodyShapeII[traits$Species=="Ogcocephalus nasutus"] <- "angular"

table(traits$BodyShapeII) #looks good
table(traits$PosofMouth) #good
table(traits$mandible) #good
table(traits$eye.position) #good
table(traits$spiracle) #good
table(traits$caudal.fin.shape) #need to fix some caudal fin shapes
traits$caudal.fin.shape[traits$Species=="Rhamphichthys hahni"] <- "reduced or absent"
traits$caudal.fin.shape[traits$Species=="Aetobatus narinari"] <- "whip-like"
traits$caudal.fin.shape[traits$Species=="Hypanus americanus"] <- "whip-like"

#fix column name and add some cols
names(traits)[names(traits)=="file_name"] <- "community"
traits$habitat <- ifelse(traits$community%in%c("traits_Caribbean","traits_Chile_reef"), "reef",
                                               ifelse(traits$community%in%c("Santa_Cruz_Channel","traits_Ythan"), "estuary", "freshwater"))
degs <- data.frame(community = c("traits_Caribbean","traits_Chile_reef","Santa_Cruz_Channel","traits_Ythan","traits_BracoMorto","traits_Little_Rock_Lake","traits_Nepean"), degsFromEquator = c(18,31.5,7,57,16,46,34))
traits <- merge(traits, degs,by="community")

#saveRDS(traits,"###/data/modern_fish_traits_RDSs/combined/together&tidied.rds")
