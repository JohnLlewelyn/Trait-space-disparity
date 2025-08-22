# Chao1 for Canowindra
library(rarestR)

# # Species counts
# Placoderm, Antiarch:
# Bothriolepis yeungae (Johanson 1998) = >1500 individuals 
# Remigolepis walkeri (Johanson 1997) = >1500 individuals
# 
# Placoderm, Arthrodire:
# Groenlandaspis sp. nov. = 70 individuals 
# 
# Sarcopterygii, Dipnoi:
# Soederberghia simpsoni (Ahlberg et al. 2001) = 2 individuals 
# 
# Sarcopterygii, Rhizodont:
# Gooloogongia loomesi (Johanson & Ahlberg 1998) = 8 individuals
# 
# Sarcopterygii, Canowindrid:
# Canowindra grossi (Thomson 1973) = 1 individual
# 
# Sarcopterygii, Tristichopterid:
# Mandageria fairfaxi (Johanson & Ahlberg 1997) = 14 individuals
# Cabonnichthys burnsi (Ahlberg & Johanson 1997) = 10 individuals
can <- data.frame(both=1500, remi = 1500, groen = 70, soed = 2, gool = 8, can = 1, man = 14, cab = 10)

# how many specimens in total
sp <- rowSums(can)

# Chao1
counts <- as.numeric(can[1, ])
f1 <- sum(counts == 1)
f2 <- sum(counts == 2)
s_obs <- sum(counts > 0)
s_obs + ifelse(f2 > 0, (f1^2) / (2 * f2), NA) #8.5 i.e., 9 species

library(iNEXT)
abund_vector <- as.numeric(can[1, ])
out <- iNEXT(abund_vector, q = 0, datatype = "abundance", endpoint = 50000)

# Then extract the row for m = 50000
out$iNextEst$size_based[40,]
out$iNextEst$coverage_based[40,] # 8.5 species


