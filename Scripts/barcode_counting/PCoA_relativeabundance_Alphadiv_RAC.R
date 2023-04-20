rbdat = read.delim("C:/Users/User/Desktop/thesis_stuff/R_working/fitness_data/html/SmeliPlant/gene_counts.tab")
# rbdat = read.delim("~/GitHub/rhizo.rb.tnseq/Data/fitness_data/html/SmeliPlant/gene_counts.tab")

fitdat = read.delim("C:/Users/User/Desktop/thesis_stuff/R_working/fitness_data/html/SmeliPlant/fit_logratios_clean.tab", row.names = 1)
# fitdat = read.delim("~/GitHub/rhizo.rb.tnseq/Data/fitness_data/html/SmeliPlant/fit_logratios.tab", row.names = 1)

cleanrbdat = t(read.delim("C:/Users/User/Desktop/thesis_stuff/R_working/fitness_data/html/SmeliPlant/gene_counts_clean.tab", row.names = 1))

goodcountdat <- subset(cleanrbdat, colnames(fitdat) %in% c("X21", "X22", "X24", "X26", "X27", "X28", "X29", "X30", "X32", "X33", "X35", "X37", "X40", "X41", "X42", "X43", "X44", "X45"))

rownames(rbdat) <- rbdat$comb
rbdat = t(rbdat[,-c(1:4)])
metadata = as.data.frame(matrix(ncol=2,nrow=21))
metadata[,1] <-c(rep("nod",9),rep("bact",7),rep("start",5))
metadata[,2] <-c(rep("fertilized",4),rep("ambient",5),rep("fertilized",4),rep("ambient",3),rep("time0",5))
colnames(metadata) <- c("compartment","Ntreatment")  #this cleans up and labels data by treatment/compartment

library(vegan)
library(ggplot2)
library(stringr)
library(tidyverse)
library(dplyr)
require("psych")

rel.dat = decostand(rbdat,method="total")

relnewnames <- substr(row.names(rel.dat), 7, 9)

relnewnames

rownames(rel.dat) <- relnewnames #rename rows to extract number ID

reldat2 <- subset(rel.dat, row.names(rel.dat) %in% c("021", "022", "023", "024", "026", "027", "028", "029", "030", "031", "032", "033", "035", "037", "039", "040"))

metadata2 <- subset(metadata, compartment!= "start") #reldat2 and metadata2 have the t0 data filtered out



#this part of the script generates alpha diversity data for the library mutant based on the data that passes the vertical line test

shannonfit <- numeric(nrow(goodcountdat))
for(j in 1:nrow(goodcountdat)) {
  shannonfit[j] <- vegan::diversity(goodcountdat[j,], index = "shannon")
}

shannonfit #shannon diversity of lib

Evar <- function(x) {
  x <- as.vector(x[x > 0])
  1 - (2/pi)*atan(var(log(x)))
} #this function calculates the Smith and Wilson's Evenness for the given site

Evarfit <- numeric(nrow(goodcountdat))
for(j in 1:nrow(goodcountdat)) {
  Evarfit[j] <- Evar(goodcountdat[j,])
}

Evarfit #evenness of lib in each pot

countzeros <- function(row) {
  numzeros <- sum(row == 0)
  percentzeros <- 100 - (numzeros / length(row) * 100)
  return(percentzeros)
} #returns the percentage of species found in each site


percentfoundfit <- apply(goodcountdat, 1, countzeros) #calculates % found for each row (site)

percentfoundfit

S_obs <- function(x = ""){
  rowSums(x>0)*1
} #observed species richness in a given site

richnessfit <-  S_obs(goodcountdat)

richnessfit

alphadiversitylib = as.data.frame(matrix(ncol=4,nrow=nrow(goodcountdat)))
colnames(alphadiversitylib) <- c('Shannon diversity','Evenness (SW)','% Found','Richness')
rownames(alphadiversitylib) <- rownames(goodcountdat)

alphadiversitylib[,1] <- shannonfit

alphadiversitylib[,2] <- Evarfit

alphadiversitylib[,3] <- percentfoundfit

alphadiversitylib[,4] <- richnessfit

alphadiversitylib #final collated dataframe of alpha diversity measures by pot

#this part of the script produces rank-abundance curves for each treatment/compartment combination

RAC <- function(x = ""){
  x = as.vector(x)
  x.ab = x[x > 0]
  x.ab.ranked = x.ab[order(x.ab, decreasing = TRUE)]
  return(x.ab.ranked)
} #returns a rank-abundance dataset for a given site to be plotted

RACmetadata <- metadata

row.names(RACmetadata) <- c("021", "022", "023", "024", "026", "027", "028", "029", "030", "031", "032", "033", "035", "037", "039", "040", "041", "042", "043", "044", "045")

RACmetadata <- subset(RACmetadata, row.names(RACmetadata) %in% c("021", "022", "024", "026", "027", "028", "029", "030", "032", "033", "035", "037", "040", "041", "042", "043", "044", "045"))

row.names(goodcountdat) <- c("021", "022", "024", "026", "027", "028", "029", "030", "032", "033", "035", "037", "040", "041", "042", "043", "044", "045")

RACsubsetter <- as.data.frame(cbind(goodcountdat, RACmetadata))

RACsubsetter <- RACsubsetter %>% mutate(group = paste(compartment, Ntreatment, sep = "-"))

RACsubsetter <- RACsubsetter[, -c(5176:5177)]

RACsubsetter <- RACsubsetter %>%
  mutate(across(1:5175, as.numeric))

RACdata <- RACsubsetter %>%
  group_by(group) %>%
  summarise(across(starts_with("SM"), mean))

RAC_bactamb <- RAC(RACdata$bact-ambient)
                                                                        
#this part of the script generates dissimilarity matrices for ranked gene fitness analyses

fitmetadata <- read.csv("C:/Users/User/Desktop/analyze_Nfitness/smeliplant_neutral_cog_clean.csv")
#fitmetadata <- read.csv("~/Github/rhizo.rb.tnseq/Scripts/reciprocal_BLAST_essentiality/rbh_essentiality_supplies/smeliplant_neutral_cog_clean.csv")

fitmetadata <- fitmetadata[, -c(8:11)] #clean up empty columns
fitmetadata <- fitmetadata[, -c(1)] #same here

newnamesfitrel <- substr(row.names(fitdat), 0, 9) #remove description info from labels
newnamesfitrel
row.names(fitdat) <- newnamesfitrel

truenames <- grep("^SMa", row.names(fitdat)) #remove trailing space and letter from SMa names
row.names(fitdat)[truenames] <- substr(row.names(fitdat)[truenames], 0, 7) 

truefitnames <- grep("^SMc", row.names(fitdat)) 
row.names(fitdat)[truefitnames] <- gsub("\\s+$", "", row.names(fitdat)[truefitnames])  #remove trailing space from SMc names

row.names(metadata) <- c("021", "022", "023", "024", "026", "027", "028", "029", "030", "031", "032", "033", "035", "037", "039", "040", "041", "042", "043", "044", "045")

metadata <- as.data.frame(t(metadata))

completefit <- merge(fitmetadata, fitdat, by.x = 1, by.y = 0)
nodESfit <- subset(completefit, essentiality=="essential")
nodNEfit <- subset(completefit, essentiality=="N") #fitness data subsetted into groups of interest

nitrogengenes <- read.delim("C:/Users/User/Desktop/blastkoala_annotation_smeli.txt", header=FALSE)
nitrogengenes <- read.delim("~/Github/rhizo.rb.tnseq/Data/fitness_data/html/SmeliPlant/blastkoala_annotation_smeli.txt", header=FALSE)

nitrogengenes$V1 <- substr(nitrogengenes$V1, 7, 15)
Nnames <- grep("^SMa", nitrogengenes$V1)
nitrogengenes$V1[Nnames] <- substr((nitrogengenes$V1)[Nnames], 0, 7)
nitronames <- grep("^SMc", nitrogengenes$V1) 
nitrogengenes$V1[nitronames] <- gsub("\\s+$", "", nitrogengenes$V1[nitronames])

metabloci <- c("SMa0228", "SMc04028", "SMc04026", "SMa1250", "SM_b20986", "SMa1182", "SMc02150", "SMa0697", "SMc04083", "SMa0045", "SMc02613", "SMc00762", "SMc00948", "SMc01594", "SMc01973", "SMc02352", "SM_b20745", "SMa1276", "SMa1236", "SMa1233", "SM_b20436", "SMa0827", "SMa0825", "SMa0829", "SMa1273", "SMc04085", "SMa0585", "SMa0583", "SMa0581", "SM_b20985", "SM_b20984") #from Nmetabolismlist.txt, which is sourced from blastKOALA

Nmetabolism <- subset(nitrogengenes, V1 %in% metabloci)

Nmetabolismfit <- subset(completefit, completefit$locus_tag %in% Nmetabolism$V1) #N metabolism fitness subset

fitmetadata <- fitmetadata[, -c(8:11)]

nodESfit <- nodESfit[, -c(2:6)]
nodNEfit <- nodNEfit[, -c(2:6)]
completefit <- completefit[, -c(2:6)]
Nmetabolismfit <- Nmetabolismfit[, -c(2:6)]

row.names(nodESfit) <- nodESfit[,1]
nodESfit <- nodESfit[, -c(1)]

row.names(nodNEfit) <- nodNEfit[,1]
nodNEfit <- nodNEfit[, -c(1)]

row.names(completefit) <- completefit[,1]
completefit <- completefit[, -c(1)]

row.names(Nmetabolismfit) <- Nmetabolismfit[,1]
Nmetabolismfit <- Nmetabolismfit[, -c(1)]

dissim_nodES <- vegdist(nodESfit, method = "euclid", diag = FALSE)
dissim_nodNE <- vegdist(nodNEfit, method = "euclid", diag = FALSE)
dissim_complete <- vegdist(fitdat, method = "euclid", diag = FALSE)

scaled_nodES <- scale(dissim_nodES)
scaled_nodNE <- scale(dissim_nodNE)
scaled_complete <- scale(dissim_complete)

cor_nodES <- cor(scaled_nodES)
cor_nodNE <- (scaled_nodNE)
cor_complete <- cor(scaled_complete)

KMO(cor_nodES)
KMO(cor_nodNE)
KMO(cor_complete)

cortest.bartlett(cor_nodES, n = nrow(cor_nodES))
cortest.bartlett(cor_nodNE, n = nrow(cor_nodNE))
cortest.bartlett(cor_complete, n = nrow(cor_complete))

pca_nodES <- rda(scaled_nodES, scale = TRUE)
pca_nodNE <- rda(scaled_nodNE, scale = TRUE)
pca_complete <- rda(scaled_complete, scale = TRUE) #PCA analysis 



#reorganizing script; here is the pcoa and rarefaction stuff, and the top 50/bottom 50 list

pcoa = capscale(reldat2~1,distance = "bray")
head(summary(pcoa))
biplot(pcoa)
pcoa.sites = scores(pcoa,display = "sites") #see loadings, takes awhile

pcoa.sites = as.data.frame(scores(pcoa,display = "sites"))

ggplot(pcoa.sites,aes(x=MDS1,y=MDS2,color=metadata2$compartment,shape=metadata2$Ntreatment)) +
  geom_point(size=8)+
  xlab("PCo 1 (18%)") +
  ylab("PCo 2 (14%)") +
  theme_bw(base_size = 30) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 1, colour = "black",fill=NA),
        panel.background = element_rect(fill="white"),
        axis.text.x=element_text(size=rel(0.9),angle=45,h=1),
        axis.text.y=element_text(size=rel(0.9))) #initial pcoa assessment to see where the groups fall

# rarefaction curve

counting <- rowSums(rbdat) #total reads per site/pot

counting

raremax <- min(rowSums(rbdat))

raremax #smallest number of reads for a site

Srare <- rarefy(rbdat, raremax)

Srare #rare genes

labelnames <- data.frame(rownames(rbdat)) #make smaller labels for rarefaction curve 

shortnames <- str_sub(labelnames$rownames.rbdat., 7, 22)

shortnames

cleannames <- gsub('\\.', ' ',shortnames)

numbersonly <- substring(labelnames$rownames.rbdat., 7, 9)

cleannames #these are meant to be passed through ordilabel() to create cleaner labels for the rarefaction curve but I'm still working on the syntax

rownames(rbdat) <- numbersonly

cols <- c("darkred", "darkred", "darkred", "darkred", "forestgreen", "forestgreen", "forestgreen", "forestgreen", "forestgreen", "blue", "blue", "blue", "blue", "purple", "purple", "purple", "orange", "orange","orange","orange","orange")

rbrar <- rarecurve(rbdat, step = 500, sample = 1000000, xlab = "Total barcodes", ylab = "Unique barcodes", col = cols, cex = 0.5, label = TRUE)

rbrar

# the following is a ChatGPT-generated script to extract the 50 most abundant genes per experiment... wow!! 

#top_bc <- apply(rbdat, 1, function(x) {
#  sorted <- sort(x, decreasing = TRUE)
#  names(sorted)[1:50]
#})

# Order the top species by their abundance
#ordered_top_bc <- apply(rbdat, 1, function(x) {
#  sorted <- sort(x, decreasing = TRUE)
#  names(sorted)[1:50][order(sorted[1:50], decreasing = TRUE)]
#})

# Print the ordered top species for each row (site)
#print(ordered_top_bc)

#end of the generated portion; however, note that I worked with + modified this script in order to produce the following code

# Extract the top 50 species for each column (species)
top_fit <- apply(fitdat, 2, function(x) {
  sorted <- sort(x, decreasing = TRUE)
  names(sorted)[1:50]
})

# Order the top species by their abundance
ordered_top_fit <- apply(fitdat, 2, function(x) {
  sorted <- sort(x, decreasing = TRUE)
  names(sorted)[1:50][order(sorted[1:50], decreasing = TRUE)]
})

orderedtopfit <- ordered_top_fit[,-1:-2]

# Print the ordered top species for each column (species)
print(orderedtopfit)

bottom50genes <- t(orderedtopfit) #the bc insertionally inactivates the gene, so more bc in the t1 sample = beneficial inactivation = lower fitness

lowest50fit <- cbind(metadata,bottom50genes) #generates list of 50 genes with lowest fitness values per treatment, with metadata

#the following code will produce a list of the 50 genes with the highest fitness

# Extract the bottom 50 species for each column (species)
bottom_fit <- apply(fitdat, 2, function(x) {
  sorted <- sort(x)
  names(sorted)[1:50]
})

# Order the bottom species by their abundance
ordered_bottom_fit <- apply(fitdat, 2, function(x) {
  sorted <- sort(x)
  names(sorted)[1:50][order(sorted[1:50])]
})

orderedbottomfit <- ordered_bottom_fit[,-1:-2]

# Print the ordered bottom species for each column (species)
print(orderedbottomfit)

top50genes <- t(orderedbottomfit)

top50fit <- cbind(metadata,top50genes) #like lowest50 but with the highest fitness

# write.csv(lowest50fit, "C:/Users/User/Desktop/lowest50fitgenes.csv", row.names = TRUE)

# write.csv(top50fit, "C:/Users/User/Desktop/top50fitgenes.csv", row.names = TRUE)