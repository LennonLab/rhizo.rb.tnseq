rbdat = read.delim("C:/Users/User/Desktop/thesis_stuff/R_working/fitness_data/html/SmeliPlant/gene_counts.tab")
# rbdat = read.delim("~/GitHub/rhizo.rb.tnseq/Data/fitness_data/html/SmeliPlant/gene_counts.tab")

rownames(rbdat) <- rbdat$comb
rbdat = t(rbdat[,-c(1:4)])
metadata = as.data.frame(matrix(ncol=2,nrow=21))
metadata[,1] <-c(rep("nod",9),rep("bact",7),rep("start",5))
metadata[,2] <-c(rep("fertilized",4),rep("ambient",5),rep("fertilized",4),rep("ambient",3),rep("time0",5))
colnames(metadata) <- c("compartment","Ntreatment")

library(vegan)
library(ggplot2)
library(stringr)
library(tidyverse)
library(dplyr)

rel.dat = decostand(rbdat,method="total")

pcoa = capscale(rel.dat~1,distance = "bray")
head(summary(pcoa))
biplot(pcoa)
pcoa.sites = scores(pcoa,display = "sites") #see loadings, takes awhile

pcoa.sites = as.data.frame(scores(pcoa,display = "sites"))

ggplot(pcoa.sites,aes(x=MDS1,y=MDS2,color=metadata$compartment,shape=metadata$Ntreatment)) +
  geom_point(size=8)+
  xlab("PCo 1 (25%)") +
  ylab("PCo 2 (14%)") +
  theme_bw(base_size = 30) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 1, colour = "black",fill=NA),
        panel.background = element_rect(fill="white"),
        axis.text.x=element_text(size=rel(0.9),angle=45,h=1),
        axis.text.y=element_text(size=rel(0.9)))

# rarefaction curve

counting <- rowSums(rbdat)

counting

raremax <- min(rowSums(rbdat))

raremax

Srare <- rarefy(rbdat, raremax)

Srare

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


fitdat = t(read.delim("C:/Users/User/Desktop/thesis_stuff/R_working/fitness_data/html/SmeliPlant/fit_logratios_clean.tab", row.names = 1))
# fitdat = read.delim("~/GitHub/rhizo.rb.tnseq/Data/fitness_data/html/SmeliPlant/fit_logratios.tab", row.names = 1)

cleanrbdat = t(read.delim("C:/Users/User/Desktop/thesis_stuff/R_working/fitness_data/html/SmeliPlant/gene_counts_clean.tab", row.names = 1))

goodcountdat <- subset(cleanrbdat, row.names(fitdat) %in% c("X21", "X22", "X24", "X26", "X27", "X28", "X29", "X30", "X32", "X33", "X35", "X37", "X40", "X41", "X42", "X43", "X44", "X45"))

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

#calculates % found for each row (site)
percentfoundfit <- apply(goodcountdat, 1, countzeros)

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