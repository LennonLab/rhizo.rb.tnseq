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
  geom_point()+
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

shortnames <- str_sub(labelnames$rownames.rbdat., 7, 20)

shortnames

cleannames <- gsub('\\.', ' ',shortnames)

cleannames #these are meant to be passed through ordilabel() to create cleaner labels for the rarefaction curve but I'm still working on the syntax

rbrar <- rarecurve(rbdat, step = 500, sample = raremax, xlab = "Total barcodes", ylab = "Unique barcodes", col = "blue", cex = 0.5, label = FALSE)

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


fitdat = read.delim("C:/Users/User/Desktop/thesis_stuff/R_working/fitness_data/html/SmeliPlant/fit_logratios.tab", row.names = 1)
# fitdat = read.delim("~/GitHub/rhizo.rb.tnseq/Data/fitness_data/html/SmeliPlant/fit_logratios.tab", row.names = 1)
                    
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

# would like to find a way to make a table of mean gene counts per treatment, then top/bottom mean gene fitness lists per treatment
# tried summarise and group_by with a cbind(metadata,rbdat) dataframe, but haven't figured out how to automate over all 5200 columns; might be able to rotate table
# once we have mean count data, should be really easy to generate meaned versions of top50fit, lowest50fit if desirable
