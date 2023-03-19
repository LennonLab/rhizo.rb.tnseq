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

gene_list <- scores(pcoa,display="species")
