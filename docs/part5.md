# Multivariate Analyses of Microbial Communities with R
## Importing multivariate data using phyloseq
Loading the required packages
We recommend checking out some of the following references:  
[GUSTA ME](https://mb3is.megx.net/gustame)  
[Phyloseq Homepage](https://joey711.github.io/phyloseq/)  
[Ecological Analysis of Ecological Communities](http://cc.oulu.fi/~jarioksa/opetus/metodi/vegantutor.pdf)  

First we'll clear our R environment of all attached objects and define the memory limit for windows systems.
```python
# Memory limit only needed if you use windows system
#clear all currently attached R objects
rm(list=ls())
memory.limit()
memory.limit(size=56000)
```

Install the following packages. You will only need to do this once on your computer.
```python
source("https://bioconductor.org/biocLite.R")
biocLite('phyloseq')
install.packages(vegan)
biocLite("Biostrings")
install.packages("biostrings")
install.packages("reshape2")
install.packages("ape")
install.packages("picante")
install.packages("ggpubr")
install.packages("gridExtra")
install.packages("ggplot2")
install.packages("devtools")
install_github("vqv/ggbiplot")
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
```

Then load the packages we need:
```python
library(phyloseq)
library(vegan)
ibrary(Biostrings) # To read fasta file
library(reshape2)
library(ape) # to read tree file
library(scales)
library(picante)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(ggbiplot)
library(pairwiseAdonis)

#For reproducibility
set.seed(123)
```

### Importing our data 
Check out this [link](https://joey711.github.io/phyloseq/import-data.html) for the authors guide on importing data into phyloseq.

```python
#First set our working directory (this will be where you placed the data files
setwd("C:/IF/microbial community course")

###### Read data tables
otu_tb <- read.table("data/otu_table.txt", sep="\t", row.names=1,header=T)
tax_tb <- read.table("data/tax_table.txt", sep="\t", row.names=1,header=T)
metadata_tb <- read.csv("data/simulated_dataset.txt", sep="\t", dec = ".", row.names=1, header=TRUE)
tree <- read.tree("data/rep_otu.phylip.tre") 
refseq <- readDNAStringSet("data/rep_otu.fasta")

# import into phyloseq object
OTU <- otu_table(as.matrix(otu_tb), taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(tax_tb)) 
METADATA <- sample_data(metadata_tb)
TREE <- phy_tree(tree)

student_data <- phyloseq(OTU, TAX, METADATA, TREE, refseq)
save(student_data, file ="student_data.phyloseq")

# Check the phyloseq file
student_data

sample_names(student_data)
head(otu_table(student_data))
head(tax_table(student_data))
str(sample_data(student_data)) # check internal structure of the metadata file
```

## Generating intial community composition figures to explore our data
First off we'll create a directory to save our figures and datasets.
```python
dir.create("output") # create a new directory
```
### Data preprocessing steps
We'll use these objects for the rest of the analyses  
#### Removing singletons
For these examples we'll remove sequences that only occur once in the dataset since we cannot be sure they are biological units or caused by sequencing errors.  
```python
student_data_rs <- prune_taxa(taxa_sums(student_data) > 1, student_data)
```

For the next plots we'll scale the species counts into a percent abundance value. We'll also use this dataframe in other places.
```python
student_data_prop <- transform_sample_counts(student_data_rs, function(x) 100*x/sum(x))
```
The second is to rarefy our data by subsampling to an even depth. Note there are more robust methods to handle uneven datasets and we recommend checking out [this publication by McMurdie and Holmes](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531)
```python
data_rarefy <- rarefy_even_depth(student_data_rs, sample.size = min(sample_sums(student_data_rs)), verbose = FALSE, replace = TRUE)
```
### Generating a barplot of the dominant taxonomic groups
To start off with we will merge all OTUs per well and the summarize them by the phylum level
```python
student_data_prop_merged <- merge_samples(student_data_prop, "Well")
```

We need to re-apply names, before after merging, phyloseq will change the wells into level 1, 2, 3, not H41, H43 and H52 as before
```python
sample_data(student_data_prop_merged)$Well <- c("H41", "H43", "H52")
```

Then averge your data again (value should be 100% per sample)
```python
student_data_prop_merged_prop <- transform_sample_counts(student_data_prop_merged, function(x) 100*x/sum(x))
```

Merges species that have the same taxonomy at a certain taxaonomic rank.
```python
RANK="Phylum" # select Phylum as rank
student_data_prop_merged_prop_rank <- tax_glom(student_data_prop_merged_prop, taxrank=RANK)
```
Generate the barplot using the phyolseq function "plot_bar"
```python
plot_bar(student_data_prop_merged_prop_rank, fill=RANK) 
```

#### Subset and show the composition of the Actinobacteria
We'll subset our dataset to grab only bacterial OTU from the phylum **Actinobacteria**
```python
actino <- subset_taxa(student_data_prop_merged_prop, Phylum == "Actinobacteria")
```
We'll define our own vector of colors to plot by:
```python
set_colors <- c(
  "#DA5724",  "#508578", "#c481fd",  "#CD9BCD","#AD6F3B", "#652926", "#6dcff6",  "#C84248",  "#1e90ff", "#8569D5", 
  "#cc0953", "#D1A33D", "grey",  "pink", "gold","#8A7C64", "#599861", "navy", "#5F7FC7" , "tomato", "#673770",  
  "#008080", "#2F4F4F", "#FAEBD7", "#ff1493", "#5e738f","#808000", "#D14285", "#ffa500", "cbd588", "wheat", 
  "#d2b48c", "cyan2","black",  "#BC8F8F", "#800000","#008B8B",  "#BC8F8F"
)
```

Then we can generate the following bar plot:
```python
plot_bar(actino , x="Class", fill="Order", facet_grid=~Well) + ylab("Relative abundance (%)") + scale_fill_manual(values=set_colors)
```

#### Improving our first bar plot
We'll extract the abundance data for the abundance phyla (e.g. Relative abundance > 1%) and save the data
```python
student_data_prop_merged_prop_rank_filter <- filter_taxa(student_data_prop_merged_prop_rank, function(x) mean(x) > 1, TRUE)
abun.phylum <- psmelt(student_data_prop_merged_prop_rank_filter)
View(abun.phylum) 
```
If everything looks good, we can save this dataset as a text file
```python
write.table(abun.phylum, "output/phylum_over1percent.txt", sep="\t", dec=".", col.names=NA)
```

Then we'll calculate the sum of all remaining groups for an "Others" category.
```python
others <- tax_glom(student_data_prop_merged_prop_rank_filter, "Kingdom" )
otu_table(others) <- 100 - otu_table(others)
tax_table(others)@.Data[,2:6] <- "Others"# Define all taxanomic levels as "Others"
taxa_names(others) <- "Others" # Define the taxa name as "Others"
```

Next we can combine the abundant taxa with the sum of rare taxa
```python
OTU1 <- otu_table(cbind(otu_table(others), otu_table(student_data_prop_merged_prop_rank_filter) ), taxa_are_rows = FALSE)
TAX1 <- tax_table(rbind(tax_table(others), tax_table(student_data_prop_merged_prop_rank_filter) )) 
METADATA1 <- sample_data(student_data_prop_merged_prop_rank_filter)

final <- phyloseq(OTU1, TAX1, METADATA1)
```
Plot the abundance phyla using bar or heatmap
```python
plot_bar(final, fill=RANK) + 
  scale_fill_manual(values=set_colors)+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
```
And we can save the figure using the ggplot function, "ggsave"
```python
ggsave("output/taxonomic composition_phylum.png", width = 4, height = 4)
```

We could also generate a heatmap of the same data
```python
plot_heatmap(final, method ="PCoA", distance = "bray", sample.label="Well", taxa.label="Phylum", low="#000033", high="#FF3300")+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
# Note: by default plot_heatmap takes a log transformation (trans = log_trans(4))
ggsave("output/heatmap_phylum.png", width = 4, height = 4)
```

#### Extracting the top 30 OTUs and showing the abundance across samples
```python
top30 <- prune_taxa(names(sort(taxa_sums(student_data_prop),TRUE)[1:30]), student_data_prop)

plot_heatmap(top30, 
             # method ="NMDS", # default method
             method = NULL,
            # distance = "bray", 
            sample.label="Well", 
            taxa.label="Class", trans = NULL)+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0.5))

ggsave("output/heatmap_top30.png", width = 4, height = 4)
```

#### Phylogenetic tree 
We can use "Ape" and the provided phylogenetic tree to plot our abundant OTUs
```python
top20 <- prune_taxa(names(sort(taxa_sums(student_data_prop_merged_prop),TRUE)[1:20]), student_data_prop_merged_prop)

plot_tree(top20, color = "Phylum", shape = "Well", 
          #label.tips = "Family", 
         # size = "abundance", 
          plot.margin = 0.5, 
         ladderize = TRUE) +
  # coord_polar(theta="y")+
  scale_shape_manual(values = c(1,2,16))

ggsave("output/phylo_tree_phylum.png", width = 4, height =4)
```

## Exploring Alpha Diversity
Many richness estimates are modeled on singletons and doubletons in the abundance data. You need to leave them in the dataset if you want a meaningful estimate.
We'll be using the rarefied dataset we created in the "Data Preprocessing steps" above.
```python
student_data_rarefy <- data_rarefy
```

```python
plot_richness(student_data_rarefy) # warnings probably because of the graphics that didn't show clearly

p_alpha <- plot_richness(student_data_rarefy, x="Well",  
                          measures=c("Simpson", "Shannon"), 
                         color="Well") + geom_point(size=5, alpha=0.7) + scale_color_manual(values = c("red", "blue", "black"))

p_alpha 
```

### estimate alpha diversity
```python
alphadiv <- estimate_richness(student_data_rarefy
                              # , measures=c("Simpson", "Shannon")
)
alphadiv 
```

### Phylogenetic Diversity
Start off with a function that transposes the OTU matrix into a form that vegan can use.
```python
vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
```
To calculate Faith's phylogenetic distance we need a community data matrix and a phylo tree object
```python
FaithPD <- pd(vegan_otu(student_data_rarefy), phy_tree(student_data_rarefy))

## Returns a dataframe of the PD and species richness (SR) values for all samples
alphadiv2 <- data.frame(alphadiv, FaithPD)

write.table(alphadiv2, "output/alpha diversity.txt", sep="\t", dec=".", col.names=NA)

# Note: SR equals Observed using different functions
identical(alphadiv2$SR, alphadiv2$Observed)
```

Now you can go on with univariate data analysis, e.g. ANOVA or plot nice figures, For instance: 
```python
# combine the diversity results with the sample data file
sampledf <- data.frame(sample_data(student_data_rarefy))
df <- data.frame(alphadiv2, sampledf)
```

### Plot boxplot with significance
Add significance levels of the diversity indices between different wells
```python
my_comparisons <- list( c("H41", "H43"), c("H43", "H52"), c("H41", "H52") )
```

Phylogenetic diversity
```python
pd_well <- ggboxplot(df, "Well", "PD",
          color = "Well", palette =c("red", "blue", "black"),
          add = "jitter", 
          shape = "Aquifer")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + # pairwise comparison
  stat_compare_means(label.y = 50)      # Add global p-value      
```

Shannon diversity
```python
shannon_well <- ggboxplot(df, "Well", "Shannon",
          color = "Well", palette =c("red", "blue", "black"),
          add = "jitter", 
          shape = "Aquifer")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + # pairwise comparison
  stat_compare_means(label.y = 5.5)      # Add global p-value       
```

Plotting both together:
```python  
grid.arrange(pd_well, shannon_well, ncol=2)

dev.copy(png, file="alpha_diversity.png", width=1500, height=600, res=144)
dev.off()
#could also use ggsave() here
```

## Principal Component Analysis of Geochemical Variables
Read the data tables
```python
metadata_tb <- read.csv("data/simulated_dataset.txt", sep="\t", dec = ".", row.names=1, header=TRUE)
str(metadata_tb) # check internal structure of the metadata file
```

```python
metadata_tb.pca <- prcomp(metadata_tb[,c(7:21)], center = TRUE, scale. = TRUE)
# center = TRUE: a logical value indicating whether the variables should be shifted to be zero centered
# center = TRUE, scale. = TRUE: the variables is adviced to be scaled to have unit variance before the analysis takes place
```

To check the results of the PCA. The print method returns the standard deviation of each of the four PCs, and their rotation (or loadings), which are the coefficients of the linear combinations of the continuous variables.
```python
print(metadata_tb.pca)
```

The summary method describe the importance of the PCs. We can see there that the first two PCs accounts for more than 51% of the variance of the data.
```python
summary(metadata_tb.pca)
```

### Plotting with base R
The plot method returns a plot of the variances (y-axis) associated with the PCs (x-axis). The Figure below is useful to decide how many PCs to retain for further analysis. We can see that the first two PCs explain most of the variability in the data.
```python
plot(metadata_tb.pca, type = "l")
plot(metadata_tb.pca)
```

A simple PCA plot
```python
biplot(metadata_tb.pca, choices = 1:2, scale = 1)
```

Using ggbiplot
```python
ggbiplot(metadata_tb.pca)
```

Nice looking plot
```python
ggbiplot(metadata_tb.pca, ellipse=TRUE, 
         choices = 1:2, 
         var.axes = TRUE, # Draw arrows for the variables?
         obs.scale = 1, var.scale = 1, 
         labels=rownames(metadata_tb), groups=metadata_tb$Well) +
  scale_colour_manual(name="Well", values= c("red", "blue", "black"))+
  ggtitle("PCA of simulated dataset")+
  theme_bw()+
  theme(legend.position = "right")

ggsave("output/PCA of simulated dataset.png", width = 8, height = 6)
```

## Beta Diversity Ordination Plots
### Clustering and Dendrogram 

Calculate the distance matrix between samples
```python
dis <- phyloseq::distance(student_data_prop, "bray") # distance calculation
```

Hierarchical clustering
```python
hc1 <- hclust(dis, method="average")
hc2 <- hclust(dis, method="ward.D")
par(mfrow = c(2,1)) 
plot(hc1) 
plot(hc2)
```

Add colors to the tips
```python
sampledf <- data.frame(sample_data(student_data_prop))
well <- sampledf$Well

# palette<-hue_pal()(length(levels(well)))
palette <- c("red", "blue", "black" )
tipColoor <- col_factor(palette, levels=levels(well))(well)
```

Save it as phylo object
```python
clust.uf<-as.phylo(hclust(dis,method="ward.D2"))
# png("output/cluster.png")
plot(clust.uf,tip.color=tipColoor, direction="downwards") #10*20
dev.off()
```
More info on different clustering methods: <http://girke.bioinformatics.ucr.edu/GEN242/pages/mydoc/Rclustering.html>

### Unconstrained ordination 
Perform an ordination on phyloseq data
```python
set.seed(123)
```

NMDS plot based on Bray-Curtis distance
```python
ord <- ordinate(student_data_prop, "NMDS", "bray")
```
plot ordination
```python
p = plot_ordination(student_data_prop, ord, 
                    # type="biplot", 
                     type="split",
                    # type="samples",
                    # type="taxa",
                     shape="Well", color="Phylum", 
                     label="Well",
                     title="NMDS based on Bray-Curtis distance")+
                    # geom_point(size=7, alpha=0.75)+ 
                     theme_bw()
print(p)
ggsave("output/NMDS_bray.png", width = 8, height = 4)
```

#### PCoA plot based on Unifrac distance
```python
ordu = ordinate(student_data_prop, "PCoA", "unifrac", weighted=TRUE)
summary(ordu)
p1 = plot_ordination(student_data_prop, ordu, 
                     axes = 1:2,
                     type="samples", shape="Season", color="Well", 
                     #label="Well",
                     # type="taxa",color="Phylum", 
                     title="PCoA based on weighted UniFrac distance")+
  geom_point(size=7, alpha=0.75)+ 
  scale_colour_manual(values= c("red", "blue", "black"))+
  scale_shape_manual( values= c(0, 1,16,17))+
  theme_bw()

print(p1)
```
remove the inner dots
```python
p1$layers<-p1$layer [-1]
print(p1)

ggsave("output/PCoA_wunifrac.png", width = 8, height = 6)
```
Scree plot: shows the fraction of total variance in the data as explained or represented by each PC. 
```python
p_scree <- plot_ordination(student_data_prop, ordu, type="scree")
print(p_scree)
```
### Constrained ordination 
#### RDA
```python
ordu2 = ordinate(student_data_prop, "RDA", "unifrac", weighted=TRUE)
summary(ordu2)
p2 = plot_ordination(student_data_prop, ordu2, 
                     type="samples", 
                     # shape="Season",
                     color="Well", 
                     #label="Well",
                     # type="taxa",color="Phylum", 
                     title="RDA based on weighted UniFrac distance")+
  geom_point(size=7, alpha=0.75)+ 
  theme_bw()

print(p2)
```

Remove the inner dots
```python
p2$layers<-p2$layer [-1]
print(p2)
ggsave("output/RDA_wunifrac.png", width = 6, height = 4)
```
Scree plot: shows the fraction of total variance in the data as explained or represented by each PC. 
```python
p2_scree <- plot_ordination(student_data_prop, ordu2, type="scree")
print(p2_scree)
```

#### ccA plot based on Unifrac distance
```python
set.seed(123)
ordu3 = ordinate(student_data_prop, "CCA", "unifrac", weighted=TRUE)
summary(ordu3)
p3 = plot_ordination(student_data_prop, ordu3, 
                     axes = 1:2,
                     type="samples", shape="Well", color="Well", 
                     #label="Well",
                     # type="taxa",color="Phylum", 
                     title="CCA based on weighted UniFrac distance")+
  geom_point(size=7, alpha=0.75)+ 
  theme_bw()

print(p3)
# remove the inner dots
p3$layers<-p3$layer [-1]
print(p3)

# scree plot: shows the fraction of total variance in the data as explained or represented by each PC. 
p3_scree <- plot_ordination(student_data_prop, ordu3, type="scree")
print(p3_scree)
```

## Beta Diversity Hypothesis Testing
### PERMANOVA using adonis function in Vegan package
First calculate bray curtis distance matrix using either weighted unifrac distance or Bray-Curtis distance
```python
# student_data_prop_bray <- phyloseq::distance(student_data_prop, method = "unifrac", weighted=TRUE)

student_data_prop_bray <- phyloseq::distance(student_data_prop, method = "bray")
```

Make a data frame from the sample_data
```python
sampledf <- data.frame(sample_data(student_data_prop))
head(sampledf)
```

Then run the Adonis test (Number of permutations: 999)
```python
ht_well <- adonis(student_data_prop_bray ~ Well, data = sampledf) #  ***
```
This output tells us that our adonis test is significant so we can reject the null hypothesis that our three countries have the same centroid.

```python
ht_well
ht_well$aov.tab$"Pr(>F)"
write.table(data.frame(ht_well$aov.tab), "output/permanova_well.txt", sep="\t", dec=".", col.names=NA)
```

#### Posthoc Pairwise Adonis
This is an R wrapper function for multilevel pairwise comparison using adonis (~Permanova) from package 'vegan'. 
The function returns adjusted p-values using p.adjust().
First, we need to convert the phyloseq object into an OTU-table

```python
vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

pairwise.adonis(vegan_otu(student_data_prop), sampledf$Well)
```

!!! Your Turn
Test if there was any effect of other factors on bacterial community composition, e.g. Season

Hint:
```python
adonis(student_data_prop_bray ~ Season, data = sampledf)
# This output tells us that our adonis test is non-significant 
```

You can do the rest....

### Analysis of Similarities using anosim function in Vegan package
```python
ano_well <- anosim(student_data_prop_bray, sampledf$Well)
summary(ano_well)
par(mfrow=c(1,1))
plot(ano_well)
```

You can test if there was any effect of other factors on bacterial community composition, e.g. Season 
```python
ano_season <- anosim(student_data_prop_bray, sampledf$Season)
```

### Homogeneity of dispersion test using betadisper function in Vegan package
This is a nonmetric test based on permutations

```python
dis_well <- betadisper(student_data_prop_bray, sampledf$Well)
permutest(dis_well) # **
dis_season <- betadisper(student_data_prop_bray, sampledf$Season)
permutest(dis_season) # ns
```
Additionally, our betadisper results are significant, meaning we cannot reject the null hypothesis that our groups have the same dispersions. 
This means we cannot be confident that our adonis result is a real result, and not due to differences in group dispersions

### Tukey's Honest Significant Differences

```python
well.HSD <- TukeyHSD(dis_well)
well.HSD 
```
Plot dispersion distances between groups
```python
par(mfrow=c(1,2))
plot(well.HSD, las=1)
plot(dis_well, las=1)
```
## Linking to environmental parameters and biplots
We'll use OTUs with mean relative abundance over 0.01% for CCA analysis
```python
student_data_prop_filter <- filter_taxa(student_data_prop, function(x) mean(x) > 0.01, TRUE)
```

Import the OTU table and geochemical data:
```python
otutable <- vegan_otu(student_data_prop_filter)
sampledf <- data.frame(sample_data(student_data_prop_filter))
```

To check whether to choose CCA or RDA, we should use Detrended correspondence analysis (DCA)
```python
dca <- decorana(otutable)
dca 
```
Check the "Axis lengths" of the DCA1 (Important)
```python
summary(dca) # to check the "Axis lengths" of the DCA1 (Important)
```
If this value < 3, it is better to use RDA
If this value > 4, it is better to use CCA
If this value is between 3 and 4, either use CCA or RDA


###  Envfit + CCA/RDA
Much more powerful when more complex factors are tested
#### CCA
Scale species to unit variance
```python
df <- data.frame(scale(sampledf[,-c(1:6)]))
```

Test all environmental factors
```python
ord_all <- cca(otutable ~ ., data=df) 
```
You can choose to test only some of the environmental factors
```python
ord_select <- cca(otutable ~ EC + DO + PH, df) # 
```

#### RDA
```python
bray <- vegdist(otutable, "bray")
ord_all <- capscale(bray~., df])
```

Plot all environmental factors
```python
plot(ord_all, type = "p", scaling = "sites") 
```

### Variance Inflation Factors
Linear dependencies between constraints can be investigated via the variance inflation factor or VIF
VIF is a measure of how much the variance of $\hat{\beta}_j$ is inflated by presence of other covariates

Lots of rules of thumb  
VIF >= 20 indicates strong collinearity in constraints  
VIF >= 10 potnetially of concern & should be looked at  

They will be completely removed from the estimation, and no biplot scores or centroids are calculated for these aliased constraints. 

```python
vif.cca(ord)
temp <- vif.cca(ord_all)
temp
select_para <- names(temp[temp < 10])
select_para
```
Keep the environmental factors with VIF less than 10
```python
ord <- cca(otutable ~ ., df[,select_para]) 
```

### Fit environmental vectors onto the ordination
The function fits environmental vectors or factors onto an ordination. 
The projections of points onto vectors have maximum correlation with corresponding environmental variables, and the factors show the averages of factor levels.

```python
fit <- envfit(ord, df[,select_para], perm = 999, display = "lc", scaling = "sites")
fit$vectors # check the significance
```

Extract the best (significant) variables (p < 0.05)
```python
spp.scrs <- as.data.frame(scores(fit, display = "vectors"))
spp.scrs
pval <-  fit$vectors$pvals
pval 
```

Data for the envfit arrows
This is necessary, see Gavin Simpson http://stackoverflow.com/questions/14711470/plotting-envfit-vectors-vegan-package-in-ggplot2/25425258#25425258
```python
fdat <- cbind(spp.scrs, Vector = rownames(spp.scrs), pval)
```

Now select only the signifant factors
```python
bestEnvVariables<-rownames(fdat)[fdat$pval<=0.05]
```
If you have NA entries, remove them
```python
# bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]
```

Redo CCA using the best environmental variables
```python
eval(parse(text=paste("ord1 <- cca(otutable ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=sampledf)",sep="")))
summary(ord1)
```
You can also use anova.cca to select the significant variables
```python
anova.cca(ord1, perm=9999)
anova.cca(ord1, by="margin", perm=9999) # marginal effects of the terms 
anova.cca(ord1, by="terms", perm=9999) # sequential
anova.cca(ord1, by="axis") # axis, slow
drop1(ord1, test="perm")
```
Now re-fit the environmental factors on ordination
```python
fit1 <- envfit(ord1,sampledf[,bestEnvVariables], perm = 999, display = "lc", scaling = "sites")
```
Simple triplot
```python
plot(ord1)
```
Simple biplot
```python
plot(ord1, type="n")
```

Plot the samples
```python
points(ord1, display = "sites", 
       col = as.numeric(Moisture),
       pch=16)
```
Choose the significant environmental factors
```python
plot(fit1, col = "red", cex=1.2, axis=TRUE, p.max = 0.05)
summary(ord1)
```
#### Plot using ggplot2
```python
spp.scrs <- data.frame(scores(fit1, display = "vectors"))
pval <-  fit1$vectors$pvals

#data for the envfit arrows

spp.scrs <- cbind(spp.scrs, Vector = rownames(spp.scrs), pval) # vector table
scrs <- as.data.frame(scores(ord1, display = "sites")) # sample table
scrs <- cbind(scrs, sampledf)

spp.scrs1<- subset(spp.scrs, pval<=0.05) #extracts relevant environment vectors from envifit
```

Biplot
```python
library(digest)
p <- ggplot(scrs) +
  geom_point(mapping = aes(x = CCA1, y = CCA2, colour = Well), alpha = 0.8, size = 8) +
   #coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = spp.scrs1,
               aes(x = 0, xend = CCA1*2.5, y = 0, yend = CCA2*2.5),
               arrow = arrow(length = unit(0.25, "cm")), colour = "blue", size=1) +
  geom_text(data = spp.scrs1, aes(x = CCA1*3, y = CCA2*3, label = Vector),
            size = 4) +
  theme_bw() +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14)) +
  #scale_color_manual(values = set_colors3, guide = guide_legend(ncol=1)) +
  scale_color_manual(values =c("red", "blue", "black"), guide = guide_legend(ncol=1)) +
  #xlim(-2,3)+
  theme(legend.position="right") 

p 

ggsave("output/cca_linking_env.png", width = 6, height = 4)
```

!!! Question  

ggplot: Triplot?
Hint:extract species scores
```python
spe <- scores(ord1, display = "species") 
tax <- tax_table(student_data_prop_filter) 
otu <- otu_table(student_data_prop_filter) 
abundance <- rowMeans(x=otu)# calculate mean abundance of each OTU 
df2 <- data.frame(spe, tax, Abundance = abundance) 
df2 <- subset(df2, abundance >=0.05) 
```
More info on model selection: https://github.com/naupaka/esa_vegan/blob/master/03-constrained-ordination/constrained-ordination.md

## Bioindicators
```python
biocLite("DESeq2")
biocLite("IRanges")
library(IRanges)
library(DESeq2)
library(ggpubr)
```

For this analysis, we use the count data. because this method does its own normalization
If we only want to identy the taxa that discriminate H43 and H52
```python
student_data_well <- subset_samples(student_data, Well!="H41")
```
### DESeq2
The following two lines actually do all the complicated DESeq2 work. 
The function phyloseq_to_deseq2 converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated using the experimental design formula, also shown (the ~Well term). 
```python
dds = phyloseq_to_deseq2(student_data_well, ~ Well)
```
Then DESeq function does the rest of the testing
    - estimation of size factors: estimateSizeFactors
    - estimation of dispersion: estimateDispersions
    - Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest
Wald statistics: nbinomWaldTest
```python
dds <- DESeq(dds, test="Wald", fitType="parametric")
```
Negative Binomial GLM fitting and 
```python
ddsLRT <- DESeq(dds, test="LRT", reduced= ~ 1)
```
Check our results
```python
res = results(dds, cooksCutoff = FALSE)
alpha = 0.001
sigtab = res[which(res$padj < alpha), ]
# sigtab = sigtab[which(abs(sigtab$log2FoldChange) > 2), ] # sometimes you want to report the taxa that had higher log2FoldChange than 2

sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(student_data_well)[rownames(sigtab), ], "matrix"))
head(sigtab)
```

#### Plot
Phylum order
```python
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
```
Genus order
```python
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
```
Ploting (Note: H43 is used as a control, FoldChange = H52/H43)
```python
sigtab$OTU <- rownames(sigtab)
ggbarplot(sigtab, x = "OTU", y = "log2FoldChange",
          fill = "Phylum",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in ascending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          ylab = "log2FoldChange",
          legend.title = "Phylum",
          rotate = TRUE,
          ggtheme = theme_classic()
)

ggsave("output/bioindicators.png", width = 6, height = 3)
```
Show the abundance difference of the interesting taxa between the two wells
```python
select <- rownames(sigtab)
student_data_well_select <- subset_taxa(student_data_well, rownames(tax_table(student_data_well)) %in% select[2])

plot_bar(student_data_well_select, "Well", fill="Genus") +
  ylab("Sum of normalized abundance") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0, hjust = 1))

ggsave("output/bioindicators_bar.png", width = 6, height = 3)

plot_heatmap(student_data_well_select, sample.label="Well", taxa.label="Genus", taxa.order="Phylum", sample.order="Well")
```

### Random Forest
<https://rstudio-pubs-static.s3.amazonaws.com/115631_7397b7cf67534479ae80f70546610eea.html>

<sub>Written by Lijuan Yan</sub>
<sub>Oct 2018</sub>