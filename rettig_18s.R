### Below are chunks of annotated code for (some) analysis of sequencing data
### Edit for your file names, modify as you see fit. 

### Also, though exploratory analysis can be helpful and interesting-- for the time constraints of the summer 
### make sure you don't use up all of your time just making preliminary figures instead of asking a specific question 
### and then answering it using the data
### Such as "How does syndiniales diversity vary over time?" then pursue how you would answer that... 
### maybe by estimating multiple diversity metrics, plotting over time, and conducting a statistical test. Just an example.

### Load Files and Set Working Directory
# Set your working directory to wherever you want all of your figures and files to save to
setwd("~/Desktop/surfo_project")

### Data and sample information
counts <- read.table("syn_asv_counts.tsv", header=T, row.names=1,check.names=F, sep="\t")
taxa <- as.matrix(read.table("syn_asv_taxonomy.tsv",header=T, row.names=1,check.names=F, sep="\t"))
samples <- read.table("syn_sample_info.tsv", header=T, row.names=1,check.names=F, sep="\t")
# View(counts) #Useful to know what you're working with and confirm the files imported correctly.
# View(taxa)
# View(samples)

##### Load Libraries ####
# You might not need all of these... install what is useful
# install.packages()

library("devtools")
library("dada2")
library("ggplot2")
library("phyloseq")
library("vegan")
library("DESeq2") # might not need
library("dendextend") #might not need
library("pBrackets") #might not need
library("tidyr")
library("viridis") #colors
library("reshape")
library("ape") # can't remember what this was for...
library("ggrepel")
library("tibble")
library("ggsci") # colors
####

#### Phyloseq Object ####
# A phyloseq object is a modifiable object that contains all of the ASV counts, taxonomy, and sample information
# information on phyloseq here: https://joey711.github.io/phyloseq/
# You can also work with just the extracted proportions table, but a phyloseq obejct allows you to switch around 
# what taxonomic level you are looking at, which can be useful.

count_tab <- otu_table(counts, taxa_are_rows=T)
tax_tab <- tax_table(taxa)
sample_info <- sample_data(samples)
ASV_physeq <- phyloseq(count_tab, tax_tab, sample_info)

#### 
#### Beta Diversity ####

#### Hierarchical Clustering Dendrogram ####
# A phylogenetic tree can be informative in looking at the relatedness of samples (in terms of their ASV composition)
# You can specify a specific subset of taxonomic lineage 
# You can inforporate sample information like depth, light level, etc. 
# See this tutorial: https://joey711.github.io/phyloseq/
# http://joey711.github.io/phyloseq/plot_tree-examples.html#example

# An example where GP is the phyloseq object and the tree produced uses the counts (abundance) for point size, labels from genus, point shape as family, and point color by sampel type
GP.chl <- subset_taxa(GP, Phylum=="Chlamydiae")
plot_tree(GP.chl, color="SampleType", shape="Family", label.tips="Genus", size="Abundance")

##### Richness and Beta Diversity ####
# "Richness" is simply the quantity of species/taxa/ASVs whereas "evenness" would look at the distribution and relative proportions of each species/taxa/ASV
# https://www.bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-analysis.html
# https://grunwaldlab.github.io/analysis_of_microbiome_community_data_in_r/07--diversity_stats.html
# https://joey711.github.io/phyloseq/plot_richness-examples.html
# Other ways to explore "diversity" --> prevalent taxa, composition (presence/absence), rare taxa

# Example code below
# using phyloseq object, julian day (a more numerical Date), sample depth, and testing different richness metrics- chao1 and shannon, google them!
plot_richness(ASV_physeq, x="julianD", color="Depth",  measures=c("Chao1", "Shannon")) + 
  geom_point(size=2) +
  scale_color_gradient(low = "#56B1F7", high = "#132B43") +
  scale_fill_gradient(low = "#56B1F7", high = "#132B43") +
  theme(legend.title = element_blank())
# usually nicer to make this into boxplots 

#### Syndiniales Community Compositions ####
# There are lots of ways to show the community composition: barplot, pie, heatmap, area plot, etc. 
# This is just how I quickly make a barplot of community composition for a subset of samples

# Below, I am using my phyloseq object "autotrophs" and pulling out only the samples that meet the 
# criteria "ship== JC214 and depth==5" (meters) and this way I can look at JUST the surface phytoplankton community 
# You can do similar things if you wanted to pull out surface parasite samples instead of looking at the whole shebang.

# JC214 Surface Samples
surface_auto_jc <- subset_samples(autotrophs, ship=="JC214" & Depth==5)
surface_auto_jc <- subset_samples(surface_auto_jc, type=="Vertical" | type=="Light") #similarly i am subsetting again to use specific types of smaples from certain experiments
surface_auto_counts <- otu_table(tax_glom(surface_auto_jc, taxrank="Genus"),taxa_are_rows) #this extracts a counts table at a specific taxonomic level, in this case genus-- so it collapses things down so you aren't looking at full taxonomy
surface_auto_tax <- as.vector(tax_table(tax_glom(surface_auto_jc, taxrank="Genus"))[,7]) #you do a similar extraction for the taxonomy table to create a taxonomic vector that will be combined with the counts table
rownames(surface_auto_counts) <- surface_auto_tax #this assigns the taxa vector as the rownmaes for the counts, so now you have one data table of samples as column names, genus as row names, and the ASV counts 

Pruned = filter_taxa(surface_auto_counts, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
# This pruning removes taxa that are not seen more than 3 times in at least 20% of samples; this is optional.
# I do this to weed out the ASVs that are so rare they're not ecologically informative

surface_auto_df <- Pruned # renaming
prop1 <- apply(surface_auto_df, 2, function(x) x/sum(x)*100) # converting the counts to relative adundance (%) 
View(prop1) #always check that it worked as intended

## Filter Taxa that compose < 1%
# This is a subsequent step to put any remaining low proportion taxa into a category called "other" 
# Mostly this keeps plots tidy without removing data

temp_filt_prop <- data.frame(prop1[apply(prop1, 1, max) > 1, ])
# this is everything above the filter

filt_prop1 <- colSums(prop1) - colSums(temp_filt_prop)
#this should leave filt_prop1 as the < 1% taxa....

prop1_surf <- rbind(temp_filt_prop, "Other"=filt_prop1)
# then recombine the "other" with the  proportions table

df1 <- prop1_surf #renaming

surf_df <- as.data.frame(df1) #ensure the table is a dataframe
surf_df <- surf_df %>% rownames_to_column("Genus") #make the taxa rownames a column
surf_df_long <- gather(surf_df, Sample, Proportion, -Genus) #this transforms data into longformat, which is easier for plotting
head(surf_df_long)
## this is the completed relative abundance table with all below threshold taxa compiled into "Other", using a 1% filter

sample_info_for_merge <- data.frame("Sample"=row.names(sample_info), "cruiseID"=sample_info$cruiseID,"julian"=sample_info$julianD,"epoch"=sample_info$epoch, stringsAsFactors=F) 
#this is merging the sample information data with the long format dataframe, you can specify what sample 
#information you want added to the dataframe and name it whatever you prefer i.e. "julian" is the new name for the julianD extracted from the sample info file
autotrophic_community <- merge(surf_df_long, sample_info_for_merge) # merging info with data and reassigning name
# for ease of plotting, make sure any metadata is read as a factor, otherwise you get errors
autotrophic_community$epoch <- as.factor(autotrophic_community$epoch) 
autotrophic_community$cruiseID <- as.factor(autotrophic_community$cruiseID)
autotrophic_community$julian <- as.factor(autotrophic_community$julian)
autotrophic_community <- autotrophic_community[order(autotrophic_community$julian),] #this is re-ording the data, so that when i plot it reads it in julian day (date) order
autotrophic_community$cruiseID <- factor(autotrophic_community$cruiseID, levels = unique(autotrophic_community$cruiseID[order(autotrophic_community$julian)]))
# this is a necessary step to re-order the actual sample IDs that show on the plot in association with the preferred data order

# Now I am makingg a barplot with the relative abundance of samples ordered by julian day (even though they aren't the x labels)
ggplot(autotrophic_community, aes(x=cruiseID, y=Proportion, fill=Genus)) +
  geom_bar(width=1, stat="identity") +
  guides(fill = guide_legend(ncol = 2))+
  scale_fill_manual(values = c(pal1))+ # I made a color palette to use, but you can use default colors if you omit this
  #scale_colour_iwanthue() +
  #scale_fill_iwanthue(hmax = 295, cmin = 15, cmax = 70, lmin = 30, lmax = 70) + #an alternative color option
  theme_bw() + 
  scale_y_continuous(limits = c(0,101), expand = c(0, 0))+
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.key.size = unit(1,"line"), 
        legend.text = element_text(size = 7), legend.title=element_text())+
  labs(x= "Sample", y="Relative 18S Abundance (%)", fill = "Genus")

#### Ordination ####
### Read up a little on using ordination methods, there are numerous and some are more 
# appropriate than others and different methods will tell you different things
#
# generating and visualizing the PCoA with phyloseq
pcoa <- ordinate(ASV_physeq, method="MDS", distance="euclidean")
eigen_vals <- pcoa$values$Eigenvalues 

plot_ordination(ASV_physeq, pcoa, color="layer_epoch") +
  scale_color_manual(values=c("#FFDB6D", "#56B4E9", "#C4961A", "#4E84C4", "#D16103", "#293352")) +
  scale_fill_manual(values=c("#FFDB6D", "#56B4E9", "#C4961A", "#4E84C4", "#D16103", "#293352")) +
  geom_point(aes(colour=factor(layer_epoch))) +
  #geom_point(aes(shape=factor(type)),size=3) + 
  labs(col="layer_epoch") + 
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(group = depth_layer))+
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + 
  theme_bw()+
  theme(legend.title = element_blank())

# generating and visualizing the PCoA with phyloseq
syn_nmds <- ordinate(syn_pruned, method="NMDS", distance="jaccard")
syn_nmds$stress #0.06 good fit
plot_ordination(syn_pruned, syn_nmds, color="layer_epoch") +
  scale_color_manual(values=c("#FFDB6D", "#56B4E9", "#C4961A", "#4E84C4", "#D16103", "#293352")) +
  scale_fill_manual(values=c("#FFDB6D", "#56B4E9", "#C4961A", "#4E84C4", "#D16103", "#293352")) +
  geom_point(aes(colour=factor(layer_epoch))) +
  geom_point(size=3)+
  labs(col="layer_epoch") + 
  theme_bw()+
  theme(legend.title = element_blank())

