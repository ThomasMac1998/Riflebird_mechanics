### Code for Analyses of Manus Angle 
### Thomas MacGillavry 
### 23.01.2024

## load packages 
library(phytools)
library(geiger)
library(dplyr)
library(viridis)
library(ggtree)
library(ggplot2)
library(scico)

## Set working directory 
setwd("")

## Import dataset 
data <- read.csv("Baliga.2019.ManusMax.csv")
head(data, 4)

## Specify row names so the data matches the tree
rownames(data) <- data$species_phylo

## Read the tree file
## The concatenated tree, already ultrametricized in Geneious 
tree <- read.nexus("ManusTree_Consesus.nex")
print(tree, printlen = 2)

## First, let's plot the tree 
plotTree(tree,fsize=0.9,ftype="i",lwd=1)

## Randomly resolve the one polytomy in Falco 

# Load the package 
library(RRphylo)

# Use the fix.poly() function to randomly resolve polytomies. 
tree <- fix.poly(tree,type="resolve")
is.ultrametric(tree)
plot(tree, type = "fan") # The tree is now bifurcating and ultrametric 

## Convert to a vector 
xx<-setNames(data$manusAngle.max,rownames(data))
xx

## Now we want to reconstruct ancestral states, and plot these results on a phylogeny. 

## Start again by plotting the tree: 
plotTree(tree,type="fan",ftype="i")

## Read in data frame again 
df <-read.csv("Baliga.2019.ManusMax.csv",row.names=1)
xx <- select(df, manusAngle.max) # Select only brain mass and species variable 

## Change data frame to a vector
xx<-as.matrix(xx)[,1]
xx

## Estimate ancestral states 
fit <- fastAnc(tree, xx, vars=TRUE, CI=TRUE)
fit

## We can also calculate 95% CIs 
fit$CI[1,]
range(xx)

### Plot again in ggtree ### 

# Fit an ancestral state character reconstruction
fit <- phytools::fastAnc(tree, xx, vars = TRUE, CI = TRUE)

# Make a dataframe with trait values at the tips
td <- data.frame(
  node = nodeid(tree, names(xx)),
  trait = xx)

# Make a dataframe with estimated trait values at the nodes
nd <- data.frame(node = names(fit$ace), trait = fit$ace)

# Combine these with the tree data for plotting with ggtree
d <- rbind(td, nd)
d$node <- as.numeric(d$node)
tree <- full_join(tree, d, by = 'node')

# Plot a Phenogram: 

# Create the plot with the entire tree but only display tip labels for selected species
ggtree(tree, aes(color = trait), continuous = 'colour', yscale = "trait") +
  theme_classic() + # Choose theme 
  scale_x_continuous(name = "Evolutionary distance", breaks = seq(0, 100, by = 25)) + # Set margins
  scale_y_continuous(name = "Maximum manus angle", limits = c(100, 240), breaks = seq(100, 240, by = 25)) + 
  theme(
    plot.margin = margin(1, 1, 1, 1, "cm"), # Adjust the margins as needed
    axis.text.x = element_text(size = 20),   # Set x-axis tick label size
    axis.text.y = element_text(size = 20),   # Set y-axis tick label size
    axis.title.x = element_text(size = 20),  # Set x-axis title size
    axis.title.y = element_text(size = 20), # Set y-axis title size
    legend.position = c(0.05, 0.95),               # Move legend to the top left
    legend.justification = c(0, 1)           # Justify legend to the top left
  ) +
  labs(color = "Manus angle") +  # Change the legend title
  scale_color_continuous(breaks = seq(0, 240, by = 25)) + # Adjust the breaks on the legend scale 
  scale_color_scico(palette = "berlin") # https://cran.r-project.org/web/packages/scico/scico.pdf 

### Now we want to plot the data normally to better visualize the values of each species 

# Remove underscores and replace them with spaces
data$species_phylo <- gsub("_", " ", data$species_phylo)

# Order the levels of species_phylo based on manusAngle.max
data$species_phylo <- reorder(data$species_phylo, data$manusAngle.max)

# Now, create your ggplot with the modified data
ggplot(data = data, aes(x = species_phylo, y = manusAngle.max, color = manusAngle.max)) +
  theme_classic() + 
  geom_segment(aes(xend = species_phylo, yend = 100), color = "grey90", size = 2) + 
  geom_hline(yintercept=180, linetype="dotted", color = "grey70", size=1) + 
  geom_point(size = 5, shape = 20) +
  ylim(100, 240) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_scico(palette = "berlin") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none") +    # Remove legends
  xlab(NULL) + 
  theme(
    plot.margin = margin(1, 1, 1, 1, "cm"), # Adjust the margins as needed
    axis.text.x = element_text(size = 10),   # Set x-axis tick label size
    axis.text.y = element_text(size = 8)) + 
  scale_y_continuous(name = "Maximum manus angle", limits = c(100, 240), breaks = seq(100, 240, by = 25)) 

##############################################################################################################################
##############################################################################################################################

## Some STDEV analyses: 

# Filter out rows for Ptiloris victoriae and Ptiloris magnificus
filtered_data <- subset(data, species_phylo != "Ptiloris victoriae" & species_phylo != "Ptiloris magnificus")

# Calculate standard deviation for manusAngle.max in the filtered data
sd(filtered_data$manusAngle.max) # 16.06149
mean(filtered_data$manusAngle.max) # 149.9533

# For Victoria's
(237.0670-149.9533)/16.06149
# 5.423762 s.d. above mean for control species! 

# For Magnificent 
(238.6150-149.9533)/16.06149
# 5.520142 s.d. above mean for control species! 







