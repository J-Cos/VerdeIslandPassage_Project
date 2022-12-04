#plot_ordination depends on distance and ordinate functions
library("phyloseq"); packageVersion("phyloseq")
#below is R inbuilt 16S bacteria data
data("GlobalPatterns")
?packageVersion
#packageVersion function tells you info about the package e.g what version
library("ggplot2"); packageVersion("ggplot2")
#old version of dplyr for data wrangling
library("plyr"); packageVersion('plyr')

theme_set(theme_bw())
#filter low occurrence OTU's as noise variables
#preprocessing shows patterns in data
#however here we focus on removing OTU's that do not appear mor than 5 times in half the samples
print(GlobalPatterns)
#first filter/ pre process OTU's in GP1
GP<-GlobalPatterns
?filterfun_sample#from bioconducter repository for taxa filtering 
#function(x) also like the apply function. Bascially applies the function to all numbers in the data set 
wh0<- genefilter_sample(GP, filterfun_sample(function(x) x>5), A=0.5*nsamples(GP))
#prune_taxa removes (prunes) unwanted OTU's from phylgenetic data
GP1<- prune_taxa(wh0, GP)

#transform to even sampling depth
#sampling depth is the number of sequenced bases for a given sample
GP1 = transform_sample_counts(GP1, function(x) 1E6 * x/sum(x))
#keep only the most abundant 5 phyla
phylum.sum = tapply(taxa_sums(GP1), tax_table(GP1)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
GP1 = prune_taxa((tax_table(GP1)[, "Phylum"] %in% top5phyla), GP1)

#this still leaves us with 204 OTU's in the data set GP1

#some of the microbiomes are human associated and some aren't so need to seperate these
human = get_variable(GP1, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue")
#store in column human the category human
sample_data(GP1)$human <- factor(human)

#four main ordination plots#############
#plot_ordination function supports 4 basic representations of an ordination
#here plot OTU's and shade points by phylum
GP.ord <- ordinate(GP1, "NMDS", "bray")
p1 = plot_ordination(GP1, GP.ord, type="taxa", color="Phylum", title="taxa")
print(p1)
#lots of occlusion here so lots of points get in the way of understanding actual data
#facet_wrap makes long rows of panels to display graphs
p1+facet_wrap(~Phylum,3)
#next plot the samples and shade the points by sample type

p2 = plot_ordination(GP1, GP.ord, type="samples", color="SampleType", shape="human") 
p2 + geom_polygon(aes(fill=SampleType)) + geom_point(size=5) + ggtitle("samples")

#####bioplot grpahic
#plot ordination crates 2 different graphics where sample type and OTU's plotted 
#together in one big biplot

p3 = plot_ordination(GP1, GP.ord, type="biplot", color="SampleType", shape="Phylum", title="biplot")
# Some stuff to modify the automatic shape scale
GP1.shape.names = get_taxa_unique(GP1, "Phylum")
GP1.shape <- 15:(15 + length(GP1.shape.names) - 1)
names(GP1.shape) <- GP1.shape.names
GP1.shape["samples"] <- 16
p3 + scale_shape_manual(values=GP1.shape)


#####split graphic
#Previous graph there was a big issuee with occlusion
#type="split" useful as OTU's and samples on different graphs in different panels
p4 = plot_ordination(GP1, GP.ord, type="split", color="Phylum", shape="human", label="SampleType", title="split") 
p4

#now make the sample colors black
gg_color_hue <- function(n){
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
color.names <- levels(p4$data$Phylum)
p4cols <- gg_color_hue(length(color.names))
names(p4cols) <- color.names
p4cols["samples"] <- "black"
p4 + scale_color_manual(values=p4cols)

#Supported ordination methods##############
#here store the plot results in a list then plot these results in 
# a combined graphic using ggplot2
#using the different method parameter functions in the plot_ordination function

dist = "bray"
ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, physeq, dist){
  ordi = ordinate(physeq, method=i, distance=dist)
  plot_ordination(physeq, ordi, "samples", color="SampleType")
}, GP1, dist)

?llply
?lapply
#lapply- retunrs list of same length of x, each result is from applying  a function to each element of x
#llply same as lapply except it preserves the labels in the list

?ord_meths
names(plist)<-ord_meths

#code below extracts data from each individual plit and puts it back together
#in one big data frame so then all the plots are included in one big plot 
#ldply takes in a list, will apply a function then export the results into a datafram
#laply does the same but not into a data frame, into a data array
pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"

#above, can see using ldply thar all ordination points combined into one data frame called pdataframe
#pdata frame then used to make a standard faceted ggplot scatterplot

p = ggplot(pdataframe, aes(Axis_1, Axis_2, color=SampleType, shape=human, fill=SampleType))
p = p + geom_point(size=4) + geom_polygon()
p = p + facet_wrap(~method, scales="free")
p = p + scale_fill_brewer(type="qual", palette="Set1")
p = p + scale_colour_brewer(type="qual", palette="Set1")
p

########MDS (principal co-ordinate analysis) on UniFrac distances######
#unifrac (unique fraction metric) measures the phylogenetic distance between sets of taxa in a phylogenetic tree
#this is as a fraction of the whole branch length of the tree

#use ordinate function to do a unifrac the do a principle co-ord analysis

ordu = ordinate(GP1, "PCoA", "unifrac", weighted=TRUE)

#plot_ordination then creates a ggplot graphic of this

plot_ordination(GP1, ordu, color="SampleType", shape="human")

#add some stuff to make the graph look nicer

p = plot_ordination(GP1, ordu, color="SampleType", shape="human")
p = p + geom_point(size=7, alpha=0.75)
p = p + scale_colour_brewer(type="qual", palette="Set1")
p + ggtitle("MDS/PCoA on weighted-UniFrac distance, GlobalPatterns")

