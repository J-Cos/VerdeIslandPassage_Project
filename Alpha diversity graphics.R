library("phyloseq"); packageVersion("phyloseq")
data("GlobalPatterns")
library("ggplot2"); packageVersion("ggplot2")
theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)
}

#first need to prune OTU's that are not present in any of the samples

GP <- prune_species(speciesSums(GlobalPatterns) > 0, GlobalPatterns)

#here plot default graphic for species richness using plot_richness function on GP
plot_richness(GP)

#instead og interpreting sample names directly like in previous, can organise using the sample type
#specify sample variable what to organise along the horizontal axis

plot_richness(GP, x="SampleType", measures=c("Chao1", "Shannon"))

#might want an external varibale thats not in the dataset e.g whether the samples are human or not
#first define new variable human as a factor 
sampleData(GP)$human <- getVariable(GP, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue")

#now tell plot_richness to map new human variable on horizontal axis, 
#shade points in different colour groups accroding to which sample type they belong to
plot_richness(GP, x="human", color="SampleType", measures=c("Chao1", "Shannon"))


#can then merge the samples using column sample type so data points are not 
#spread in the graph

GPst = merge_samples(GP, "SampleType")
# repair variables that were damaged during merge (coerced to numeric)
sample_data(GPst)$SampleType <- factor(sample_names(GPst))
sample_data(GPst)$human <- as.logical(sample_data(GPst)$human)

#now we can plot the environment merged data
#first store the default ggplot graphic as p then additional geom_point layer

p = plot_richness(GPst, x="human", color="SampleType", measures=c("Chao1", "Shannon"))
p + geom_point(size=5, alpha=0.7)


#more details about ggplot2

#check the list which is present in P
p$layers

#add new geom plot with the following point size

p$layers <- p$layers[-1]
p + geom_point(size=5, alpha=0.7)



