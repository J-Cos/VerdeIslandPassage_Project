# 1. Packages and path to project directory ######
    set.seed(0.1)
    library(tidyverse)
    library(DECIPHER)
    library(phyloseq)

    #computational parameters
        nproc<-4
    #1.1 load functions
        #functions we write for this project
            function_files<-list.files("Functions")
            sapply(file.path("Functions",function_files),source)
        #functions we need from the Bioinformatic Pipeline Repository (don't worry about this repository you won't have to work directly with it)    
            function_files<-list.files("../../BioinformaticPipeline_Env/BioinformaticPipeline/SupportFunctions")
            sapply(file.path("../../BioinformaticPipeline_Env/BioinformaticPipeline/SupportFunctions",function_files),source)

#2. load data
    #load main ps
        ps<-SeqDataTable2Phyloseq(  SeqDataTablePath="../Data/Dummy_COI_SeqDataTable.RDS", #this is the final function from the bioinformatic pipeline
                                    clustering="ESV", 
                                    Metadata=NULL, 
                                    assignment="Idtaxa", 
                                    ClusterAssignment=NULL)

    #load metadata once we have it
        metadata<- read.csv(file.path("../Data/Metadata.csv") , header=TRUE,  sep=",") %>%
            as_tibble %>%
            mutate_if(is.character, as.factor) %>% 
            as.data.frame
        #repair rownames as they were removed by making into a tibble
        rownames(metadata)<-metadata$SampleID

    #incorporate metadata into physeq object
        sample_data(ps)<-metadata

#3. Clean data
    #prune samples without metadata or with too few reads
        ps<-ps %>%
            prune_samples(sample_sums(.) > 50000, .) %>% # no samples with fewer than 50,000 reads
            prune_taxa(taxa_sums(.)>0, .) # remove any ESVs that now don't appear in any samples

#4. Save data
    saveRDS(ps, file="../Outputs/ps.RDS")
