#------------------------------
# MAKE 18S VIP PHYSEQ
#------------------------------
# 1. Dependencies and Packages and path to project directory ######
    set.seed(0.1)
    #requires bioinformatic pipeline functions and crosspacific data
    library(tidyverse)
    library(DECIPHER)
    library(phyloseq)

    #1.1 load functions
        #functions we need from the Bioinformatic Pipeline Repository (don't worry about this repository you won't have to work directly with it)    
            function_files<-list.files("../../BioinformaticPipeline_Env/BioinformaticPipeline/SupportFunctions")
            sapply(file.path("../../BioinformaticPipeline_Env/BioinformaticPipeline/SupportFunctions",function_files),source)

#2. Parameters
    nproc<-4
    MakeTree<-FALSE

#3. load data from cross pacific folder
    ps18<-readRDS(file=file.path("/home/j/Dropbox/CrossPacific_Paper/Outputs/ps18.RDS"))

#4. subset to VIP
    psVIP<-ps18 %>%             
        prune_samples( sample_data(.)$Region=="Philippines" , .) %>%
        prune_taxa(taxa_sums(.)>0, .)

#5. create additional physeq object components if desired
    if (MakeTree){
        # align reference sequences for easier downstream sequence-based analyses
            AlignedSeqs<-refseq(psVIP) %>%
                DECIPHER::AlignSeqs(., processors = nproc)
        #create phylogenetic tree
            tree<-DECIPHER::TreeLine(AlignedSeqs, reconstruct=TRUE, maxTime=2) # default is method="ML"
        # add tree to physeq
            #write code if we decide we need a tree
    }


#6. Clean data for use in VIP project
    #prune samples with too few reads
        ps<-psVIP %>%
            prune_samples(sample_sums(.) > 50000, .) %>% # no samples with fewer than 50,000 reads
            prune_taxa(taxa_sums(.)>0, .) # remove any ESVs that now don't appear in any samples
    #select relevant columns from metadata
        newmetadata<-sample_data(ps) %>%
            as_tibble %>%
            select(ARMS, Fraction, Site_Depth_ft, Year, Site, Site_Latitude, Site_Longitude) %>%
            sample_data
        rownames(newmetadata)<-sample_names(ps)
        sample_data(ps)<-newmetadata

#5. Save data
    saveRDS(ps, file="../VIP_Project/Outputs/ps.RDS")
