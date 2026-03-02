#loading libraries
library(tidyverse)
library(phyloseq)
#library(ape)
#library(DECIPHER)
library(microbiome)
library(ggthemes) # additional themes fro ggplot2
library(ggpubr)
library(vegan)
library(repr)
library(ggpmisc) #to use stat_poly_eq
library(RColorBrewer) # nice color options
library(gridExtra) # gridding plots
library(viridis)
library(ggrepel)
library(wesanderson) #new palettes http://www.sthda.com/english/wiki/ggplot2-colors-how-to-change-colors-automatically-and-manually
library(rioja) # plotting poackages for tabular bubbleplots
library(reshape2) 
library(patchwork) #To construct composite plots
options(repr.plot.width=16, repr.plot.height=12)
set.seed(10000)

theme_glab <- function(base_size = 30,
                    base_family = "",
                    base_line_size = base_size / 180,
                    base_rect_size = base_size / 180) {
   
    font <- "Helvetica" #assign font family up front
   
    theme_bw(base_size = base_size,
                base_family = base_family,
                base_line_size = base_line_size) %+replace%
    theme(
        legend.background =  element_blank(),
        legend.title =       element_text(color = rgb(100, 100, 100, maxColorValue = 255),
                                          size = rel(0.65),
                                         hjust = 0),
        legend.text =        element_text(color = rgb(100, 100, 100, maxColorValue = 255),
                                          size = rel(0.65)),
        legend.key.size =    unit(0.8, "lines"),
     
      plot.title = element_text(
        color = rgb(100, 100, 100, maxColorValue = 255),
        hjust = 0),
       
      axis.title = element_text(
        color = rgb(100, 100, 100, maxColorValue = 255),
        size = rel(0.65)),
      axis.text = element_text(
        color = rgb(100, 100, 100, maxColorValue = 255),
        size = rel(0.65)),
       
      plot.caption = element_text(
        color = rgb(100, 100, 100, maxColorValue = 255),
        size = rel(0.7),
        hjust = 1),
       
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(),  
      panel.border = element_rect(fill = NA, colour = rgb(100, 100, 100, maxColorValue = 255)),

     
      complete = TRUE
    )
}

#uploading he env_dataset
vl_prj_env_data_full<-read.csv("../dataset/gracil_dataset_full.csv",heade=T, sep=',', row.names = 1)
head(vl_prj_env_data_full)

#upload 16S rDNA data for phyloseq
dna_prok_data<-phyloseq(sample_data(vl_prj_env_data_full),
                        otu_table(dna_otu_table,taxa_are_rows = F),
                        tax_table(dna_tax_table))
dna_prok_data

#upload 16S rDNA data for phyloseq
rna_otutable<-readRDS("../rds/rna/rna_otu_table.rds")
rna_taxtable<-readRDS("../rds/rna/rna_tax_table.rds")

rna_prok_data<-phyloseq(sample_data(vl_prj_env_data_full),
                        otu_table(rna_otutable,taxa_are_rows = F),
                        tax_table(rna_taxtable))
rna_prok_data

vl_prj_prok_data_raw_all<-merge_phyloseq(dna_prok_data,rna_prok_data)
vl_prj_prok_data_raw_all<-phyloseq(sample_data(vl_prj_env_data_full),
                                  otu_table(otu_table(vl_prj_prok_data_raw_all),taxa_are_rows = F),
                                  tax_table(tax_table(vl_prj_prok_data_raw_all))
                                  )
vl_prj_prok_data_raw_all

# Cleaning up unwanted sequences from Eukarya, mitochrondria and chloroplast

prok_data_euk_all <- subset_taxa(vl_prj_prok_data_raw_all,  (Kingdom != "Eukaryota") | is.na(Kingdom))
prok_data_euk_all
sum(readcount(prok_data_euk_all))
1 - (sum(readcount(prok_data_euk_all))/sum(readcount(vl_prj_prok_data_raw_all))) # sanity check percentage of reads after Eukaryota removal

prok_data_chl_all <- subset_taxa(prok_data_euk_all, (Order!="Chloroplast") | is.na(Order))
prok_data_chl_all
sum(readcount(prok_data_chl_all))
1- (sum(readcount(prok_data_chl_all))/sum(readcount(vl_prj_prok_data_raw_all))) # sanity check percentage of reads after Chloroplast removal

prok_data_mit_all <- subset_taxa(prok_data_chl_all, (Family!="Mitochondria") | is.na(Family))
prok_data_mit_all
sum(readcount(prok_data_mit_all))
1- (sum(readcount(prok_data_mit_all))/sum(readcount(vl_prj_prok_data_raw_all))) # sanity check percentage of reads after Chloroplast removal

## Removing the known DNA Extraction contaminants from Sheik et al., 2018
# Assuming you are starting from a phyloseq object called prok_data_raw
prok_data_contam <- subset_taxa(prok_data_mit_all,  (Genus != "Afipia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Aquabacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Asticcacaulis") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Aurantimonas") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Beijerinckia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Bosea") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Bradyrhizobium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Brevundimonas") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Caulobacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Craurococcus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Devosia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Hoefleae") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Mesorhizobium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Methylobacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Novosphingobium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Ochrobactrum") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Paracoccus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Pedomicrobium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Phyllobacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Rhizobium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Sphingobium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Sphingomonas") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Sphingopyxis") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Acidovorax") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Azoarcus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Azospira") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Burkholderia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Comamonas") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Cupriavidus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Curvibacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Delftiae") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Duganella") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Herbaspirillum") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Janthinobacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Kingella") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Leptothrix") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Limnobacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Massilia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Methylophilus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Methyloversatilis") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Oxalobacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Pelomonas") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Polaromonas") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Neisseria") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Ralstonia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Schlegelella") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Sulfuritalea") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Undibacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Variovorax") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Acinetobactera") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Enhydrobacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Enterobacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Escherichia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Nevskia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Pasteurella") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Pseudoxanthomonas") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Psychrobacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Stenotrophomonas") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Xanthomonas") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Aeromicrobium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Actinomyces") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Arthrobacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Beutenbergia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Brevibacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Corynebacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Curtobacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Dietzia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Janibacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Kocuria") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Microbacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Micrococcus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Microlunatus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Patulibacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Propionibacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Rhodococcus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Tsukamurella") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Chryseobacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Dyadobacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Flavobacterium") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Hydrotalea") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Niastella") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Parabacteroides") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Pedobacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Prevotella") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Wautersiella") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Deinococcus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Abiotrophia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Bacillus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Brevibacillus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Brochothrix") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Facklamia") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Olivibacter") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Lactobacillus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Paenibacillus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Ruminococcus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Staphylococcus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Streptococcus") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Veillonella") | is.na(Genus))
prok_data_contam <- subset_taxa(prok_data_contam,  (Genus != "Fusobacterium") | is.na(Genus))
prok_data_contam
message("reads saved from prok_data_raw after contaminants removal")
sum(readcount(prok_data_contam))/sum(readcount(vl_prj_prok_data_raw_all)) # reads saved from prok_data_raw after contaminants removal 
message("reads saved from prok_data_prune_mit after Eukaryota, Chloroplast, and Mitochondria removal")
sum(readcount(prok_data_contam))/sum(readcount(prok_data_mit_all)) # reads saved from prok_data_prune_mit after Eukaryota, Chloroplast, and Mitochondria removal

# Removing the potential human pathogens and contaminants
# This step needs to be evaluated with attention since many of these genera
# might be relevant in many environmental settings. Usually it is better to compare
# before-after removal to see what and how much you are removing. Feel free to
# experiment with the different groups and evaluate the results.

prok_data_humpath <- subset_taxa(prok_data_contam, (Genus != "Abiotrophia") |  is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Achromobacter") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Acinetobacter") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Actinobacillus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Arcanobacterium") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Babesia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Bifidobacterium") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Bartonella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Bordetella") |  is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Borrelia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Brodetella") |  is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Brucella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Capnocytophaga") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Chlamydia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Citrobacter") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Comamonas") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Corynebacterium_1") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Corynebacterium") | is.na(Genus)) # anche marino
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Coxiella") | is.na(Genus)) # anche marino
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Cronobacter") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Cutibacterium") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Dermatophilus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Ehrlichia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Order != "Enterobacteriales") | is.na(Order))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Enterococcus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Erysipelothrix") |  is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Escherichia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Escherichia/Shigella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Francisella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Gardnerella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Granulicatella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Haemophilus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Hafnia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Helicobacter") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Klebsiella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Kocuria") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Lactococcus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Lactobacillus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Lawsonia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Legionella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Leptospira") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Listeria") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Merkel_cell") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Micrococcus") | is.na(Genus)) # anche marino
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Morganella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Mycoplasma") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Neisseria") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Nocardia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Pasteurella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Plesiomonas") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Propionibacterium") | is.na(Genus))
#prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Proteus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Providencia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Pseudomonas") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Rhodococcus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Rickettsiae") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Roseomonas") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Rothia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Salmonella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Serratia") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Shewanella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Shigella") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Sphaerophorus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Staphylococcus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Stenotrophomonas") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Streptococcus") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Treponema") | is.na(Genus))
#prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Vibrio") | is.na(Genus))
prok_data_humpath <- subset_taxa(prok_data_humpath, (Genus != "Yersinia") | is.na(Genus))
message("reads saved from prok_data_raw after contaminants removal")
sum(readcount(prok_data_humpath))/sum(readcount(vl_prj_prok_data_raw_all)) # reads saved from prok_data_raw after contaminants removal 
message("reads saved from prok_data_contam after DNA KIT EXTRAC. contaminants removal")
sum(readcount(prok_data_humpath))/sum(readcount(prok_data_contam)) # reads saved from prok_data_contam after DNA KIT EXTRAC. contaminants removal
prok_data_humpath
sum(readcount(prok_data_humpath))

cleaning_table1<-t(data.frame(
c((1 - (sum(readcount(prok_data_euk_all))/sum(readcount(vl_prj_prok_data_raw_all))))*100,'%'),
c((1 - (sum(readcount(prok_data_chl_all))/sum(readcount(vl_prj_prok_data_raw_all))))*100,'%'),
c((1 - (sum(readcount(prok_data_mit_all))/sum(readcount(vl_prj_prok_data_raw_all))))*100,'%'),
c((1 - (sum(readcount(prok_data_contam))/sum(readcount(prok_data_mit_all))))*100,'%'),
c((1 - (sum(readcount(prok_data_humpath))/sum(readcount(prok_data_contam))))*100,'%')    
    )
 )
rownames(cleaning_table1)<-c("Eukaryota","Chloroplast","Mitochondria","Kit contam.","Human contam.")
cleaning_table1 #tiny organization of the contamination results

#sanity check
prok_data_all2 = filter_taxa(prok_data_humpath, function(x) sum(x) > 0, TRUE)
prok_data_all2
sum(readcount(prok_data_all2))
100*(1-(sum(readcount(prok_data_all2))/sum(readcount(vl_prj_prok_data_raw_all))))

#preparing for figure S2

#figure S2_A
#svg("rarefaction_curve_dna.svg", width=12,height=12)
rarecurve(as.matrix(data.frame(otu_table(subset_samples(vl_prj_prok_data_raw_all,type=='dna')))), step=100 , lwd=3, ylab="ASVs", label=T)
#dev.off() #using INKSCAPE to adjust the colors and other details

#figure S2_B
#svg("rarefaction_curve_dna.svg", width=12,height=12)
rarecurve(as.matrix(data.frame(otu_table(subset_samples(vl_prj_prok_data_raw_all,type=='rna')))), step=100 , lwd=3, ylab="ASVs", label=T)
#dev.off() #using INKSCAPE to adjust the colors and other details

save.image()

#preparing for figure S3

#figure S3
alpha_div2<-plot_richness(subset_samples(vl_prj_prok_data_raw_all, time != "T0"),
              measures = c("Observed", "Shannon"), x = "substrate3") +
  geom_point(aes(
    fill = ifelse(type == "rna", "white", substrate2),
    shape = type
                  ), size = 5) +
  scale_fill_manual("Sample type",
                    values = c("Bare" = "#dccd7dc2", "Macroalgae" = "#7e2954c0",
                               "Seagrass" = "#337538ff", "white" = "white"),
                    labels = c("Bare", "Macroalgae", "Seagrass")) +
  scale_shape_manual(values = c("dna" = 21, "rna" = 24)) +
  scale_x_discrete(labels = c("T1", "T1C", "T2", "T2C",
                              "T1", "T1C", "T2", "T2C",
                              "T1", "T1C", "T2", "T2C")) +
  xlab("") +
  theme_glab() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 0.5),
        axis.text = element_text(size = 20),
        legend.position = "right",
        aspect.ratio=0.7)

alpha_div2$data$substrate3 <- factor(alpha_div2$data$substrate3, levels = c("Bare_T1", "Bare_1C", "Bare_T2", "Bare_2C",
                                                    "Macroalgae_T1", "Macroalgae_1C", "Macroalgae_T2", "Macroalgae_2C",
                                                    "Seagrass_T1", "Seagrass_1C", "Seagrass_T2", "Seagrass_2C"))

alpha_div2

#preparing for figure 2

#figure 2B
plot_richness(subset_samples(vl_prj_prok_data_raw_all, group != "control"),
              measures = c("Observed", "Shannon"), x = "substrate3") +
  geom_point(aes(
    fill = ifelse(type == "rna", "white", substrate2),
    shape = type
                  ), size = 5) +
  scale_fill_manual("Sample type",
                    values = c("Bare" = "#dccd7dc2", "Macroalgae" = "#7e2954c0",
                               "Seagrass" = "#337538ff", "white" = "white"),
                    labels = c("Bare", "Macroalgae", "Seagrass")) +
  scale_shape_manual(values = c("dna" = 21, "rna" = 24)) +
  scale_x_discrete(labels = c("T0", "T1", "T2", "T0", "T1", "T2", "T0", "T1", "T2")) +
  xlab("") +
  theme_glab() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.25, hjust = 0.5),
        axis.text = element_text(size = 20),
        legend.position = "right",
        aspect.ratio=0.7)

library(rempsyc)
library(flextable)
library(broom)

alpha_diversity_table <-
                        cbind(
                        readcount(vl_prj_prok_data_raw_all),
                        estimate_richness(vl_prj_prok_data_raw_all, split = TRUE, measures = "Observed"),
                        estimate_richness(vl_prj_prok_data_raw_all, split = TRUE, measures = "Shannon"),
                        estimate_richness(vl_prj_prok_data_raw_all, split = TRUE, measures = "Chao1"),
                        estimate_richness(vl_prj_prok_data_raw_all, split = TRUE, measures = "Observed"))

alpha_diversity_table<-alpha_diversity_table[,-4]
colnames(alpha_diversity_table) <- c("Reads","Observed", "Shannon", "Chao1", "ASV")
alpha_diversity_table

save.image()

#preparing for figure 2

#Uploading average dataset
env_data_avg<-read.csv("../dataset/gracil_df.csv",header=T,sep=',',row.names = 1)
head(env_data_avg)

#preparing otu_table for average calc transformation

otu_dna<-data.frame(otu_table(subset_samples(prok_data_all2,type=="dna")))
head(otu_dna)

dna_sample_names<-c("Bare_T0","Bare_T1","Bare_T1","Bare_T1","Seagrass_T1","Seagrass_T1","Seagrass_T1","Macroalgae_T1",
  "Macroalgae_T1","Macroalgae_T1","Bare_1C","Bare_T0","Bare_1C","Bare_1C","Seagrass_1C","Seagrass_1C",
  "Seagrass_1C","Macroalgae_1C","Macroalgae_1C","Macroalgae_1C","Bare_T2","Bare_T2","Bare_T0","Bare_T2",
  "Seagrass_T2","Seagrass_T2","Seagrass_T2","Macroalgae_T2","Macroalgae_T2","Macroalgae_T2","Bare_2C",
  "Bare_2C","Bare_2C","Macroalgae_T0","Macroalgae_2C","Macroalgae_2C","Macroalgae_2C","Seagrass_2C",
  "Seagrass_2C","Seagrass_2C","Macroalgae_T0","Macroalgae_T0","Seagrass_T0","Seagrass_T0","Seagrass_T0")

# moving column replics in the first position
otu_dna2<-cbind(otu_dna,dna_sample_names)
otu_dna2<-otu_dna2 %>% select(dna_sample_names, everything())
head(otu_dna2)

# calc the average
otu_dna3 <- aggregate(otu_dna2[,2:58788], by=list(otu_dna2$dna_sample_names), FUN=mean, na.rm=F)
names(otu_dna3)[names(otu_dna3) == 'Group.1'] <- 'name'
otu_dna3

#naming the rows by the replica column
row.names(otu_dna3)<-otu_dna3[,1]
otu_dna3<-otu_dna3[,-1]
head(otu_dna3)

#phyloseq object with mean values of replicas
prok_data_dna_av<-phyloseq(sample_data(env_data_avg),
                           otu_table(otu_dna3,taxa_are_rows=F),
                           tax_table(tax_table(dna_prok_data)))
prok_data_dna_av

prok_ndata_dna_avg <- transform_sample_counts(prok_data_dna_av, function(x) ((x / sum(x))*median(readcount(prok_data_dna_av))))

prok_ra_dna_avg = transform_sample_counts(prok_ndata_dna_avg, function(x){x / sum(x)})

## Agglomerate at a specific taxonomic level at the Genus level
prok_ra_dna_genus_avg = tax_glom(prok_ra_dna_avg, "Genus", NArm = FALSE)
prok_ra_dna_family_avg = tax_glom(prok_ra_dna_avg, "Family", NArm = FALSE)
prok_ra_dna_order_avg = tax_glom(prok_ra_dna_avg, "Order", NArm = FALSE)
prok_ra_dna_class_avg = tax_glom(prok_ra_dna_avg, "Class", NArm = FALSE)
prok_ra_dna_phyla_avg = tax_glom(prok_ra_dna_avg, "Phylum", NArm = FALSE)
prok_ra_dna_kingdom_avg = tax_glom(prok_ra_dna_avg, "Kingdom", NArm = FALSE)

library(MicEco)

#Barplot Class level: DNA

prok_ra_dna_class_avg_noCTRL<-subset_samples(prok_ra_dna_class_avg,group != 'control')
prok_ra_dna_class_avg_noCTRL = transform_sample_counts(prok_ra_dna_class_avg_noCTRL, function(x){x / sum(x)})

prok_ra_dna_class_avg_noCTRL_perc1 <- ps_prune(prok_ra_dna_class_avg_noCTRL, min.abundance = 0.01)
#prok_ra_dna_class_avg_noCTRL_perc1 <-ps_tax_clean(prok_ra_dna_class_avg_noCTRL_perc1)

#figure 2A
plot_bar(prok_ra_dna_class_avg_noCTRL_perc1, fill ="Class") +
  ggtitle("Class Abundance of >1% ASVs") +
  theme(axis.text.x = element_text(face="bold",size=10))+
  theme_glab() + theme(legend.position = "bottom") +
  scale_fill_manual(na.value="black","Class",values = colorRampPalette(brewer.pal(11, "Spectral"))(15)) + 
  scale_x_discrete(limits=rev) +
  coord_flip()

#Class barplot for RNA samples

otu_rna<-(otu_table(subset_samples(prok_data_all2,type=="rna")))
taxa_rna<-(tax_table(subset_samples(prok_data_all2,type=="rna")))

prok_data_rna2<-phyloseq(sample_data(subset_samples(prok_data_all2,type=="rna")),
                         otu_table(otu_rna,taxa_are_rows = F),
                         tax_table(taxa_rna))
prok_data_rna2

prok_ndata_rna <- transform_sample_counts(prok_data_rna2, function(x) ((x / sum(x))*median(readcount(prok_data_rna2))))

prok_ra_rna = transform_sample_counts(prok_ndata_rna, function(x){x / sum(x)})

## Agglomerate at a specific taxonomic level at the Genus level
prok_ra_rna_genus = tax_glom(prok_ra_rna, "Genus", NArm = FALSE)
prok_ra_rna_family = tax_glom(prok_ra_rna, "Family", NArm = FALSE)
prok_ra_rna_order = tax_glom(prok_ra_rna, "Order", NArm = FALSE)
prok_ra_rna_class = tax_glom(prok_ra_rna, "Class", NArm = FALSE)
prok_ra_rna_phyla = tax_glom(prok_ra_rna, "Phylum", NArm = FALSE)
prok_ra_rna_kingdom = tax_glom(prok_ra_rna, "Kingdom", NArm = FALSE)

prok_ra_rna_class_noCTRL<-subset_samples(prok_ra_rna_class,group != 'control')
prok_ra_rna_class_noCTRL = transform_sample_counts(prok_ra_rna_class_noCTRL, function(x){x / sum(x)})

prok_ra_rna_class_noCTRL_perc1 <- ps_prune(prok_ra_rna_class_noCTRL, min.abundance = 0.01)
#prok_ra_rna_class_noCTRL_perc1 <-ps_tax_clean(prok_ra_rna_class_noCTRL_perc1)

plot_bar(prok_ra_rna_class_noCTRL_perc1, fill ="Class") +
  ggtitle("Class Abundance of >1% ASVs") +
  theme(axis.text.x = element_text(face="bold",size=10))+
  theme_glab() + theme(legend.position = "bottom") +
  scale_fill_manual(na.value="black","Class",values = colorRampPalette(brewer.pal(11, "Spectral"))(13)) + 
  scale_x_discrete(limits=rev) +
  coord_flip()

#use inkscape to give the get the same colors for each class in dna and rna samples

env_data_avg

library(svglite)

ggarrange(
ggplot(subset(env_data_avg,group=="sample"), aes(x=type,y=mca)) + 
geom_bar(stat = "identity",fill="#440154",linewidth=0.446,color="black",) +
scale_x_discrete(limits=rev) + ylab("uM") + xlab("") + ylab(expression(nmol~g^{-1}~h^{-1})) + xlab("") +
#coord_flip() +
theme_glab() + theme(axis.text = element_text(size = 8),
                     axis.text.x = element_text(angle = 90, vjust = 0.5),
                           legend.position = "none", aspect.ratio=0.7) ,
    
ggplot(subset(env_data_avg,group=="sample"), aes(x=type,y=dic.conc)) + 
geom_bar(stat = "identity",fill="#3b528b",linewidth=0.446,color="black") +
scale_x_discrete(limits=rev) + ylab("uM") + xlab("") +
#coord_flip() +
theme_glab() + theme(axis.text = element_text(size = 8),
                     axis.text.x = element_text(angle = 90, vjust = 0.5),
                           legend.position = "none", aspect.ratio=0.7, axis.text.y=element_blank(),
                    axis.ticks.y=element_blank()) ,
    
ggplot(subset(env_data_avg,group=="sample"), aes(x=type,y=no3.conc)) + 
geom_bar(stat = "identity",fill="#21918c",linewidth=0.446,color="black") +
scale_x_discrete(limits=rev) + ylab("uM") + xlab("") +
#coord_flip() +
theme_glab() + theme(axis.text = element_text(size = 8),
                     axis.text.x = element_text(angle = 90, vjust = 0.5),
                           legend.position = "none", aspect.ratio=0.7, axis.text.y=element_blank(),
                    axis.ticks.y=element_blank()) ,
    
ggplot(subset(env_data_avg,group=="sample"), aes(x=type,y=nh4.conc)) + 
geom_bar(stat = "identity",fill="#5ec962",linewidth=0.446,color="black") +
scale_x_discrete(limits=rev) + ylab("uM") + xlab("") +
#coord_flip() +
theme_glab() + theme(axis.text = element_text(size = 8),
                     axis.text.x = element_text(angle = 90, vjust = 0.5),
                           legend.position = "none", aspect.ratio=0.7, axis.text.y=element_blank(),
                    axis.ticks.y=element_blank()) ,
    
ggplot(subset(env_data_avg,group=="sample"), aes(x=type,y=o2.conc)) + 
geom_bar(stat = "identity",fill="#fde725",linewidth=0.446,color="black") +
scale_x_discrete(limits=rev) + ylab("uM") + xlab("") +
#coord_flip() +
theme_glab() + theme(axis.text = element_text(size = 8),
                     axis.text.x = element_text(angle = 90, vjust = 0.5),
                           legend.position = "none", aspect.ratio=0.7, axis.text.y=element_blank(),
                    axis.ticks.y=element_blank()), 
ncol=1,nrow=5, align = "hv")
#ggsave("~/Downloads/prova.svg",height=18,width=12)

save.image()

#preparing for figure 3

prok_ndata_dna_noCTRL <- transform_sample_counts(subset_samples(subset_samples(prok_ndata_all,type=='dna'),group=='sample'),
                                          function(x) ((x / sum(x))*median(readcount(subset_samples(subset_samples(prok_ndata_all,type=='dna'),group=='sample')))))

#calc distance matrix: Weighted Jaccard index
prok_dist_uwjac_dna <- distance(prok_ndata_dna_noCTRL, method = "jaccard",binary=TRUE)

# NMDS Jaccard Weighted
prok_nmds_juw_dna <- ordinate(prok_ndata_dna_noCTRL,prok_dist_uwjac_dna, method = "NMDS", trymax=9999)

p_nmds_uwj_dna<-(plot_ordination(prok_ndata_dna_noCTRL, prok_nmds_juw_dna,type="samples"))$data

# build ggplot dataframe with points (x,y) and corresponding groups (cluster)
nmds_uwj_dna_polygon <- data.frame(cluster=factor(p_nmds_uwj_dna$time), 
                                  x=p_nmds_uwj_dna$NMDS1, y=p_nmds_uwj_dna$NMDS2)


# calculate group centroid locations
centroids_dna <- aggregate(cbind(x,y)~cluster,data=nmds_uwj_dna_polygon,mean)


# merge centroid locations into ggplot dataframe
nmds_uwj_dna_polygon2 <- merge(nmds_uwj_dna_polygon,centroids_dna,by="cluster",suffixes=c("",".centroid"))

#figure 3
ggplot(p_nmds_uwj_dna, aes(x=NMDS1, y=NMDS2)) +
geom_point(aes(fill=substrate2,color=substrate2,shape=time),size=8,stroke=0.3) +
#geom_text(aes(label=name,),color="black", size=4,hjust=0.5, vjust=-0.6) + 
scale_fill_manual("Substrate",values = c("Bare" = "#dccd7dc2", "Macroalgae" = "#7e2954c0",
                               "Seagrass" = "#337538ff")) +
scale_color_manual(values = c("Bare" = "#dccd7dc2", "Macroalgae" = "#7e2954c0",
                               "Seagrass" = "#337538ff")) +
scale_shape_manual("Time",values = c("T0" = 21, "T1" = 22, "T2" = 23))  +
stat_chull(aes(color = substrate3, fill = substrate2),alpha = 0.1, geom = "polygon") + 
 

geom_point(data=centroids_dna, aes(x=x, y=y),fill="black",shape=21, size=4) +

geom_segment(data=nmds_uwj_dna_polygon2,aes(x=x.centroid, y=y.centroid, xend=x, yend=y),
             linetype = "dashed",color="black",alpha=.5) +

stat_conf_ellipse(aes(color=time)) +


theme_glab() + theme(legend.position = "right") + 
guides(fill = guide_legend(override.aes = list(shape = 21)))

#ggsave("plot/nmds_jw_polygon.svg",width=12,height=8)

adonis2(prok_dist_uwjac_dna~time,permutations=999,
    data = p_nmds_uwj_dna)

adonis2(prok_dist_uwjac_dna~substrate,permutations=999,
    data = p_nmds_uwj_dna)

library(pairwiseAdonis)

pairwise.adonis(x = prok_dist_uwjac_dna, 
                          groups = time, 
                          p.adjust.method = "bonferroni", 
                          perm = 9999)



save.image()

#Preparing Figure 4

prok_ndata_dna <- transform_sample_counts(subset_samples(prok_ndata_all,type=='dna'),
                  function(x) ((x / sum(x))*median(readcount(subset_samples(prok_ndata_all,type=='dna')))))

prok_ndata_bare_dna<-subset_samples(prok_ndata_dna, substrate2== "Bare")
prok_ndata_seagr_dna<-subset_samples(prok_ndata_dna, substrate2== "Seagrass")
prok_ndata_macroalg_dna<-subset_samples(prok_ndata_dna, substrate2== "Macroalgae")

prok_wjac_bare <- distance(prok_ndata_bare_dna, method = "jaccard")
prok_wjac_macroalg <- distance(prok_ndata_macroalg_dna, method = "jaccard")
prok_wjac_seagr <- distance(prok_ndata_seagr_dna, method = "jaccard")

prok_uwjac_bare <- distance(prok_ndata_bare_dna, method = "jaccard", binary = TRUE)
prok_uwjac_macroalg <- distance(prok_ndata_macroalg_dna, method = "jaccard", binary = TRUE)
prok_uwjac_seagr <- distance(prok_ndata_seagr_dna, method = "jaccard", binary = TRUE)

# NMDS FOR BARE SUBSTRATE

# NMDS Jaccard Weighted - BARE SUBSTRATE
prok_nmds_jw_bare_dna <- ordinate(prok_ndata_bare_dna,prok_wjac_bare, method = "NMDS", weighted=T, trymax=9999)

p_nmds_wj_bare_dna<-(plot_ordination(prok_ndata_bare_dna, prok_nmds_jw_bare_dna,type="samples"))$data

# build ggplot dataframe with points (x,y) and corresponding groups (cluster)
p_nmds_wj_bare_dna_polygon <- data.frame(cluster=factor(p_nmds_wj_bare_dna$substrate3), 
                                  x=p_nmds_wj_bare_dna$NMDS1, y=p_nmds_wj_bare_dna$NMDS2)


# calculate group centroid locations
centroids_bare_dna <- aggregate(cbind(x,y)~cluster,data=p_nmds_wj_bare_dna_polygon,mean)


# merge centroid locations into ggplot dataframe
p_nmds_wj_bare_dna_polygon2 <- merge(p_nmds_wj_bare_dna_polygon,centroids_bare_dna,by="cluster",suffixes=c("",".centroid"))

fig_s4a<-ggplot(p_nmds_wj_bare_dna, aes(x=NMDS1, y=NMDS2)) +
geom_point(aes(fill=substrate2,color=substrate2,shape=substrate3),size=8,stroke=0.3) +
#geom_text(aes(label=name,),color="black", size=4,hjust=0.5, vjust=-0.6) + 
scale_fill_manual("Substrate",values = c("Bare" = "#dccd7dc2")) +
scale_color_manual(values = c("Bare" = "#dccd7dc2")) +
scale_shape_manual("Time",values = c(24,25,21,22,23))  +
stat_chull(aes(color = substrate3, fill = substrate2),alpha = 0.1, geom = "polygon") + 
 
geom_point(data=centroids_bare_dna, aes(x=x, y=y),fill="black",shape=21, size=4) +
geom_segment(data=p_nmds_wj_bare_dna_polygon2,aes(x=x.centroid, y=y.centroid, xend=x, yend=y),
             linetype = "dashed",color="black",alpha=.5) +

theme_glab() + theme(legend.position = "right") + 
guides(fill = guide_legend(override.aes = list(shape = 21)))
fig_s4a
#ggsave("plot/nmds_jw_polygon.svg",width=12,height=8)

adonis2(prok_wjac_bare~time,permutations=999,
    data = p_nmds_wj_bare_dna)

# NMDS Jaccard Unweighted - BARE SUBSTRATE
prok_nmds_juw_bare_dna <- ordinate(prok_ndata_bare_dna,prok_uwjac_bare, method = "NMDS", trymax=9999)

p_nmds_uwj_bare_dna<-(plot_ordination(prok_ndata_bare_dna, prok_nmds_juw_bare_dna,type="samples"))$data

# build ggplot dataframe with points (x,y) and corresponding groups (cluster)
p_nmds_uwj_bare_dna_polygon <- data.frame(cluster=factor(p_nmds_uwj_bare_dna$substrate3), 
                                  x=p_nmds_uwj_bare_dna$NMDS1, y=p_nmds_uwj_bare_dna$NMDS2)


# calculate group centroid locations
centroids_bare_dna_uwj <- aggregate(cbind(x,y)~cluster,data=p_nmds_uwj_bare_dna_polygon,mean)


# merge centroid locations into ggplot dataframe
p_nmds_uwj_bare_dna_polygon2 <- merge(p_nmds_uwj_bare_dna_polygon,centroids_bare_dna_uwj,by="cluster",suffixes=c("",".centroid"))

fig_s4b<-ggplot(p_nmds_uwj_bare_dna, aes(x=NMDS1, y=NMDS2)) +
geom_point(aes(fill=substrate2,color=substrate2,shape=substrate3),size=8,stroke=0.3) +
#geom_text(aes(label=name,),color="black", size=4,hjust=0.5, vjust=-0.6) + 
scale_fill_manual("Substrate",values = c("Bare" = "#dccd7dc2")) +
scale_color_manual(values = c("Bare" = "#dccd7dc2")) +
scale_shape_manual("Time",values = c(24,25,21,22,23))  +
stat_chull(aes(color = substrate3, fill = substrate2),alpha = 0.1, geom = "polygon") + 
 
geom_point(data=centroids_bare_dna_uwj, aes(x=x, y=y),fill="black",shape=21, size=4) +
geom_segment(data=p_nmds_uwj_bare_dna_polygon2,aes(x=x.centroid, y=y.centroid, xend=x, yend=y),
             linetype = "dashed",color="black",alpha=.5) +

theme_glab() + theme(legend.position = "right") + 
guides(fill = guide_legend(override.aes = list(shape = 21)))
fig_s4b
#ggsave("plot/nmds_jw_polygon.svg",width=12,height=8)

adonis2(prok_uwjac_bare~time,permutations=999,
    data = p_nmds_uwj_bare_dna)

# NMDS FOR MACROALGAE SUBSTRATE

# NMDS Jaccard Weighted - MACROALGAE SUBSTRATE
prok_nmds_jw_macroalg_dna <- ordinate(prok_ndata_macroalg_dna,prok_wjac_macroalg, method = "NMDS", weighted=T, trymax=9999)

p_nmds_wj_macroalg_dna<-(plot_ordination(prok_ndata_macroalg_dna, prok_nmds_jw_macroalg_dna,type="samples"))$data

# build ggplot dataframe with points (x,y) and corresponding groups (cluster)
p_nmds_wj_macroalg_dna_polygon <- data.frame(cluster=factor(p_nmds_wj_macroalg_dna$substrate3), 
                                  x=p_nmds_wj_macroalg_dna$NMDS1, y=p_nmds_wj_macroalg_dna$NMDS2)


# calculate group centroid locations
centroids_macroalg_dna <- aggregate(cbind(x,y)~cluster,data=p_nmds_wj_macroalg_dna_polygon,mean)


# merge centroid locations into ggplot dataframe
p_nmds_wj_macroalg_dna_polygon2 <- merge(p_nmds_wj_macroalg_dna_polygon,centroids_macroalg_dna,by="cluster",suffixes=c("",".centroid"))

fig_s4c<-ggplot(p_nmds_wj_macroalg_dna, aes(x=NMDS1, y=NMDS2)) +
geom_point(aes(fill=substrate2,color=substrate2,shape=substrate3),size=8,stroke=0.3) +
#geom_text(aes(label=name,),color="black", size=4,hjust=0.5, vjust=-0.6) + 
scale_fill_manual("Substrate",values = c("Macroalgae" = "#7e2954c0")) +
scale_color_manual(values = c("Macroalgae" = "#dccd7dc2")) +
scale_shape_manual("Time",values = c(24,25,21,22,23))  +
stat_chull(aes(color = substrate3, fill = substrate2),alpha = 0.1, geom = "polygon") + 
 
geom_point(data=centroids_macroalg_dna, aes(x=x, y=y),fill="black",shape=21, size=4) +
geom_segment(data=p_nmds_wj_macroalg_dna_polygon2,aes(x=x.centroid, y=y.centroid, xend=x, yend=y),
             linetype = "dashed",color="black",alpha=.5) +

theme_glab() + theme(legend.position = "right") + 
guides(fill = guide_legend(override.aes = list(shape = 21)))
fig_s4c
#ggsave("plot/nmds_jw_polygon.svg",width=12,height=8)

adonis2(prok_wjac_macroalg~time,permutations=999,
    data = p_nmds_wj_macroalg_dna)

# NMDS Jaccard Unweighted - macroalg SUBSTRATE
prok_nmds_juw_macroalg_dna <- ordinate(prok_ndata_macroalg_dna,prok_uwjac_macroalg, method = "NMDS", trymax=9999)

p_nmds_uwj_macroalg_dna<-(plot_ordination(prok_ndata_macroalg_dna, prok_nmds_juw_macroalg_dna,type="samples"))$data

# build ggplot dataframe with points (x,y) and corresponding groups (cluster)
p_nmds_uwj_macroalg_dna_polygon <- data.frame(cluster=factor(p_nmds_uwj_macroalg_dna$substrate3), 
                                  x=p_nmds_uwj_macroalg_dna$NMDS1, y=p_nmds_uwj_macroalg_dna$NMDS2)


# calculate group centroid locations
centroids_macroalg_dna_uwj <- aggregate(cbind(x,y)~cluster,data=p_nmds_uwj_macroalg_dna_polygon,mean)


# merge centroid locations into ggplot dataframe
p_nmds_uwj_macroalg_dna_polygon2 <- merge(p_nmds_uwj_macroalg_dna_polygon,centroids_macroalg_dna_uwj,by="cluster",suffixes=c("",".centroid"))

fig_s4d<-ggplot(p_nmds_uwj_macroalg_dna, aes(x=NMDS1, y=NMDS2)) +
geom_point(aes(fill=substrate2,color=substrate2,shape=substrate3),size=8,stroke=0.3) +
#geom_text(aes(label=name,),color="black", size=4,hjust=0.5, vjust=-0.6) + 
scale_fill_manual("Substrate",values = c("Macroalgae" = "#7e2954c0")) +
scale_color_manual(values = c("Macroalgae" = "#7e2954c0")) +
scale_shape_manual("Time",values = c(24,25,21,22,23))  +
stat_chull(aes(color = substrate3, fill = substrate2),alpha = 0.1, geom = "polygon") + 
 
geom_point(data=centroids_macroalg_dna_uwj, aes(x=x, y=y),fill="black",shape=21, size=4) +
geom_segment(data=p_nmds_uwj_macroalg_dna_polygon2,aes(x=x.centroid, y=y.centroid, xend=x, yend=y),
             linetype = "dashed",color="black",alpha=.5) +

theme_glab() + theme(legend.position = "right") + 
guides(fill = guide_legend(override.aes = list(shape = 21)))
fig_s4d
#ggsave("plot/nmds_jw_polygon.svg",width=12,height=8)

adonis2(prok_uwjac_bare~time,permutations=999,
    data = p_nmds_wj_bare_dna)

# NMDS FOR SEAGRASS SUBSTRATE

# NMDS Jaccard Weighted - MACROALGAE SUBSTRATE
prok_nmds_jw_seagr_dna <- ordinate(prok_ndata_seagr_dna,prok_wjac_seagr, method = "NMDS", weighted=T, trymax=9999)

p_nmds_wj_seagr_dna<-(plot_ordination(prok_ndata_seagr_dna, prok_nmds_jw_seagr_dna,type="samples"))$data

# build ggplot dataframe with points (x,y) and corresponding groups (cluster)
p_nmds_wj_seagr_dna_polygon <- data.frame(cluster=factor(p_nmds_wj_seagr_dna$substrate3), 
                                  x=p_nmds_wj_seagr_dna$NMDS1, y=p_nmds_wj_seagr_dna$NMDS2)


# calculate group centroid locations
centroids_seagr_dna <- aggregate(cbind(x,y)~cluster,data=p_nmds_wj_seagr_dna_polygon,mean)


# merge centroid locations into ggplot dataframe
p_nmds_wj_seagr_dna_polygon2 <- merge(p_nmds_wj_seagr_dna_polygon,centroids_seagr_dna,by="cluster",suffixes=c("",".centroid"))

fig_s4e<-ggplot(p_nmds_wj_seagr_dna, aes(x=NMDS1, y=NMDS2)) +
geom_point(aes(fill=substrate2,color=substrate2,shape=substrate3),size=8,stroke=0.3) +
#geom_text(aes(label=name,),color="black", size=4,hjust=0.5, vjust=-0.6) + 
scale_fill_manual("Substrate",values = c("Seagrass" = "#337538ff")) +
scale_color_manual(values = c("Seagrass" = "#337538ff")) +
scale_shape_manual("Time",values = c(24,25,21,22,23))  +
stat_chull(aes(color = substrate3, fill = substrate2),alpha = 0.1, geom = "polygon") + 
 
geom_point(data=centroids_seagr_dna, aes(x=x, y=y),fill="black",shape=21, size=4) +
geom_segment(data=p_nmds_wj_seagr_dna_polygon2,aes(x=x.centroid, y=y.centroid, xend=x, yend=y),
             linetype = "dashed",color="black",alpha=.5) +

theme_glab() + theme(legend.position = "right") + 
guides(fill = guide_legend(override.aes = list(shape = 21)))
fig_s4e
#ggsave("plot/nmds_jw_polygon.svg",width=12,height=8)

adonis2(prok_wjac_seagr~time,permutations=999,
    data = p_nmds_wj_seagr_dna)

# NMDS Jaccard Unweighted - seagr SUBSTRATE
prok_nmds_juw_seagr_dna <- ordinate(prok_ndata_seagr_dna,prok_uwjac_seagr, method = "NMDS", trymax=9999)

p_nmds_uwj_seagr_dna<-(plot_ordination(prok_ndata_seagr_dna, prok_nmds_juw_seagr_dna,type="samples"))$data

# build ggplot dataframe with points (x,y) and corresponding groups (cluster)
p_nmds_uwj_seagr_dna_polygon <- data.frame(cluster=factor(p_nmds_uwj_seagr_dna$substrate3), 
                                  x=p_nmds_uwj_seagr_dna$NMDS1, y=p_nmds_uwj_seagr_dna$NMDS2)


# calculate group centroid locations
centroids_seagr_dna_uwj <- aggregate(cbind(x,y)~cluster,data=p_nmds_uwj_seagr_dna_polygon,mean)


# merge centroid locations into ggplot dataframe
p_nmds_uwj_seagr_dna_polygon2 <- merge(p_nmds_uwj_seagr_dna_polygon,centroids_seagr_dna_uwj,by="cluster",suffixes=c("",".centroid"))

fig_s4f<-ggplot(p_nmds_uwj_seagr_dna, aes(x=NMDS1, y=NMDS2)) +
geom_point(aes(fill=substrate2,color=substrate2,shape=substrate3),size=8,stroke=0.3) +
#geom_text(aes(label=name,),color="black", size=4,hjust=0.5, vjust=-0.6) + 
scale_fill_manual("Substrate",values = c("Seagrass" = "#337538ff")) +
scale_color_manual(values = c("Seagrass" = "#337538ff")) +
scale_shape_manual("Time",values = c(24,25,21,22,23))  +
stat_chull(aes(color = substrate3, fill = substrate2),alpha = 0.1, geom = "polygon") + 
 
geom_point(data=centroids_seagr_dna_uwj, aes(x=x, y=y),fill="black",shape=21, size=4) +
geom_segment(data=p_nmds_uwj_seagr_dna_polygon2,aes(x=x.centroid, y=y.centroid, xend=x, yend=y),
             linetype = "dashed",color="black",alpha=.5) +

theme_glab() + theme(legend.position = "right") + 
guides(fill = guide_legend(override.aes = list(shape = 21)))
fig_s4f
#ggsave("plot/nmds_jw_polygon.svg",width=12,height=8)

adonis2(prok_uwjac_seagr~time,permutations=999,
    data = p_nmds_wj_seagr_dna)

save.image()

#Preparing figure S5

library(VennDiagram)
library(gplots)
library(venn)
library(UpSetR)
library(MicrobiotaProcess)
library(ggvenn)

# using prok_data_all2 that contain both dna and rna filtered samples: normalize by median and transform into relative abundance
prok_ndata_all2<-transform_sample_counts(prok_data_all2,
                  function(x) ((x / sum(x))*median(readcount(prok_data_all2))))

prok_ra_all2 = transform_sample_counts(prok_ndata_all2, function(x){x / sum(x)})                                                       

#DNA
#Bare

prok_ra_bare_dna_noCtrl<-subset_samples(subset_samples(subset_samples(subset_samples(prok_ra_all2, type %in% 'dna'), 
                                             group %in% 'sample'),
                                             rep %in% 'R1'),
                                             substrate2 %in% 'Bare')

prok_ra_merged_bare_dna <- merge_samples(prok_ra_bare_dna_noCtrl,"time")
prok_ra_merged_bare_dna

# remove controls, divide dna and rna samples and get the list for the venn diagrams
prok_ra_bare_dna_list<-get_vennlist(prok_ra_bare_dna_noCtrl,factorNames='time')

ggvenn(
  prok_ra_bare_dna_list, 
  fill_color = c("#dccd7dc2","#dccd7dc2","#dccd7dc2"),
  stroke_size = 0.5, set_name_size = 10, text_size=8
  ) + theme_glab()

#Macroalgae

prok_ra_macroalg_dna_noCtrl<-subset_samples(subset_samples(subset_samples(subset_samples(prok_ra_all2, type %in% 'dna'), 
                                             group %in% 'sample'),
                                             rep %in% 'R1'),
                                            substrate2 %in% 'Macroalgae')

# remove controls, divide dna and rna samples and get the list for the venn diagrams
prok_ra_macroalg_dna_list<-get_vennlist(prok_ra_macroalg_dna_noCtrl,factorNames='time')

ggvenn(
  prok_ra_macroalg_dna_list, 
  fill_color = c("#7e2954c0","#7e2954c0","#7e2954c0"),
  stroke_size = 0.5, set_name_size = 10, text_size=8
  ) + theme_glab()

prok_ra_merged_macroalg_dna <- merge_samples(prok_ra_macroalg_dna_noCtrl,"time")
prok_ra_merged_macroalg_dna

#Seagrass

prok_ra_seagr_dna_noCtrl<-subset_samples(subset_samples(subset_samples(subset_samples(prok_ra_all2, type %in% 'dna'), 
                                             group %in% 'sample'),
                                             rep %in% 'R1'),
                                             substrate2 %in% 'Seagrass')

# remove controls, divide dna and rna samples and get the list for the venn diagrams
prok_ra_seagr_dna_list<-get_vennlist(prok_ra_seagr_dna_noCtrl,factorNames='time')

ggvenn(
  prok_ra_seagr_dna_list, 
  fill_color = c("#337538ff","#337538ff","#337538ff"),
  stroke_size = 0.5, set_name_size = 10, text_size=8
  ) + theme_glab()

prok_ra_merged_seagr_dna <- merge_samples(prok_ra_seagr_dna_noCtrl,"time")
prok_ra_merged_seagr_dna

#RNA
#Bare

prok_ra_bare_rna_noCtrl<-subset_samples(subset_samples(subset_samples(subset_samples(prok_ra_all2, type %in% 'rna'), 
                                             group %in% 'sample'),
                                             substrate2 %in% 'Bare')

# remove controls, divide dna and rna samples and get the list for the venn diagrams
prok_ra_bare_rna_list<-get_vennlist(prok_ra_bare_rna_noCtrl,factorNames='time')

ggvenn(
  prok_ra_bare_rna_list, 
  fill_color = c("#dccd7dc2","#dccd7dc2","#dccd7dc2"),
  stroke_size = 0.5, set_name_size = 10, text_size=8
  ) + theme_glab()

prok_ra_merged_bare_rna <- merge_samples(prok_ra_bare_rna_noCtrl,"time")
prok_ra_merged_bare_rna

#Macroalgae

prok_ra_macroalg_rna_noCtrl<-subset_samples(subset_samples(subset_samples(prok_ra_all2, type %in% 'rna'), 
                                             group %in% 'sample'),
                                             substrate2 %in% 'Macroalgae')

# remove controls, divide dna and rna samples and get the list for the venn diagrams
prok_ra_macroalg_rna_list<-get_vennlist(prok_ra_macroalg_rna_noCtrl,factorNames='time')

ggvenn(
  prok_ra_macroalg_rna_list, 
  fill_color = c("#7e2954c0","#7e2954c0","#7e2954c0"),
  stroke_size = 0.5, set_name_size = 10, text_size=8
  ) + theme_glab()

prok_ra_merged_macroalg_rna <- merge_samples(prok_ra_macroalg_rna_noCtrl,"time")
prok_ra_merged_macroalg_rna

#Seagrass

prok_ra_seagr_rna_noCtrl<-subset_samples(subset_samples(subset_samples(prok_ra_all2, type %in% 'rna'), 
                                             group %in% 'sample'),
                                             substrate2 %in% 'Seagrass')

# remove controls, divide dna and rna samples and get the list for the venn diagrams
prok_ra_seagr_rna_list<-get_vennlist(prok_ra_seagr_rna_noCtrl,factorNames='time')

ggvenn(
  prok_ra_seagr_rna_list, 
  fill_color = c("#337538ff","#337538ff","#337538ff"),
  stroke_size = 0.5, set_name_size = 10, text_size=8
  ) + theme_glab()

prok_ra_merged_seagr_rna <- merge_samples(prok_ra_seagr_rna_noCtrl,"time")
prok_ra_merged_seagr_rna

save.image()

#Preparing figure 5

#DNA Core
pseq.core_bare_dna <- core(prok_ra_merged_bare_dna, detection = 0, prevalence = .75)
pseq.core_macroalg_dna <- core(prok_ra_merged_macroalg_dna, detection = 0, prevalence = .75)
pseq.core_seagr_dna <- core(prok_ra_merged_seagr_dna, detection = 0, prevalence = .75)

#RNA Core
pseq.core_bare_rna <- core(prok_ra_merged_bare_rna, detection = 0, prevalence = .75)
pseq.core_macroalg_rna <- core(prok_ra_merged_macroalg_rna, detection = 0, prevalence = .75)
pseq.core_seagr_rna <- core(prok_ra_merged_seagr_rna, detection = 0, prevalence = .75)

library(MicEco)

sample_data(pseq.core_bare_dna)$substrate2<-c("Bare_DNA","Bare_DNA","Bare_DNA")
sample_data(pseq.core_macroalg_dna)$substrate2<-c("Macroalgae_DNA","Macroalgae_DNA","Macroalgae_DNA")
sample_data(pseq.core_seagr_dna)$substrate2<-c("Seagrass_DNA","Seagrass_DNA","Seagrass_DNA")

pseq.core_bare_dna_ra = transform_sample_counts(merge_samples(pseq.core_bare_dna,"substrate2"), function(x){x / sum(x)})                                                       
pseq.core_macroalg_dna_ra = transform_sample_counts(merge_samples(pseq.core_macroalg_dna,"substrate2"), function(x){x / sum(x)})                                                       
pseq.core_seagr_dna_ra = transform_sample_counts(merge_samples(pseq.core_seagr_dna,"substrate2"), function(x){x / sum(x)})    

pseq.core_bare_dna_class = tax_glom(pseq.core_bare_dna_ra, "Class", NArm = FALSE)
pseq.core_macroalg_dna_class = tax_glom(pseq.core_macroalg_dna_ra, "Class", NArm = FALSE)
pseq.core_seagr_dna_class = tax_glom(pseq.core_seagr_dna_ra, "Class", NArm = FALSE)

#DNA: Subset taxa with relative abundance > 1%
pseq.core_bare_dna_class_perc1 <- ps_prune(pseq.core_bare_dna_class, min.abundance = 0.01)
pseq.core_macroalg_dna_class_perc1 <- ps_prune(pseq.core_macroalg_dna_class, min.abundance = 0.01)
pseq.core_seagr_dna_class_perc1 <- ps_prune(pseq.core_seagr_dna_class, min.abundance = 0.01)

#RNA

sample_data(pseq.core_bare_rna)$substrate2<-c("Bare_RNA","Bare_RNA","Bare_RNA")
sample_data(pseq.core_macroalg_rna)$substrate2<-c("Macroalgae_RNA","Macroalgae_RNA","Macroalgae_RNA")
sample_data(pseq.core_seagr_rna)$substrate2<-c("Seagrass_RNA","Seagrass_RNA")

pseq.core_bare_rna_ra = transform_sample_counts(merge_samples(pseq.core_bare_rna,"substrate2"), function(x){x / sum(x)})                                                       
pseq.core_macroalg_rna_ra = transform_sample_counts(merge_samples(pseq.core_macroalg_rna,"substrate2"), function(x){x / sum(x)})                                                       
pseq.core_seagr_rna_ra = transform_sample_counts(merge_samples(pseq.core_seagr_rna,"substrate2"), function(x){x / sum(x)})                                                       

pseq.core_bare_rna_class = tax_glom(pseq.core_bare_rna_ra, "Class", NArm = FALSE)
pseq.core_macroalg_rna_class = tax_glom(pseq.core_macroalg_rna_ra, "Class", NArm = FALSE)
pseq.core_seagr_rna_class = tax_glom(pseq.core_seagr_rna_ra, "Class", NArm = FALSE)

#RNA: Subset taxa with relative abundance > 1%
pseq.core_bare_rna_class_perc1 <- ps_prune(pseq.core_bare_rna_class, min.abundance = 0.01)
pseq.core_macroalg_rna_class_perc1 <- ps_prune(pseq.core_macroalg_rna_class, min.abundance = 0.01)
pseq.core_seagr_rna_class_perc1 <- ps_prune(pseq.core_seagr_rna_class, min.abundance = 0.01)

#Clean taxa names
pseq.core_bare_rna_class_perc1 <-ps_tax_clean(pseq.core_bare_rna_class_perc1)
pseq.core_macroalg_rna_class_perc1 <-ps_tax_clean(pseq.core_macroalg_rna_class_perc1)
pseq.core_seagr_rna_class_perc1 <-ps_tax_clean(pseq.core_seagr_rna_class_perc1)

#Donuts plot

# Saving a .csv with combined taxonomy and ASVs abundances for the selected taxonomic level
#DNA BARE
summary_object <- pseq.core_bare_dna_class_perc1
summary_otu <- as.matrix(t(otu_table(summary_object)))
summary_tax <- as.matrix(tax_table(summary_object))
bare_core_dna_1perc_summary <- cbind(summary_otu, summary_tax)
bare_core_dna_1perc_summary <- data.frame(bare_core_dna_1perc_summary[,c(1,4)])
bare_core_dna_1perc_summary<-bare_core_dna_1perc_summary %>% arrange(is.na(Class), Class)
bare_core_dna_1perc_summary$ymax = cumsum(bare_core_dna_1perc_summary$Bare_DNA)
bare_core_dna_1perc_summary$ymin = c(0, head(bare_core_dna_1perc_summary$ymax, n=-1))
bare_core_dna_1perc_summary[,2] <- na_if(bare_core_dna_1perc_summary[,2], "NA_NA")
bare_core_dna_1perc_summary <- bare_core_dna_1perc_summary %>% arrange(is.na(.[[2]]), .[[2]])                                                                      
bare_core_dna_1perc_summary

#DNA MACROALGAE
summary_object <- pseq.core_macroalg_dna_class_perc1
summary_otu <- as.matrix(t(otu_table(summary_object)))
summary_tax <- as.matrix(tax_table(summary_object))
macroalg_core_dna_1perc_summary <- cbind(summary_otu, summary_tax)
macroalg_core_dna_1perc_summary <- data.frame(macroalg_core_dna_1perc_summary[,c(1,4)])
macroalg_core_dna_1perc_summary<-macroalg_core_dna_1perc_summary %>% arrange(is.na(Class), Class)
macroalg_core_dna_1perc_summary$ymax = cumsum(macroalg_core_dna_1perc_summary$Macroalgae_DNA)
macroalg_core_dna_1perc_summary$ymin = c(0, head(macroalg_core_dna_1perc_summary$ymax, n=-1))
macroalg_core_dna_1perc_summary[,2] <- na_if(macroalg_core_dna_1perc_summary[,2], "NA_NA")
macroalg_core_dna_1perc_summary <- macroalg_core_dna_1perc_summary %>% arrange(is.na(.[[2]]), .[[2]])                                                                      
macroalg_core_dna_1perc_summary

#DNA SEAGRASS
summary_object <- pseq.core_seagr_dna_class_perc1
summary_otu <- as.matrix(t(otu_table(summary_object)))
summary_tax <- as.matrix(tax_table(summary_object))
seagr_core_dna_1perc_summary <- cbind(summary_otu, summary_tax)
seagr_core_dna_1perc_summary <- data.frame(seagr_core_dna_1perc_summary[,c(1,4)])
seagr_core_dna_1perc_summary<-seagr_core_dna_1perc_summary %>% arrange(is.na(Class), Class)
seagr_core_dna_1perc_summary$ymax = cumsum(seagr_core_dna_1perc_summary$Seagrass_DNA)
seagr_core_dna_1perc_summary$ymin = c(0, head(seagr_core_dna_1perc_summary$ymax, n=-1))
seagr_core_dna_1perc_summary[,2] <- na_if(seagr_core_dna_1perc_summary[,2], "NA_NA")
seagr_core_dna_1perc_summary <- seagr_core_dna_1perc_summary %>% arrange(is.na(.[[2]]), .[[2]])                                                                      
seagr_core_dna_1perc_summary

#RNA BARE
summary_object <- pseq.core_bare_rna_class_perc1
summary_otu <- as.matrix(t(otu_table(summary_object)))
summary_tax <- as.matrix(tax_table(summary_object))
bare_core_rna_1perc_summary <- cbind(summary_otu, summary_tax)
bare_core_rna_1perc_summary <- data.frame(bare_core_rna_1perc_summary[,c(1,4)])
bare_core_rna_1perc_summary[,2] <- na_if(bare_core_rna_1perc_summary[,2], "NA_NA")
bare_core_rna_1perc_summary <- bare_core_rna_1perc_summary %>% arrange(is.na(.[[2]]), .[[2]])                                                                      
bare_core_rna_1perc_summary <- bare_core_rna_1perc_summary[c(1, 3, 2, 4:nrow(bare_core_rna_1perc_summary)), ]
bare_core_rna_1perc_summary$ymax = cumsum(bare_core_rna_1perc_summary$Bare_RNA)
bare_core_rna_1perc_summary$ymin = c(0, head(bare_core_rna_1perc_summary$ymax, n=-1))
bare_core_rna_1perc_summary

#RNA MACROALGAE
summary_object <- pseq.core_macroalg_rna_class_perc1
summary_otu <- as.matrix(t(otu_table(summary_object)))
summary_tax <- as.matrix(tax_table(summary_object))
macroalg_core_rna_1perc_summary <- cbind(summary_otu, summary_tax)
macroalg_core_rna_1perc_summary <- data.frame(macroalg_core_rna_1perc_summary[,c(1,4)])
macroalg_core_rna_1perc_summary[,2] <- na_if(macroalg_core_rna_1perc_summary[,2], "NA_NA")
macroalg_core_rna_1perc_summary <- macroalg_core_rna_1perc_summary %>% arrange(is.na(.[[2]]), .[[2]])                                                                      
macroalg_core_rna_1perc_summary$ymax = cumsum(macroalg_core_rna_1perc_summary$Macroalgae_RNA)
macroalg_core_rna_1perc_summary$ymin = c(0, head(macroalg_core_rna_1perc_summary$ymax, n=-1))
macroalg_core_rna_1perc_summary

#RNA SEAGRASS
summary_object <- pseq.core_seagr_rna_class_perc1
summary_otu <- as.matrix(t(otu_table(summary_object)))
summary_tax <- as.matrix(tax_table(summary_object))
seagr_core_rna_1perc_summary <- cbind(summary_otu, summary_tax)
seagr_core_rna_1perc_summary <- data.frame(seagr_core_rna_1perc_summary[,c(1,4)])
seagr_core_rna_1perc_summary[,2] <- na_if(seagr_core_rna_1perc_summary[,2], "NA_NA")
seagr_core_rna_1perc_summary <- seagr_core_rna_1perc_summary %>% arrange(is.na(.[[2]]), .[[2]])                                                                      
seagr_core_rna_1perc_summary$ymax = cumsum(seagr_core_rna_1perc_summary$Seagrass_RNA)
seagr_core_rna_1perc_summary$ymin = c(0, head(seagr_core_rna_1perc_summary$ymax, n=-1))
seagr_core_rna_1perc_summary

unique_classes <-unique(c(
                          bare_core_dna_1perc_summary$Class,
                          bare_core_rna_1perc_summary$Class,
                          macroalg_core_dna_1perc_summary$Class,
                          macroalg_core_rna_1perc_summary$Class,
                          seagr_core_dna_1perc_summary$Class,
                          seagr_core_rna_1perc_summary$Class
                          )
                       )

sorted_classes <- sort(unique_classes[unique_classes != "NA"])
sorted_classes <- c(sorted_classes, "NA")
sorted_classes

# Generate Spectral palette for real classes
n_colors <- length(sorted_classes)
has_na <- any(is.na(sorted_classes))

# If more than 11 classes, interpolate (we have 20 classes and 1 NA that stay for Others)
if (n_colors > 11) {
  palette_colors <- colorRampPalette(brewer.pal(11, "Spectral"))(n_colors)
}

# Assign names
class_colors <- setNames(palette_colors, sorted_classes)

# Add NA as black if present
if (has_na) {
  class_colors <- c(class_colors, "black")
  names(class_colors)[length(class_colors)] <- "NA"
}

ggplot(bare_core_dna_1perc_summary, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Class)) +
     geom_rect() +
     scale_fill_manual(na.value="black",values = class_colors) +
     coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
     xlim(c(2, 4)) + # Try to remove that to see how to make a pie chart
     theme_glab()

ggplot(macroalg_core_dna_1perc_summary, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Class)) +
     geom_rect() +
     scale_fill_manual(na.value="black",values = class_colors) +
     coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
     xlim(c(2, 4)) + # Try to remove that to see how to make a pie chart
     theme_glab()

ggplot(seagr_core_dna_1perc_summary, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Class)) +
     geom_rect() +
     scale_fill_manual(na.value="black",values = class_colors) +
     coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
     xlim(c(2, 4)) + # Try to remove that to see how to make a pie chart
     theme_glab()

ggplot(bare_core_rna_1perc_summary, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Class)) +
     geom_rect() +
     scale_fill_manual(na.value="black",values = class_colors) +
     coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
     xlim(c(2, 4)) + # Try to remove that to see how to make a pie chart
     theme_glab()

ggplot(macroalg_core_rna_1perc_summary, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Class)) +
     geom_rect() +
     scale_fill_manual(na.value="black",values = class_colors) +
     coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
     xlim(c(2, 4)) + # Try to remove that to see how to make a pie chart
     theme_glab()

ggplot(seagr_core_rna_1perc_summary, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Class)) +
     geom_rect() +
     scale_fill_manual(na.value="black",values = class_colors) +
     coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
     xlim(c(2, 4)) + # Try to remove that to see how to make a pie chart
     theme_glab()

save.image()



#Preparing dataframe for RNA:DNA plot

bare_rna_dna_df<-na.omit(left_join(bare_core_rna_1perc_summary[,1:2],bare_core_dna_1perc_summary[,1:2], by='Class'))  %>%
    mutate(rna_dna = as.numeric(Bare_RNA)/as.numeric(Bare_DNA))

macroalg_rna_dna_df<-na.omit(left_join(macroalg_core_rna_1perc_summary[,1:2],macroalg_core_dna_1perc_summary[,1:2], by='Class'))  %>%
    mutate(rna_dna = as.numeric(Macroalgae_RNA)/as.numeric(Macroalgae_DNA))

seagr_rna_dna_df<-na.omit(left_join(seagr_core_rna_1perc_summary[,1:2],seagr_core_dna_1perc_summary[,1:2], by='Class'))  %>%
    mutate(rna_dna = as.numeric(Seagrass_RNA)/as.numeric(Seagrass_DNA))


bare_rna_dna_df

macroalg_rna_dna_df

seagr_rna_dna_df

ggplot(bare_rna_dna_df, aes(x=Class,y=sort(log10(rna_dna)))) + 
geom_bar(aes(fill=Class),stat = "identity",linewidth=0.446,color="black",) +
ylab(expression(log[10](RNA:DNA))) + xlab("") +
scale_fill_manual(na.value="black",values = class_colors) +
theme_glab() + theme(legend.position = "none", aspect.ratio=0.7, axis.text.x = element_text(angle = 45, hjust=0.95))

ggplot(macroalg_rna_dna_df, aes(x=Class,y=sort(log10(rna_dna)))) + 
geom_bar(aes(fill=Class),stat = "identity",linewidth=0.446,color="black",) +
ylab(expression(log[10](RNA:DNA))) + xlab("") +
scale_fill_manual(na.value="black",values = class_colors) +
theme_glab() + theme(legend.position = "none", aspect.ratio=0.7, axis.text.x = element_text(angle = 45, hjust=0.95))

ggplot(seagr_rna_dna_df, aes(x=Class,y=sort(log10(rna_dna)))) + 
geom_bar(aes(fill=Class),stat = "identity",linewidth=0.446,color="black",) +
ylab(expression(log[10](RNA:DNA))) + xlab("") +
scale_fill_manual(na.value="black",values = class_colors) +
theme_glab() + theme(legend.position = "none", aspect.ratio=0.7, axis.text.x = element_text(angle = 45, hjust=0.95))

save.image()





#supplementary figure 4

#RNA BETA DIVERSITY

prok_ndata_rna <- transform_sample_counts(subset_samples(prok_ndata_all,type=='rna'),
                                          function(x) ((x / sum(x))*median(readcount(subset_samples(prok_ndata_all,type=='rna')))))

#calc distance matrix: Weighted Jaccard index
prok_dist_jac_rna <- distance(prok_ndata_rna, method = "jaccard")
prok_dist_uwjac_rna <- distance(prok_ndata_rna, method = "jaccard",binary=TRUE)

# NMDS Jaccard weighted
prok_nmds_jac_rna <- ordinate(prok_ndata_rna,prok_dist_jac_rna, method = "NMDS", trymax=9999)

# NMDS Jaccard Unweighted
prok_nmds_juw_rna <- ordinate(prok_ndata_rna,prok_dist_uwjac_rna, method = "NMDS", trymax=9999)

p_nmds_wj_rna<-(plot_ordination(prok_ndata_rna,prok_nmds_jac_rna, type="samples"))$data
p_nmds_uwj_rna<-(plot_ordination(prok_ndata_rna,prok_nmds_juw_rna, type="samples"))$data

p_nmds_wj_rna[c(2,3,14),7] <- "T1C"
p_nmds_wj_rna[c(4,5,6),7] <- "T2C"
#p_nmds_wj_rna


p_nmds_uwj_rna[c(2,3,14),7] <- "T1C"
p_nmds_uwj_rna[c(4,5,6),7] <- "T2C"
#p_nmds_uwj_rna

#figure s4_a_b

ggarrange(

####weighted

ggplot(p_nmds_wj_rna, aes(x=NMDS1, y=NMDS2)) +
geom_point(aes(fill=substrate2,color=substrate2,shape=time),size=8,stroke=0.3) +
#geom_text(aes(label=name),color="black", size=4,hjust=0.5, vjust=-0.6) + 
scale_fill_manual("Substrate",values = c("Bare" = "#dccd7dc2", "Macroalgae" = "#7e2954c0",
                               "Seagrass" = "#337538ff")) +
scale_color_manual(values = c("Bare" = "#dccd7dc2", "Macroalgae" = "#7e2954c0",
                               "Seagrass" = "#337538ff")) +
scale_shape_manual("Time",values = c("T0" = 21, "T1" = 22, "T2" = 23, "T1C"=24 , "T2C" =25))  +
    
theme_glab() + theme(legend.position = "bottom") #+ guides(fill = guide_legend(override.aes = list(shape = 21)))

,
    
####unweighted

ggplot(p_nmds_uwj_rna, aes(x=NMDS1, y=NMDS2)) +
geom_point(aes(fill=substrate2,color=substrate2,shape=time),size=8,stroke=0.3) +
#geom_text(aes(label=name,),color="black", size=4,hjust=0.5, vjust=-0.6) + 
scale_fill_manual("Substrate",values = c("Bare" = "#dccd7dc2", "Macroalgae" = "#7e2954c0",
                               "Seagrass" = "#337538ff")) +
scale_color_manual(values = c("Bare" = "#dccd7dc2", "Macroalgae" = "#7e2954c0",
                               "Seagrass" = "#337538ff")) +
scale_shape_manual("Time",values = c("T0" = 21, "T1" = 22, "T2" = 23, "T1C"=24 , "T2C" =25))  +
    
theme_glab() + theme(legend.position = "bottom") #+ guides(fill = guide_legend(override.aes = list(shape = 21)))

,
ncol=2, nrow=2, legend="right", 
    labels= c('(A)', '(B)'), font.label = list(size = 18,face="plain",color = "black"),

    align = c("hv"), common.legend = TRUE)
    
#ggsave("../../Venice_Lagoon_figures_tables/Venice_Lagoon_figures_tables//figures/nmds_w_uw_cdna.svg",width=12,height=8)

adonis2(prok_dist_jac_rna~time,permutations=999,
    data = p_nmds_wj_rna)

adonis2(prok_dist_uwjac_rna~time,permutations=999,
    data = p_nmds_uwj_rna)

save.image()

#RNA NO CTRL

prok_ndata_rna_noCTRL <- transform_sample_counts(subset_samples(subset_samples(prok_ndata_all,type=='rna'), group=='sample'),
                                          function(x) ((x / sum(x))*median(readcount(subset_samples(subset_samples(prok_ndata_all,type=='rna'), 
                                                                                                    group=='sample')))))

#calc distance matrix: Weighted Jaccard index
prok_dist_jac_rna_noCTRL <- distance(prok_ndata_rna_noCTRL, method = "jaccard")
prok_dist_uwjac_rna_noCTRL <- distance(prok_ndata_rna_noCTRL, method = "jaccard",binary=TRUE)

# NMDS Jaccard weighted
prok_nmds_jac_rna_noCTRL <- ordinate(prok_ndata_rna,prok_dist_jac_rna_noCTRL, method = "NMDS", trymax=9999)

# NMDS Jaccard Unweighted
prok_nmds_juw_rna_noCTRL <- ordinate(prok_ndata_rna,prok_dist_uwjac_rna_noCTRL, method = "NMDS", trymax=9999)

p_nmds_wj_rna_noCTRL<-(plot_ordination(prok_ndata_rna,prok_nmds_jac_rna_noCTRL, type="samples"))$data
p_nmds_uwj_rna_noCTRL<-(plot_ordination(prok_ndata_rna,prok_nmds_juw_rna_noCTRL, type="samples"))$data

####weighted

ggplot(p_nmds_wj_rna_noCTRL, aes(x=NMDS1, y=NMDS2)) +
geom_point(aes(fill=substrate2,color=substrate2,shape=time),size=8,stroke=0.3) +
geom_text(aes(label=name,),color="black", size=4,hjust=0.5, vjust=-0.6) + 
scale_fill_manual("Substrate",values = c("Bare" = "#dccd7dc2", "Macroalgae" = "#7e2954c0",
                               "Seagrass" = "#337538ff")) +
scale_color_manual(values = c("Bare" = "#dccd7dc2", "Macroalgae" = "#7e2954c0",
                               "Seagrass" = "#337538ff")) +
scale_shape_manual("Time",values = c("T0" = 21, "T1" = 22, "T2" = 23))  +

theme_glab() + theme(legend.position = "right") + 
guides(fill = guide_legend(override.aes = list(shape = 21)))


####unweighted

ggplot(p_nmds_uwj_rna_noCTRL, aes(x=NMDS1, y=NMDS2)) +
geom_point(aes(fill=substrate2,color=substrate2,shape=time),size=8,stroke=0.3) +
geom_text(aes(label=name,),color="black", size=4,hjust=0.5, vjust=-0.6) + 
scale_fill_manual("Substrate",values = c("Bare" = "#dccd7dc2", "Macroalgae" = "#7e2954c0",
                               "Seagrass" = "#337538ff")) +
scale_color_manual(values = c("Bare" = "#dccd7dc2", "Macroalgae" = "#7e2954c0",
                               "Seagrass" = "#337538ff")) +
scale_shape_manual("Time",values = c("T0" = 21, "T1" = 22, "T2" = 23))  +

theme_glab() + theme(legend.position = "right") + 
guides(fill = guide_legend(override.aes = list(shape = 21)))

#ggsave("plot/nmds_jw_polygon.svg",width=12,height=8)

save.image()

#supplementary figure 1 

library(polygonPlot)

head(vl_prj_env_data_full)

bare_polygon_df<-prepare_dataframe(vl_prj_env_data_full %>%
  filter(time == "T0",
         type == "dna",
         substrate2 == "Bare") %>%
  select(9:13))

macroalg_polygon_df<-prepare_dataframe(vl_prj_env_data_full %>%
  filter(time == "T0",
         type == "dna",
         substrate2 == "Macroalgae") %>%
  select(9:13))

seagr_polygon_df<-prepare_dataframe(vl_prj_env_data_full %>%
  filter(time == "T0",
         type == "dna",
         substrate2 == "Seagrass") %>%
  select(9:13))

df_list <- list(bare_polygon_df,macroalg_polygon_df,seagr_polygon_df)

all_data <- do.call(
  rbind,
  lapply(df_list, function(df) {
    df[!df$info %in% c("axis_min", "axis_max"), ]
  })
)

limits <- sapply(
  all_data[ , -1],   # remove "info"
  function(x) {
    r <- range(x, na.rm = TRUE)
    span <- diff(r)
    c(
      axis_min = floor(r[1] - 0.05 * span),
      axis_max = ceiling(r[2] + 0.05 * span)
    )
  }
)


colnames(limits)[colnames(limits) == "phaeopigments"] <- "Phaeo"
colnames(limits)[colnames(limits) == "chlorophyll"]   <- "Chl-a"
colnames(limits)[colnames(limits) == "prt"] <- "Prt"
colnames(limits)[colnames(limits) == "cho"]   <- "Cho"
colnames(limits)[colnames(limits) == "lip"]   <- "Lip"

limits

update_axes <- function(df, limits) {

  min_row <- which(!is.na(df$info) & df$info == "axis_min")
  max_row <- which(!is.na(df$info) & df$info == "axis_max")

  df[min_row, -1] <- limits["axis_min", ]
  df[max_row, -1] <- limits["axis_max", ]

  df
}

df_list <- lapply(df_list, update_axes, limits = limits)

df_list <- lapply(df_list, function(df) {
  num_cols <- sapply(df, is.numeric)
  df[num_cols] <- lapply(df[num_cols], round, 2)
  df
})

df_list <- lapply(df_list, function(df) {
  names(df)[names(df) == "prt"] <- "Prt"
  names(df)[names(df) == "cho"] <- "Cho"
  names(df)[names(df) == "lip"] <- "Lip"
  names(df)[names(df) == "phaeopigments"] <- "Pheo"
  names(df)[names(df) == "chlorophyll"]   <- "Chl-a"
  df
})

bare_polyg_df <- df_list[[1]]
head(bare_polyg_df)
macroalg_polyg_df <- df_list[[2]]
head(macroalg_polyg_df)
seagr_polyg_df <- df_list[[3]]
head(seagr_polyg_df)

#figure s1_a

ggarrange(

polygonPlot(bare_polyg_df, shape = 5, fillcolor = "#dccd7dc2",
            labels_axis = c(expression(Prt~(mg~g^{-1})), expression(Cho~(mg~g^{-1})), expression(Lip~(mg~g^{-1})), 
                            expression(Pheo~(mu*g~g^{-1})), expression(Chl-a~(mu*g~g^{-1}))), 
            title="Bare") ,
    
polygonPlot(macroalg_polyg_df, shape = 5, fillcolor = "#7e2954c0",
            labels_axis = c(expression(Prt~(mg~g^{-1})), expression(Cho~(mg~g^{-1})), expression(Lip~(mg~g^{-1})), 
                            expression(Pheo~(mu*g~g^{-1})), expression(Chl-a~(mu*g~g^{-1}))), 
            title="Macroalgae") ,

polygonPlot(seagr_polyg_df, shape = 5, fillcolor = "#337538ff",
            labels_axis = c(expression(Prt~(mg~g^{-1})), expression(Cho~(mg~g^{-1})), expression(Lip~(mg~g^{-1})), 
                            expression(Pheo~(mu*g~g^{-1})), expression(Chl-a~(mu*g~g^{-1}))), 
            title="Seagrass") ,

ncol=3, nrow=1, align="hv", common.legend = TRUE, legend="bottom")

#ggsave("../plot/polygonplot.svg", width=12, height=12)

save.image()

#figure s1_b

ggarrange(
ggplot(vl_prj_env_data_full, aes(x=substrate2,y=bpc)) +
geom_boxplot(aes(fill=substrate2)) +
scale_fill_manual("Substrate",values = c("Bare" = "#dccd7dc2", "Macroalgae" = "#7e2954c0",
                               "Seagrass" = "#337538ff")) + xlab("") + ylab("BPC") + scale_x_discrete(limits=rev) +
theme_glab() + theme(aspect.ratio=0.15) + coord_flip() ,

ggplot(vl_prj_env_data_full, aes(x=substrate2,y=prt/cho)) +
geom_boxplot(aes(fill=substrate2)) +
scale_fill_manual("Substrate",values = c("Bare" = "#dccd7dc2", "Macroalgae" = "#7e2954c0",
                               "Seagrass" = "#337538ff")) + xlab("") + ylab("Prt:Cho") + scale_x_discrete(limits=rev) +
theme_glab() + theme(aspect.ratio=0.15) + coord_flip() ,

ggplot(vl_prj_env_data_full, aes(x=substrate2,y=chlorophyll/phaeopigments)) +
geom_boxplot(aes(fill=substrate2)) +
scale_fill_manual("Substrate",values = c("Bare" = "#dccd7dc2", "Macroalgae" = "#7e2954c0",
                               "Seagrass" = "#337538ff")) + xlab("") + ylab("Chl-a:Pheo") + scale_x_discrete(limits=rev) +
theme_glab() + theme(aspect.ratio=0.15) + coord_flip() ,
ncol=1, nrow=3, align="hv", common.legend = TRUE, legend="bottom")
#ggsave("../plot/om_boxplots.svg", width=16, height=16)

save.image()
