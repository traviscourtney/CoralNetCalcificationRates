######################################################################
##’ @title Area-normalized scaling of ReefBudget calcification and bioerosion rates for use with CoralNet
##’
##’ @author Travis A Courtney
##’ @contact traviscourtney@gmail.com
##’ @date 2021-06–28
##’ @reference Courtney TA, Chan S, Lange ID, Perry CT, Andersson AJ (2021) Area-normalized ReefBudget calcification and bioerosion rates for use with CoralNet
##’ @log Version 1

###########include or exclude bioerosion##########
#bioerosion rates are included by default, set include_bioerosion to FALSE to calculate gross carbonate production only
include_bioerosion=TRUE

###########require packages and scripts##########
library(fs)
library(tidyverse)
library(readxl)
library(sn)
source("ReefBudgetRandCalc.R")

#turn off scientific notation
options(scipen=999)

#########CoralNet Indo-Pacific Calcification Rates#########
#Import Required Data
setwd("Data")
setwd("NOAA_Pacific_Coral_Demography") #Import NOAA Pacific Coral Demograhy Data from multiple files
NOAA_demography_raw=
  dir_ls(regexp = "\\.csv$") %>% 
  map_dfr(read_csv, .id = "source")
setwd("..")
NOAA_CRED_CoralNet_labelset <- read_excel("NOAA_CRED_CoralNet_labelset.xlsx") #Import NOAA CRED CoralNet labelset
Gonzalez_Barrios <- read_csv("Gonzalez-Barios_Alvarez-Filip_2019_Rugosity.csv") #Import Gonzalez-Barrios taxa-level estimates of colony structural complexity
ReefBudget_Pacific <- read_excel("Indo-Pacific_Carbonate_Production_v1.2.xlsx", sheet = "Calcification Rates", skip = 9) #import Indo-Pacific ReefBudget rates
macromicro_bioerosion=read_excel("Indo-Pacific_Carbonate_Production_v1.2.xlsx", sheet = "Macro- & Microbioerosion") #import Indo_Pacific ReefBudget bioerosion rates
duplicates=read_csv("CoralNetDuplicates.csv") #import coralnet duplicate files
setwd("..")

#remove missing data points
NOAA_demography = subset(NOAA_demography_raw,COLONYLENGTH!=-9)

#re-assign labels to match CoralNet ID's
NOAA_demography$MORPHOLOGY=str_replace(NOAA_demography$MORPHOLOGY,"Mounding","Massive")
NOAA_demography$MORPHOLOGY=str_replace(NOAA_demography$MORPHOLOGY,"Encrusting \\(columnar\\)","Encrusting")
NOAA_demography$MORPHOLOGY=str_replace(NOAA_demography$MORPHOLOGY,"Encrusting \\(mounding\\)","Encrusting")
NOAA_demography$MORPHOLOGY=str_replace(NOAA_demography$MORPHOLOGY,"Encrusting \\(flat\\)","Encrusting")
NOAA_demography$MORPHOLOGY=str_replace(NOAA_demography$MORPHOLOGY,"Disc - free living","Free-living")
NOAA_demography$MORPHOLOGY=tolower(NOAA_demography$MORPHOLOGY)
unique(NOAA_demography$MORPHOLOGY)

#find mean colony lengths for each genus
NOAA_genus_lengths=
  NOAA_demography %>%
  group_by(GENUS) %>%
  summarise(length_mean=quantile(COLONYLENGTH, probs = 0.5, na.rm=TRUE),length_lwr=quantile(COLONYLENGTH, probs = 0.025, na.rm=TRUE),length_upr=quantile(COLONYLENGTH, probs = 0.975, na.rm=TRUE))
NOAA_genus_lengths=NOAA_genus_lengths[!is.na(NOAA_genus_lengths$GENUS),]

#find mean colony lengths for each genus, morphology
NOAA_genus_morphology_lengths=
  NOAA_demography %>%
  group_by(GENUS,MORPHOLOGY) %>%
  summarise(length_mean=quantile(COLONYLENGTH, probs = 0.5, na.rm=TRUE),length_lwr=quantile(COLONYLENGTH, probs = 0.025, na.rm=TRUE),length_upr=quantile(COLONYLENGTH, probs = 0.975, na.rm=TRUE))
NOAA_genus_morphology_lengths=NOAA_genus_morphology_lengths[!is.na(NOAA_genus_morphology_lengths$GENUS),]

#find mean colony lengths for each morphology
NOAA_morphology_lengths=
  NOAA_demography %>%
  group_by(MORPHOLOGY) %>%
  summarise(length_mean=quantile(COLONYLENGTH, probs = 0.5, na.rm=TRUE),length_lwr=quantile(COLONYLENGTH, probs = 0.025, na.rm=TRUE),length_upr=quantile(COLONYLENGTH, probs = 0.975, na.rm=TRUE))
NOAA_morphology_lengths=NOAA_morphology_lengths[!is.na(NOAA_morphology_lengths$MORPHOLOGY),]

#Import NOAA CRED labels
NOAA_CRED_CoralNet_labelset$TAXA = substr(NOAA_CRED_CoralNet_labelset$Name, start = 6, stop=1000)
NOAA_CRED_CoralNet_labelset$GENUS = word(NOAA_CRED_CoralNet_labelset$TAXA,1)
NOAA_CRED_CoralNet_labelset$MORPHOLOGY = substr(word(NOAA_CRED_CoralNet_labelset$TAXA,2),start=5,stop=1000)
colnames(NOAA_CRED_CoralNet_labelset)=c("NAME","CODE","GROUP","TAXA","GENUS","MORPHOLOGY")
NOAA_CRED_morph_lengths=merge(subset(NOAA_CRED_CoralNet_labelset,GROUP=="Hard coral"), NOAA_genus_morphology_lengths, by=c("GENUS","MORPHOLOGY"),all.x=TRUE)

#fill in missing species data with genera means
NOAA_CRED_morph_lengths_all=NOAA_CRED_morph_lengths[!is.na(NOAA_CRED_morph_lengths$length_mean),]
NOAA_CRED_morph_lengths_na=NOAA_CRED_morph_lengths[is.na(NOAA_CRED_morph_lengths$length_mean),]
NOAA_CRED_morph_genus_lengths=merge(NOAA_CRED_morph_lengths_na[,c(1:6)],NOAA_genus_lengths,by="GENUS",all.x=TRUE)

#use consistent labels throughout
NOAA_CRED_morph_genus_lengths[NOAA_CRED_morph_genus_lengths$GENUS=="Foliose",c(2,7:9)]=NOAA_morphology_lengths[NOAA_morphology_lengths$MORPHOLOGY=="foliose",]
NOAA_CRED_morph_genus_lengths[NOAA_CRED_morph_genus_lengths$GENUS=="Branching",c(2,7:9)]=NOAA_morphology_lengths[NOAA_morphology_lengths$MORPHOLOGY=="branching",]
NOAA_CRED_morph_genus_lengths[NOAA_CRED_morph_genus_lengths$GENUS=="Columnar",c(2,7:9)]=NOAA_morphology_lengths[NOAA_morphology_lengths$MORPHOLOGY=="columnar",]
NOAA_CRED_morph_genus_lengths[NOAA_CRED_morph_genus_lengths$GENUS=="Encrusting",c(2,7:9)]=NOAA_morphology_lengths[NOAA_morphology_lengths$MORPHOLOGY=="encrusting",]
NOAA_CRED_morph_genus_lengths[NOAA_CRED_morph_genus_lengths$GENUS=="Free-living",c(2,7:9)]=NOAA_morphology_lengths[NOAA_morphology_lengths$MORPHOLOGY=="free-living",]
NOAA_CRED_morph_genus_lengths[NOAA_CRED_morph_genus_lengths$GENUS=="Massive",c(2,7:9)]=NOAA_morphology_lengths[NOAA_morphology_lengths$MORPHOLOGY=="massive",]
#count whether Goniopora or Alveopora has more records and use most frequent genus
nrow(subset(NOAA_demography,GENUS=="Goniopora")) #199 colonies
nrow(subset(NOAA_demography,GENUS=="Alveopora")) #1 colony
NOAA_CRED_morph_genus_lengths[NOAA_CRED_morph_genus_lengths$GENUS=="Goniopora/Alveopora",c(7:9)]=NOAA_genus_lengths[NOAA_genus_lengths$GENUS=="Goniopora",c(2:4)] #Goniopora is more common so use Goniopora
NOAA_CRED_morph_genus_lengths[NOAA_CRED_morph_genus_lengths$GENUS=="Diploastrea",c(2,7:9)]=NOAA_morphology_lengths[NOAA_morphology_lengths$MORPHOLOGY=="massive",] #assign massive length to Diploastrea

#find most common morphology for each genus
NOAA_typical_morphology=
  NOAA_demography %>%
  group_by(GENUS) %>%
  summarise(MORPHOLOGY=names(which.max(table(MORPHOLOGY))))

#pair most common morphology for each genera with CRED labels
NOAA_CRED_morph_genus_lengths_filled=NOAA_CRED_morph_genus_lengths[NOAA_CRED_morph_genus_lengths$MORPHOLOGY != "",]
NOAA_CRED_morph_genus_lengths_na=NOAA_CRED_morph_genus_lengths[NOAA_CRED_morph_genus_lengths$MORPHOLOGY == "",c(1,3:9)]
NOAA_CRED_morph_genus_lengths_typical=merge(NOAA_CRED_morph_genus_lengths_na,NOAA_typical_morphology,by="GENUS",all.x=TRUE)

#combine all morphology and length data
NOAA_CRED_morph_genus_lengths_all=rbind(NOAA_CRED_morph_lengths_all,NOAA_CRED_morph_genus_lengths_filled,NOAA_CRED_morph_genus_lengths_typical)

#substitute typical morphologies for Goniopora
NOAA_CRED_morph_genus_lengths_all[NOAA_CRED_morph_genus_lengths_all$GENUS=="Goniopora/Alveopora",2]="encrusting" #Goniopora is more common so use Goniopora

#calculate mean±SD rugosities for each morphology
GB_morphology_rugosites=
  Gonzalez_Barrios %>%
  group_by(morphology) %>%
  summarise(rugosity_m=mean(rugosity_mean),rugosity_sd=sd(rugosity_mean))
colnames(GB_morphology_rugosites)=c("MORPHOLOGY","rugosity_m","rugosity_sd")
GB_morphology_rugosites$MORPHOLOGY=tolower(GB_morphology_rugosites$MORPHOLOGY)

#add rugosities to dataframe
NOAA_mgr=merge(NOAA_CRED_morph_genus_lengths_all,GB_morphology_rugosites,all.x=TRUE)

#replace missing encrusting morphologies with mean±SD from all Gonzalez-Barrios rugosities
NOAA_mgr[NOAA_mgr$MORPHOLOGY=="encrusting",]$rugosity_m=mean(Gonzalez_Barrios$rugosity_mean)
NOAA_mgr[NOAA_mgr$MORPHOLOGY=="encrusting",]$rugosity_sd=sd(Gonzalez_Barrios$rugosity_mean)

#replace missing massive lobate morphologies with massive rugosities
NOAA_mgr[NOAA_mgr$MORPHOLOGY=="massive lobate",]$rugosity_m=subset(GB_morphology_rugosites,MORPHOLOGY=="massive")$rugosity_m
NOAA_mgr[NOAA_mgr$MORPHOLOGY=="massive lobate",]$rugosity_sd=subset(GB_morphology_rugosites,MORPHOLOGY=="massive")$rugosity_sd

#replace missing free-living morphologies with massive rugosities
NOAA_mgr[NOAA_mgr$MORPHOLOGY=="free-living",]$rugosity_m=subset(GB_morphology_rugosites,MORPHOLOGY=="massive")$rugosity_m
NOAA_mgr[NOAA_mgr$MORPHOLOGY=="free-living",]$rugosity_sd=subset(GB_morphology_rugosites,MORPHOLOGY=="massive")$rugosity_sd

#replace missing columnar morphologies with digitate rugosities
NOAA_mgr[NOAA_mgr$MORPHOLOGY=="columnar",]$rugosity_m=subset(GB_morphology_rugosites,MORPHOLOGY=="digitate")$rugosity_m
NOAA_mgr[NOAA_mgr$MORPHOLOGY=="columnar",]$rugosity_sd=subset(GB_morphology_rugosites,MORPHOLOGY=="digitate")$rugosity_sd

#replace missing plating morphologies with foliose rugosities
NOAA_mgr[NOAA_mgr$MORPHOLOGY=="plating",]$rugosity_m=subset(GB_morphology_rugosites,MORPHOLOGY=="foliose")$rugosity_m
NOAA_mgr[NOAA_mgr$MORPHOLOGY=="plating",]$rugosity_sd=subset(GB_morphology_rugosites,MORPHOLOGY=="foliose")$rugosity_sd

#replace missing tabulate rugosity from literature
NOAA_mgr[NOAA_mgr$MORPHOLOGY=="tabulate",]$rugosity_m=1.918 #Estimated rugosty @ 1cm chain length from WebplotDigitizer of Figure 3 from Knudby and LeDrew (2007) Measuring structural complexity on coral reefs 
NOAA_mgr[NOAA_mgr$MORPHOLOGY=="tabulate",]$rugosity_sd=0

#Import ReefBudget Indo-Pacific data
colnames(ReefBudget_Pacific)=c("code","genera","morphology","ext_mean","ext_sd","dens_mean","dens_sd","cf","coeff_mean","coeff_lwr","coeff_upr","int_mean","int_lwr","int_upr","notes")
ReefBudget_Pacific=ReefBudget_Pacific[!is.na(ReefBudget_Pacific$coeff_mean),]

NOAA_mgr$GENUS_MORPHOLOGY=paste(NOAA_mgr$GENUS,NOAA_mgr$MORPHOLOGY) #combine genus and morphology labels for consistent use
ReefBudget_Pacific$GENUS_MORPHOLOGY=paste(ReefBudget_Pacific$genera,ReefBudget_Pacific$morphology) #combine genus and morphology labels for consistent use
ReefBudget_Pacific_short = ReefBudget_Pacific[,c(16,9:14)] #select only columns needed

#Merge NOAA size and rugosity data with ReefBudget rates
NOAA_mgr_rates= merge(NOAA_mgr,ReefBudget_Pacific_short,by=c("GENUS_MORPHOLOGY"),all.x=TRUE)

#use consistent labels throughout
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Acropora tabulate",13:18]=ReefBudget_Pacific_short[ReefBudget_Pacific_short$GENUS_MORPHOLOGY=="Acropora table",2:7]
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Branching branching",13:18]=ReefBudget_Pacific_short[ReefBudget_Pacific_short$GENUS_MORPHOLOGY=="Hard coral branching",2:7]
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Columnar columnar",13:18]=ReefBudget_Pacific_short[ReefBudget_Pacific_short$GENUS_MORPHOLOGY=="Hard coral columnar",2:7]
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Foliose foliose",13:18]=ReefBudget_Pacific_short[ReefBudget_Pacific_short$GENUS_MORPHOLOGY=="Hard coral foliose",2:7]
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Encrusting encrusting",13:18]=ReefBudget_Pacific_short[ReefBudget_Pacific_short$GENUS_MORPHOLOGY=="Hard coral encrusting",2:7]
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Massive massive",13:18]=ReefBudget_Pacific_short[ReefBudget_Pacific_short$GENUS_MORPHOLOGY=="Hard coral massive",2:7]
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Fungia free-living",13:18]=ReefBudget_Pacific_short[ReefBudget_Pacific_short$GENUS_MORPHOLOGY=="Fungia freeliving",2:7]
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Goniopora/Alveopora encrusting",13:18]=ReefBudget_Pacific_short[ReefBudget_Pacific_short$GENUS_MORPHOLOGY=="Goniopora encrusting",2:7]
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Merulina foliose",13:18]=ReefBudget_Pacific_short[ReefBudget_Pacific_short$GENUS_MORPHOLOGY=="Merulina plating",2:7]
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Symphyllia massive lobate",13:18]=ReefBudget_Pacific_short[ReefBudget_Pacific_short$GENUS_MORPHOLOGY=="Symphyllia massive",2:7]
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Turbinaria foliose",13:18]=ReefBudget_Pacific_short[ReefBudget_Pacific_short$GENUS_MORPHOLOGY=="Turbinaria plating",2:7]

#create mean rates for all freeliving corals
ReefBudget_Pacific_freeliving_rates =
  subset(ReefBudget_Pacific,morphology=="freeliving") %>% #take all ReefBudget data
  summarize(coeff_mean=mean(coeff_mean),
            coeff_lwr=mean(coeff_lwr),
            coeff_upr=mean(coeff_upr),
            int_mean=mean(int_mean),
            int_lwr=mean(int_lwr), 
            int_upr=mean(int_upr)) 
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Free-living free-living",13:18]=ReefBudget_Pacific_freeliving_rates

#fill in missing ReefBudget rates with substitutions
Heliopora_branching=c(coeff_mean=0.6987,coeff_lwr=0.6384,coeff_upr=0.7614,int_mean=0,int_lwr=0,int_upr=0) #Fill in Heliopora growth info using Ext and Dens from Courtney et al. (in prep) filled into Pocillopora branching with submassive conversion factor in ReefBudget sheet
Millepora_columnar=c(coeff_mean=0.8174,coeff_lwr=0.4948,coeff_upr=1.1399,int_mean=0,int_lwr=0,int_upr=0) #substitute Millepora extension and density into "Hard Coral Columnar" cells in ReefBudget sheet
Porites_foliose=c(coeff_mean=0.1795,coeff_lwr=0.1112,coeff_upr=0.2554,int_mean=3.9329,int_lwr=2.3863,int_upr=5.7092) #substitute mean±SD of all Porites extension and density into "Hard Coral Foliose" cells in ReefBudget sheet
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Heliopora branching",13:18]=Heliopora_branching
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Millepora columnar",13:18]=Millepora_columnar
NOAA_mgr_rates[NOAA_mgr_rates$GENUS_MORPHOLOGY=="Porites foliose",13:18]=Porites_foliose

#apply conversion factor to all branching morphologies to account for open space
NOAA_mgr_rates$cf_lwr=ifelse( NOAA_mgr_rates$MORPHOLOGY=="branching" , 0.168 , 1 ) #Doszpot et al (2019) minimum proportion of 3D space occupied
NOAA_mgr_rates$cf_upr=ifelse( NOAA_mgr_rates$MORPHOLOGY=="branching" , 0.598 , 1 ) #Doszpot et al (2019) maximum proportion of 3D space occupied

#select only necessary columns
NOAA_coral_rates=NOAA_mgr_rates[,4:20]

#add in CCA, hard substrate, sponges, other growth forms besides CRED labels?
NOAA_CCA_rates=cbind(NOAA_CRED_CoralNet_labelset[str_detect(NOAA_CRED_CoralNet_labelset$NAME,"CCA"),1:4],length_mean=1,length_lwr=1,length_upr=1,NOAA_mgr[NOAA_mgr$TAXA=="Encrusting hard coral",10:11],ReefBudget_Pacific[str_detect(ReefBudget_Pacific$code,"CCA"),9:14],cf_lwr=1,cf_upr=1)
NOAA_all_rates=rbind(NOAA_coral_rates,NOAA_CCA_rates)

#combine macro and microbioerosion rates
macrobioerosion_pacific=ifelse(include_bioerosion==TRUE,-1,0)*as.numeric(macromicro_bioerosion[2,5])
microbioerosion_pacific=ifelse(include_bioerosion==TRUE,-1,0)*as.numeric(macromicro_bioerosion[20,5])
bioerosion_pacific=macrobioerosion_pacific+microbioerosion_pacific

#determines mean, lower, and upper rate in kg m^-2 yr^-1 for 100% cover of the respective coral using CARICOMP colony sizes and species-specific rugosity from Gonzalez-Barrios
NOAA_calc_rates=NULL
for (k in 1:nrow(NOAA_all_rates))
{
  nsim=10000
  NAME=NOAA_all_rates$NAME[k]
  CODE=NOAA_all_rates$CODE[k]
  GROUP=NOAA_all_rates$GROUP[k]
  TAXA=NOAA_all_rates$TAXA[k]
  s=ReefBudgetRandCalc(x=NOAA_all_rates$length_mean[k],lwr=NOAA_all_rates$length_lwr[k],upr=NOAA_all_rates$length_upr[k],n=nsim)
  c=ReefBudgetRandCalc(x=NOAA_all_rates$coeff_mean[k],lwr=NOAA_all_rates$coeff_lwr[k],upr=NOAA_all_rates$coeff_upr[k],n=nsim)
  i=ReefBudgetRandCalc(x=NOAA_all_rates$int_mean[k],lwr=NOAA_all_rates$int_lwr[k],upr=NOAA_all_rates$int_upr[k],n=nsim)
  r=runif(nsim,NOAA_all_rates$rugosity_m[k]-NOAA_all_rates$rugosity_sd[k],NOAA_all_rates$rugosity_m[k]+NOAA_all_rates$rugosity_sd[k])
  cf=runif(nsim,NOAA_all_rates$cf_lwr[k],NOAA_all_rates$cf_upr[k])
  b=microbioerosion_pacific/10
  n=100/s
  g=(n*cf*((c+b)*s*r+i))/10
  calc=round(quantile(g,probs=0.5),digits=2)
  calc_lower=round(quantile(g,probs=0.25),digits=2)
  calc_upper=round(quantile(g,probs=0.75),digits=2)
  NOAA_calc_rates=rbind(NOAA_calc_rates,data.frame(NAME,CODE,GROUP,TAXA,calc,calc_lower,calc_upper))
}

#find hard substrate that does not include CCA
NOAA_CRED_hard_substrate=NOAA_CRED_CoralNet_labelset[str_detect(NOAA_CRED_CoralNet_labelset$GROUP,"Hard Substrate"),]
NOAA_CRED_hard_nonCCA=NOAA_CRED_hard_substrate[!str_detect(NOAA_CRED_hard_substrate$NAME,"CCA"),]

#Import micro-bioerosion and macro-bioerosion data and apply to "rock" substrate labels
NOAA_CRED_bioerosion=cbind(NOAA_CRED_hard_nonCCA[,1:4],calc=bioerosion_pacific,calc_lower=bioerosion_pacific,calc_upper=bioerosion_pacific)

#assign NA values to all other labels
NOAA_CRED_other_label_rates=cbind(subset(NOAA_CRED_CoralNet_labelset,GROUP!="Hard coral"&GROUP!="Hard Substrate")[,1:4],calc=NA,calc_lower=NA,calc_upper=NA)

#select only necessary columns
NOAA_CRED_all=rbind(NOAA_calc_rates,NOAA_CRED_bioerosion,NOAA_CRED_other_label_rates)
NOAA_CRED_all$region="Indo-Pacific"
colnames(NOAA_CRED_all)=c("name","label","group","taxa","calc","calc_lower","calc_upper","region")
NOAA_CRED_all$label=str_remove(NOAA_CRED_all$label, "[*]")
calc_rates_indopacific=bind_cols('Region'=NOAA_CRED_all$region,'Name'=NOAA_CRED_all$name,'Mean'=NOAA_CRED_all$calc,'Lower bound'=NOAA_CRED_all$calc_lower,'Upper bound'=NOAA_CRED_all$calc_upper)

#import duplicate labels for NOAA CRED to general labels
colnames(duplicates)=c("Name","Duplicate")
calc_rates_indopacific_duplicates=merge(calc_rates_indopacific,duplicates,by="Name")
duplicate_labels_indopacific=calc_rates_indopacific_duplicates%>%select('Region','Duplicate','Mean','Lower bound','Upper bound')
duplicate_labels_indopacific=rename(duplicate_labels_indopacific,Name=Duplicate)

#generate genus level means from genus+morphology NOAA CRED labels
Acropora=calc_rates_indopacific[str_which(calc_rates_indopacific$Name,"Acropora"),]
Acropora_rates=bind_cols(Region="Indo-Pacific",Name="Acropora",'Mean'=round(mean(Acropora$'Mean'),digits=2),'Lower bound'=round(mean(Acropora$'Lower bound'),digits=2),'Upper bound'=round(mean(Acropora$'Upper bound'),digits=2))
Montipora=calc_rates_indopacific[str_which(calc_rates_indopacific$Name,"Montipora"),]
Montipora_rates=bind_cols(Region="Indo-Pacific",Name="Montipora",'Mean'=round(mean(Montipora$Mean),digits=2),'Lower bound'=round(mean(Montipora$'Lower bound'),digits=2),'Upper bound'=round(mean(Montipora$'Upper bound'),digits=2))
Porites=calc_rates_indopacific[str_which(calc_rates_indopacific$Name,"Porites"),]
Porites_rates=bind_cols(Region="Indo-Pacific",Name="Porites",'Mean'=round(mean(Porites$Mean),digits=2),'Lower bound'=round(mean(Porites$'Lower bound'),digits=2),'Upper bound'=round(mean(Porites$'Upper bound'),digits=2))

#combine NOAA CRED labels with duplicate labels and genus level labels
calc_rates_indopacific_all=bind_rows(calc_rates_indopacific,duplicate_labels_indopacific,Acropora_rates,Montipora_rates,Porites_rates) %>% arrange(Name)
calc_rates_indopacific_all=calc_rates_indopacific_all[complete.cases(calc_rates_indopacific_all),]

#########CoralNet Western Atlantic Calcification Rates##########
#Import Required Data
setwd("Data")
ReefBudget_Atlantic <- read_excel("Caribbean_carbonate_production_template_V2.1.xlsx", sheet = "Calcification Rates", skip = 9) #import Indo-Atlantic ReefBudget rates
Gonzalez_Barrios <- read_csv("Gonzalez-Barios_Alvarez-Filip_2019_Rugosity.csv") #Import Gonzalez-Barrios taxa-level estimates of colony structural complexity
CARICOMP <- read_delim("CARICOMP_Coral_Links.txt","\t", escape_double = FALSE, trim_ws = TRUE) #import CARICOMP data
CARICOMP_Species <- read_delim("CARICOMP_Species.txt","\t", escape_double = FALSE, trim_ws = TRUE) #import CARICOMP species names
CoralNet_Atlantic <- read_csv("CoralNet_labelset.csv") #Import CoralNet labels
microbioerosion=read_excel("Caribbean_carbonate_production_template_V2.1.xlsx", sheet = "Microbioerosion") #import Caribbean ReefBudget microbioerosion data
macrobioerosion=read_excel("Caribbean_carbonate_production_template_V2.1.xlsx", sheet = "Macrobioerosion", skip = 3) #import Caribbean ReefBudget macrobioerosion data
parrotbioerosion=read_excel("Caribbean_Parrotfish_template_V2.xlsx", sheet = "Equations") #import Caribbean ReefBudget Parrotfish bitemark data
setwd("..")

#Reorganize ReefBudget dataframe
colnames(ReefBudget_Atlantic)=c("code","name","morphology","ext_mean","ext_sd","dens_mean","dens_sd","cf","coeff_mean","coeff_lwr","coeff_upr","int_mean","int_lwr","int_upr","notes")
ReefBudget_Atlantic <- ReefBudget_Atlantic[!is.na(ReefBudget_Atlantic$coeff_mean),]

#Merge rates with taxa-specific rugosity index from Gonzalez-Barrios data
rugosity_rates= merge(ReefBudget_Atlantic,Gonzalez_Barrios,by="name",all.x=TRUE)

#find taxa with complete and missing rugosity data
complete_rugosity=rugosity_rates[!is.na(rugosity_rates$rugosity_mean),]
missing_rugosity=rugosity_rates[is.na(rugosity_rates$rugosity_mean),]

#replace missing SD data with zeros
complete_rugosity[!is.na(str_match(complete_rugosity$rugosity_sd,"nd")),]$rugosity_sd=0

#replace missing branching morphologies with mean±SD from Gonzalez-Barrios branching rugosities
missing_rugosity[!is.na(str_match(missing_rugosity$morphology.x,"bra")),]$rugosity_mean=mean(subset(Gonzalez_Barrios,morphology=="Massive")$rugosity_mean)
missing_rugosity[!is.na(str_match(missing_rugosity$morphology.x,"bra")),]$rugosity_sd=sd(subset(Gonzalez_Barrios,morphology=="Massive")$rugosity_mean)
#replace missing plating morphologies with mean±SD from Gonzalez-Barrios foliose rugosities
missing_rugosity[!is.na(str_match(missing_rugosity$morphology.x,"plat")),]$rugosity_mean=mean(subset(Gonzalez_Barrios,morphology=="Foliose")$rugosity_mean)
missing_rugosity[!is.na(str_match(missing_rugosity$morphology.x,"plat")),]$rugosity_sd=sd(subset(Gonzalez_Barrios,morphology=="Foliose")$rugosity_mean)
#replace missing massive and submassive morphologies with mean±SD from Gonzalez-Barrios massive rugosities
missing_rugosity[!is.na(str_match(missing_rugosity$morphology.x,"mass")),]$rugosity_mean=mean(subset(Gonzalez_Barrios,morphology=="Massive")$rugosity_mean)
missing_rugosity[!is.na(str_match(missing_rugosity$morphology.x,"mass")),]$rugosity_sd=sd(subset(Gonzalez_Barrios,morphology=="Massive")$rugosity_mean)
#replace missing encrusting morphologies with mean±SD from all Gonzalez-Barrios rugosities
missing_rugosity[is.na(missing_rugosity$rugosity_mean),]$rugosity_mean=mean(Gonzalez_Barrios$rugosity_mean)
missing_rugosity[is.na(missing_rugosity$rugosity_sd),]$rugosity_sd=sd(Gonzalez_Barrios$rugosity_mean)

#merge complete and missing rugosities with filled in rugosity data
rug_rates=rbind(complete_rugosity,missing_rugosity)
rug_rates$rugosity_sd=as.numeric(rug_rates$rugosity_sd)
rug_rates$genera=word(rug_rates$name, 1)

#load in mean colony sizes from the CARICOMP data set
colnames(CARICOMP)=c("station","date","transect","start_link","end_link","category","species")
CARICOMP$length_cm=with(CARICOMP,end_link-start_link) #manual states link size = 1 cm
CARICOMP$SPECIES=str_to_upper(CARICOMP$species, locale = "en")

#compute typical lengths for all encrusting CCA
mean_CCA_lengths=
  subset(CARICOMP,category=="EALG") %>%
  summarise(length_mean=quantile(length_cm, probs = 0.5),length_lwr=quantile(length_cm, probs = 0.025),length_upr=quantile(length_cm, probs = 0.975))

#compute typical lengths for all corals
mean_species_lengths=
  CARICOMP %>%
  group_by(SPECIES) %>%
  summarise(length_mean=quantile(length_cm, probs = 0.5),length_lwr=quantile(length_cm, probs = 0.025),length_upr=quantile(length_cm, probs = 0.975))

#read in species corresponding to CARICOMP codes
colnames(CARICOMP_Species)=c("label","species")
CARICOMP_Species$SPECIES=str_to_upper(CARICOMP_Species$label, locale = "en")

#Update taxa with more recent names
CARICOMP_Species$species=str_replace(CARICOMP_Species$species,"Montastrea annularis","Orbicella annularis")
CARICOMP_Species$species=str_replace(CARICOMP_Species$species,"Montastrea faveolata","Orbicella faveolata")
CARICOMP_Species$species=str_replace(CARICOMP_Species$species,"Montastrea franksi","Orbicella franksi")
CARICOMP_Species$species=str_replace(CARICOMP_Species$species,"Diploria strigosa","Pseudodiploria strigosa")
CARICOMP_Species$species=str_replace(CARICOMP_Species$species,"Diploria clivosa","Pseudodiploria clivosa")

#merge species lengths with their respective names
CARICOMP_length_species_raw=merge(mean_species_lengths,CARICOMP_Species,by="SPECIES")
CARICOMP_length_species=as_tibble(with(CARICOMP_length_species_raw,cbind(species,length_mean,length_lwr,length_upr)))
CARICOMP_length_species$name=str_remove(CARICOMP_length_species$species, " sp.")
CARICOMP_length_species$genera=word(CARICOMP_length_species$name, 1)

#compute genera means for all taxa
CARICOMP_length_genera = CARICOMP_length_species %>%
  group_by(genera) %>%
  summarize(length_mean=mean(as.numeric(length_mean)),length_lwr=mean(as.numeric(length_lwr)),length_upr=mean(as.numeric(length_upr)))

#fill in missing species data with genera means
CARICOMP_species_rates=merge(rug_rates,CARICOMP_length_species,by="name",all.x=TRUE)
CARICOMP_species_rates_na=CARICOMP_species_rates[is.na(CARICOMP_species_rates$length_mean),]
CARICOMP_species_rates_na$genera=CARICOMP_species_rates_na$genera.x
CARICOMP_species_rates_na=CARICOMP_species_rates_na[,c(1:18,25)]
CARICOMP_genera_rates=merge(CARICOMP_species_rates_na,CARICOMP_length_genera,by="genera",all.x=TRUE)
CARICOMP_genera_rates2=CARICOMP_genera_rates[!is.na(CARICOMP_genera_rates$length_mean),2:22]

#find missing genera lengths
CARICOMP_genera_rates_na=CARICOMP_genera_rates[is.na(CARICOMP_genera_rates$length_mean),2:19]

#finding missing species lengths after filling in genera lengths
CARICOMP_species_rates2=CARICOMP_species_rates[!is.na(CARICOMP_species_rates$length_mean),]
CARICOMP_species_rates2=CARICOMP_species_rates2[,c(1:18,21:23)]

#find typical lengths for each morphology to substitute in for missing species and genera data
CARICOMP_length_morphology = CARICOMP_species_rates2 %>%
  group_by(morphology.x) %>%
  summarize(length_mean=mean(as.numeric(length_mean)),length_lwr=mean(as.numeric(length_lwr)),length_upr=mean(as.numeric(length_upr)))

#merge missing data by morphologies
CARICOMP_morphology=merge(CARICOMP_genera_rates_na,CARICOMP_length_morphology,by="morphology.x",all.x=TRUE)
CARICOMP_morphology[7,19:21]=CARICOMP_length_morphology[5,2:4] #substitute massive morphology data for Dendrogyra
CARICOMP_morphology[3,19:21]=mean_CCA_lengths
CARICOMP_morphology[4,19:21]=mean_CCA_lengths
CARICOMP_morphology[5,19:21]=mean_CCA_lengths
CARICOMP_morphology[6,19:21]=mean_CCA_lengths

#merge all length data and clean up to merge with rugosity and calcification rate data
CARICOMP_rates=rbind(CARICOMP_species_rates2,CARICOMP_genera_rates2,CARICOMP_morphology)
CARICOMP_rates=CARICOMP_rates[order(CARICOMP_rates$name),]
CARICOMP_rates$length_mean=as.numeric(CARICOMP_rates$length_mean)
CARICOMP_rates$length_lwr=as.numeric(CARICOMP_rates$length_lwr)
CARICOMP_rates$length_upr=as.numeric(CARICOMP_rates$length_upr)

#apply open space conversion factor to branching taxa
CARICOMP_rates$genera=word(CARICOMP_rates$name, 1)
CARICOMP_rates$cf_lwr=ifelse(CARICOMP_rates$morphology.x=="branching" | CARICOMP_rates$morphology.y=="Digitate" & !is.na(CARICOMP_rates$morphology.y=="Digitate"), 0.168 , 1 ) #Doszpot et al (2019) minimum proportion of 3D space occupied in branching colonies
CARICOMP_rates$cf_upr=ifelse(CARICOMP_rates$morphology.x=="branching" | CARICOMP_rates$morphology.y=="Digitate" & !is.na(CARICOMP_rates$morphology.y=="Digitate"), 0.598 , 1 ) #Doszpot et al (2019) maximum proportion of 3D space occupied in branching colonies

#apply open space conversion factor specific to Acropora cervicornis
CARICOMP_rates[CARICOMP_rates$code=="ACC",]$cf_lwr=0.2 #Million et al. (2021) 25th percentile 12 month V_coral/V_convex-hull
CARICOMP_rates[CARICOMP_rates$code=="ACC",]$cf_upr=0.35 #Million et al. (2021) 75th percentile 12 month V_coral/V_convex-hull

#create genera means
genera_rates = CARICOMP_rates %>%
  group_by(genera) %>%
  summarize(coeff_mean=mean(coeff_mean),coeff_lwr=mean(coeff_lwr),coeff_upr=mean(coeff_upr),int_mean=mean(int_mean),int_lwr=mean(int_lwr),int_upr=mean(int_upr),rugosity_mean=mean(rugosity_mean),rugosity_sd=mean(rugosity_sd),length_mean=mean(length_mean),length_lwr=mean(length_lwr),length_upr=mean(length_upr),cf_lwr=mean(cf_lwr),cf_upr=mean(cf_upr))

#remove macroalgae genera because this only applies to Macroalgae with CCA
genera_rates=genera_rates[is.na(str_match(genera_rates$genera,"Macroalgae")),]
#remove Peysonellid genera to avoid duplicates
genera_rates=genera_rates[is.na(str_match(genera_rates$genera,"Peysonellid")),]

#merge genera means and taxa-specific rates
rug_means=as_tibble(with(CARICOMP_rates,cbind(name,coeff_mean,coeff_lwr,coeff_upr,int_mean,int_lwr,int_upr,rugosity_mean,rugosity_sd,length_mean,length_upr,length_lwr,cf_lwr,cf_upr)))
genera_means=as_tibble(with(genera_rates,cbind(name=genera,coeff_mean,coeff_lwr,coeff_upr,int_mean,int_lwr,int_upr,rugosity_mean,rugosity_sd,length_mean,length_upr,length_lwr,cf_lwr,cf_upr)))
rug_all=as_tibble(rbind(genera_means,rug_means))
rug_all=transform(rug_all, coeff_mean = as.numeric(coeff_mean),coeff_lwr = as.numeric(coeff_lwr),coeff_upr = as.numeric(coeff_upr),int_mean = as.numeric(int_mean),int_lwr = as.numeric(int_lwr),int_upr = as.numeric(int_upr),rugosity_mean=as.numeric(rugosity_mean),rugosity_sd=as.numeric(rugosity_sd),length_mean=as.numeric(length_mean),length_lwr=as.numeric(length_lwr),length_upr=as.numeric(length_upr),cf_lwr=as.numeric(cf_lwr),cf_upr=as.numeric(cf_upr))

#remove erroneous "other" from genera data
rug_all=rug_all[is.na(str_match(rug_all$name,"Other")),]

#include microbioerosion rate based on initial logical
microbioerosion_atlantic=ifelse(include_bioerosion==TRUE,-1,0)*as.numeric(microbioerosion[2,5])

#determines mean, lower, and upper rate in kg m^-2 yr^-1 for 100% cover of the respective coral using CARICOMP colony sizes and species-specific rugosity from Gonzalez-Barrios
scaled_calcification_rates=NULL
for (k in 1:nrow(rug_all))
{
  nsim=10000
  name=rug_all$name[k]
  s=ReefBudgetRandCalc(x=rug_all$length_mean[k],lwr=rug_all$length_lwr[k],upr=rug_all$length_upr[k],n=nsim)
  c=ReefBudgetRandCalc(x=rug_all$coeff_mean[k],lwr=rug_all$coeff_lwr[k],upr=rug_all$coeff_upr[k],n=nsim)
  i=ReefBudgetRandCalc(x=rug_all$int_mean[k],lwr=rug_all$int_lwr[k],upr=rug_all$int_upr[k],n=nsim)
  r=runif(nsim,rug_all$rugosity_mean[k]-rug_all$rugosity_sd[k],rug_all$rugosity_mean[k]+rug_all$rugosity_sd[k])
  cf=runif(nsim,rug_all$cf_lwr[k],rug_all$cf_upr[k])
  b=microbioerosion_atlantic/10
  n=100/(s/r)
  g=(n*cf*((c+b)*s+i))/10
  calc=round(quantile(g,probs=0.5),digits=2)
  calc_lower=round(quantile(g,probs=0.25),digits=2)
  calc_upper=round(quantile(g,probs=0.75),digits=2)
  scaled_calcification_rates=rbind(scaled_calcification_rates,data.frame(name,calc,calc_lower,calc_upper))
}

#Modify CoralNet label columns
colnames(CoralNet_Atlantic)=c("name","functional_group","label")

#merge calcification rates with CoralNet Labels
CoralNet_Atlantic_Rates = merge(scaled_calcification_rates,CoralNet_Atlantic,by="name",all=TRUE)

#Correct mismatched names (including typos) with CoralNet to ensure correct mapping of calcification rate data
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Crustose coralline algae")),]$name="CCA (crustose coralline algae)"
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Dichocoenia stokesii")),]$name="Dichocoenia stokesi"
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"^Hard$")),]$name="Hard coral"
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Hard coral \\(branched\\)")),]$name="Hard Coral (branching)"
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Hard coral \\(encrusting\\)")),]$name="Hard Coral (encrusting)"
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Hard coral \\(massive\\)")),]$name="Hard Coral (massive)"
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Hard coral \\(plate/foliose\\)")),]$name="Hard Coral (foliose)"
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Mycetophyllia danae")),]$name="Mycetophyllia danaana"
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Mycetophyllia lamarckiana")),]$name="Mycetophyllia lamarckana"
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Peysonellid")),]$name="Macroalgae: Laminate: red: peysonnelia"

#Import microbioerosion data and apply to "rock" substrate labels
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Rock")),]$calc=microbioerosion_atlantic
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Rock")),]$calc_lower=microbioerosion_atlantic
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Rock")),]$calc_upper=microbioerosion_atlantic

#Replace "Rock Crustose Coralline Algae" label with "Crustose Coralline Algae" rate
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Rock Crustose Coralline Algae")),]$calc=scaled_calcification_rates[!is.na(str_match(scaled_calcification_rates$name,"CCA (crustose coralline algae)")),]$calc
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Rock Crustose Coralline Algae")),]$calc_lower=scaled_calcification_rates[!is.na(str_match(scaled_calcification_rates$name,"CCA (crustose coralline algae)")),]$calc_lower
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Rock Crustose Coralline Algae")),]$calc_upper=scaled_calcification_rates[!is.na(str_match(scaled_calcification_rates$name,"CCA (crustose coralline algae)")),]$calc_upper

#Import macrobioerosion data and apply to boring sponge substrate labels
macrobiorates=macrobioerosion[1:8,c(3,7)]
colnames(macrobiorates)=c("name","rate")
macrobiorates$rate=round(ifelse(include_bioerosion==TRUE,-1,0)*as.numeric(macrobiorates$rate),digits=2)
cliona_delitrix_rate=round(macrobiorates[!is.na(str_match(macrobiorates$name,"delitrix")),]$rate,digits=2)
cliona_mean=round(mean(macrobiorates[!is.na(str_match(macrobiorates$name,"C.")),]$rate),digits=2)
cliona_sd=round(sd(macrobiorates[!is.na(str_match(macrobiorates$name,"C.")),]$rate),digits=2)
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Clio")),]$calc=cliona_mean
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Clio")),]$calc_lower=(cliona_mean-cliona_sd)
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Clio")),]$calc_upper=(cliona_mean+cliona_sd)
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Cliona delitrix")),]$calc=cliona_delitrix_rate
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Cliona delitrix")),]$calc_lower=cliona_delitrix_rate
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$name,"Cliona delitrix")),]$calc_upper=cliona_delitrix_rate

#Import CaCO3 removed per bite and apply to BiteScar label
parrot_bites=as_tibble(parrotbioerosion[45:54,2:10])
all_parrot_bites=as.numeric(stack(parrot_bites[1:ncol(parrot_bites)])$values)
parrot_mean = round(ifelse(include_bioerosion==TRUE,-1,0)*mean(all_parrot_bites,na.rm=TRUE),digits=2)
parrot_sd = round(ifelse(include_bioerosion==TRUE,-1,0)*sd(all_parrot_bites,na.rm=TRUE),digits=2)
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$label,"BiteScar")),]$calc=parrot_mean
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$label,"BiteScar")),]$calc_lower=parrot_mean-parrot_sd
CoralNet_Atlantic_Rates[!is.na(str_match(CoralNet_Atlantic_Rates$label,"BiteScar")),]$calc_upper=parrot_mean+parrot_sd

#add atlantic region label
CoralNet_Atlantic_Rates$region="Atlantic"
#select only names that match CoralNet, remove extra columns, and drop missing calcification rates
CoralNet_Atlantic_Rates2=merge(CoralNet_Atlantic_Rates,CoralNet_Atlantic,by="name",all.y=TRUE) %>% 
  select('Region'=region,'Name'=name,'Mean'=calc,'Lower bound'=calc_lower,'Upper bound'=calc_upper) %>% 
  drop_na()
#order calcification rates by name
calc_rates_atlantic=CoralNet_Atlantic_Rates2[order(CoralNet_Atlantic_Rates2$Name),]

#########Merge Indo-Pacific and Western Atlantic Calcification Rate Tables##########
calc_rates=bind_rows(calc_rates_indopacific_all,calc_rates_atlantic) %>% drop_na()
write.csv(calc_rates,"CoralNet_Calcification_Bioerosion_Rates_v1.csv",row.names=FALSE)
