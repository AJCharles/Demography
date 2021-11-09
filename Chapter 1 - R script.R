#***************PLEASE READ***************

#This script wrangles the data, producing survival and mortality plots and a simple additive and interactive cox mixed effect models. Originally written by Mirre, and worked on my most since (mainly by Gracie), I've stripped out unnecessary sections whilst adding my own specific to my data. 

# Notable changes / improvements; 
# My code below will assign the diet condition to the fly respective of fly actual age so that the diet at switch date is still the previous, but the diet at the first recorded demography count is the new one (switch accurate). It also creates a variable 'relative to switch age' respective to differences in fly DOBs so that the data may be wrangled appropriately for modelling. 


library(reshape)
library(lme4)
library(splitstackshape)
library(PropCIs)
library(ggplot2)
library(coxme)
library(GGally) #for ggsurv
library(dplyr)
library(gridExtra)
library(ggfortify)
library(multcomp)
library(scales)
library(tidyr)
library(survminer)
library(knitr)
library(kableExtra)
library(magick)
library(broom)
library(ehahelper)

  #To clear workspace:
rm(list = ls())

#setwd and load data
wd = "C:/Users/alexc/Documents/THESIS/Chapter 1 - Demography"

setwd(wd)
filetoread="Corrected-FINAL-Modelling.csv"
filetoread="Corrected-FINAL-Plotting.csv"

data<-read.csv(filetoread)
columns=4
#Setting 4 column-wide format (for each cage there is dead, at risk, age, treat)


##------------------- Drop unwanted cages by filtering values - Currently only want first 40!! --------------------#
cages=dim(data)[2]/columns-20  # drop untested cages (last 20) || What we decided on for now.
#cages=dim(data)[2]/columns-40  # drop untested and transcriptomics only (last 40)
#cages=dim(data)[2]/columns     # don't drop any cages and assess all demographic data available.

# Conveniently, I only want to drop the end-most cages from the df (41-60) as I don't analyse them for proteomics OR transcriptomics. Therefore if I simply index the first 40 only (so 41-60 are out of range) then they won't even be integrated into datalong; no faffing with individual df's downstream as they don't make the cut from the raw data here.


#How many cages is defined by how wide data sheet is - number of columns divided by 4
entries=dim(data)[1]
#The max number of rows, how long the longest entry for a specific cage is - ie. the longest a cage survived
#Dimensions of data of the longest

#Section below is to generate long file matrix
tofill=matrix(nrow=cages*dim(data)[1],ncol=columns)
#Generate an empty matrix from a wide format to long format. Amount of cages by amount of entries

#Complicated while loops for the matrix
k=1
while(k<(cages+1))
{
  tofill[(1+((k-1)*entries)):(k*entries),]=	as.matrix(data[,(1+((k-1)*columns)):(k*columns)])	
  k=k+1
}
datalong=as.data.frame(tofill)

#Matrix assignment
names(datalong)=names(data[1:columns])
names(datalong)[1]="dead"
datalong3=as.data.frame(datalong)
datalong3$dead=as.numeric(as.character(datalong$dead))
datalong3$age=as.numeric(as.character(datalong$age))
datalong3$at.risk=as.numeric(as.character(datalong$at.risk))
datalongf=datalong3
datalong=datalongf
#Cage numbers need to be kept as text


##-------------------CENSORING
#Clean n/a for statistical analysis
datalong=datalong[-which(is.na(datalong$dead)==T),] 
datalong$censor=datalong$at.risk-datalong$at.risk[2:(length(datalong$at.risk)+1)]-datalong$dead
datalong$censor[which(datalong$censor<0 | is.na(datalong$censor)==T)]=0

#Censored events like food deaths or alive need to be calculated for the output data for cox proportional hazard stats
#Subtract two matrices from one another to get a value for how many individuals are censored at each date

#Censor all <0
datalong=cSplit(datalong, 'treat', sep="@", type.convert=FALSE)
datalong=cSplit(datalong, 'treat_1', sep="_", type.convert=FALSE)
datalong=cSplit(datalong, 'treat_2', sep="_", type.convert=FALSE)
#This is to separate your variables such as tray, cage number, treatment that are in one cell of the 'treat' column in the csv file
#Two variables, treat 1 and treat 2. Treat 1 is genotype


##-------------------PREPARING DATA FOR COX MODELS
#Cox proportional likes individual based data, one row per death. So need to expand data using frequency of death as number of
#repeats and assigning 1 to each row for death. The following is to get data in order, this is to make repeats of the rows for 
#dead and censored flies
dataset<-as.data.frame(cbind(rep(datalong$treat_1_1,datalong$dead),
                             rep(datalong$treat_1_2,datalong$dead),
                             rep(datalong$age,datalong$dead),
                             rep(datalong$at.risk,datalong$dead),
                             rep(datalong$treat_2_1, datalong$dead),
                             rep(datalong$treat_2_2, datalong$dead)))
names(dataset)=c("Genotype","Cage","Age", "At.Risk", "Transfer", "Treatment") #Column names should be correct here, but if you alter the script to output different column names you will need to change the names in "" here

dataset$Age<-as.numeric(as.character(dataset$Age))
dataset$Dead=1

#rep - repeat - gives you a repeat of the data sheet for as many times as there are censors
datasetcensor=as.data.frame(cbind(rep(datalong$treat_1_1,datalong$censor),
                                  rep(datalong$treat_1_2,datalong$censor),
                                  rep(datalong$age,datalong$censor),
                                  rep(datalong$at.risk,datalong$censor),
                                  rep(datalong$treat_2_1, datalong$censor),
                                  rep(datalong$treat_2_2, datalong$censor)))
names(datasetcensor)=c("Genotype","Cage","Age", "At.Risk", "Transfer", "Treatment") #Column names should be correct here, but if you alter the
#script to output different column names you will need to change the names in "" here
datasetcensor$Age<-as.numeric(as.character(datasetcensor$Age))
datasetcensor$Dead=0 #censor = 0

#Reason for this is for the cox models that require 1 row for each observation. This isn't necessary just for plotting, but it is for the models.

#Combine dataset met dead and censor
dataset=rbind(dataset,datasetcensor) ##------------##-----------##------------##------------ END MODELLING
unique(dataset$Cage) # check that this final dataset only contains cage entries you selected for at the start.. essential! 



##-------------------PLOTTING TO CHECK FOR ERRORS

hist(datalong$age) #Check for normal distrib AND that not got ridiculous numbers for age
plot(datalong$dead) #Check not ridiculous number died on a day
plot(datalong$censor) #X axis is each cage but random number basically. Check how many censored per cage
hist(dataset$Age) #Check not any ridiculous ages
hist(dataset$Dead) #Check they are either dead or alive


##-------------------SURVIVAL PLOT ####
finaldata <- dataset

#Change what's after the ~ to change the survival plot produced
kaplan<-survfit(Surv(Age,Dead)~Treatment,data=finaldata) #kaplan mire plot, survfit is function you'd use to plot a kaplan mire
#Explain (number of dead by age) then explain all of that by treatment

plot.data = fortify(kaplan)
#Uses ggfortify turning the kaplan model into a dataframe

## CHANGE - plot.data to diet groups, and strata = diet. ####

# change 'strata' to be 'Treatment
names(plot.data)[names(plot.data)=='strata'] <- "Treatment"

# check current treatment names, substitute with unique descriptors & fix any changes if necessary.
plot.data$Treatment <- gsub("2-8", "Switch to rich", plot.data$Treatment)
plot.data$Treatment <- gsub("8-2", "Switch to DR", plot.data$Treatment) # run the hypens first
plot.data$Treatment <- gsub("2", "Continuous DR", plot.data$Treatment)
plot.data$Treatment <- gsub("8", "Continuous rich", plot.data$Treatment)
plot.data$Treatment <- as.factor(plot.data$Treatment)
str(plot.data)

kaplan_modelling_fortified <- plot.data # keep separate in the environment for use later... 
## hangon, the kaplan pdf is giving actual median ci's here. TF?!? page 12...
kaplan ##
       # well then... this is what I wanted from the start. 
kaplan_modelling_raw <- kaplan



## Specify manual aesthetics for survival plots ####
cols <-     c('Continuous DR' = 'black', 
              'Continuous rich' = 'red', 
              'Switch to DR' = 'goldenrod', 
              'Switch to rich' = 'cornflowerblue')

linetype <- c('Continuous DR' = 'solid', 
              'Continuous rich' = 'solid', 
              'Switch to DR' = 'dashed', 
              'Switch to rich' = 'dashed')

## Plot survival curve####
#If you have >13 treatments, delete 'linetype =  strata' because it is limiting. Use colours only
QC_Survival_scoring_intervals <- ggplot(plot.data, aes(time, surv)) +
  geom_point(size = 2, aes(colour = Treatment)) +
  scale_color_manual(values=cols) +
  scale_linetype_manual(values=linetype) +
  labs(x="Age (Days)", y="Survival") +
  theme_classic() +
  theme(axis.title.x = element_text(face="bold", size=22),
        axis.text.x = element_text(face= 'bold', size=13, colour = 'gray30')) +
  theme(axis.title.y = element_text(face="bold", size=22),
        axis.text.y = element_text(face= 'bold', size=13, colour = 'gray30')) +
  theme(axis.line = element_line(colour = 'black', size = 1)) +
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  theme(legend.text=element_text(size=13)) +
  theme(legend.key.width = unit(2.5, "line")) +
  theme(legend.position = c(0.82, 0.90)) +    # TO CHANGE LEGEND POSITION
  theme(legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank()) +
  scale_y_continuous(expand = c(0.002, 0.002), #0.002 for both
                     breaks = c(0, .25, .5, .75, 1),
                     labels = c('0', '.25', '.5', '.75', '1'),
                     limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, max(plot.data$time)+10, by = 10),
                     labels = seq(0, max(plot.data$time)+10, by = 10),
                     expand = c(0.005, 0.005),
                     limits = c(0, 90)) +
  theme(aspect.ratio = 1) +
  geom_vline(xintercept = seq(1,89,2), linetype = "dashed", colour = "dark grey") 


## Plot survival QC ####
QC_Survival_scoring_intervals # need to keep colour codes for 2% and 8% consistent! 

SurvivalCurve <- ggplot(plot.data, aes(time, surv)) +
  geom_step(size = 1.4, aes(colour = Treatment, linetype = Treatment)) +
  scale_color_manual(values=cols) +
  scale_linetype_manual(values=linetype) +
  labs(x="Age (Days)", y="Survival", size = 10) +
  theme_classic() +
  theme(axis.title.x = element_text(face="bold", size=22),
        axis.text.x = element_text(face= 'bold', size=17, colour = 'gray30')) +
  theme(axis.title.y = element_text(face="bold", size=22),
        axis.text.y = element_text(face= 'bold', size=17, colour = 'gray30')) +
  theme(axis.line = element_line(colour = 'black', size = 1.5)) +
  theme(axis.ticks = element_line(colour = "black", size = 1.75)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  theme(legend.text=element_text(size=17)) +
  theme(legend.title = element_text(size=20)) +
  theme(legend.key.width = unit(2.5, "line")) +
  theme(legend.position = c(0.76, 0.90)) +    # TO CHANGE LEGEND POSITION
  theme(legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank()) +
  scale_y_continuous(expand = c(0.002, 0.002), #0.002 for both
                     breaks = c(0, .25, .5, .75, 1),
                     labels = c('0', '.25', '.5', '.75', '1'),
                     limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, max(plot.data$time)+10, by = 10),
                     labels = seq(0, max(plot.data$time)+10, by = 10),
                     expand = c(0.005, 0.005),
                     limits = c(0, 90)) +
  theme(aspect.ratio = 1) +
  geom_vline(xintercept = 21, linetype = "dotted", colour = "gray30", size = 1.25)  ## OPTIONAL intercept LINE


SurvivalCurve # need to keep colour codes for 2% and 8% consistent! 

##-------------------MORTALITY PLOT CALCULATIONS ####
datalong <- dplyr::rename(datalong, Transfer = treat_2_1, Genotype = treat_1_1, 
                          Cage = treat_1_2, Treatment = treat_2_2, Dead = dead,
                          At.Risk = at.risk, 
                          Age = age, Censor = censor)

dataplot=aggregate(cbind(At.Risk,Dead)~Age+Treatment,data=datalong,FUN=sum) #cbind = columnbind. Gets (for whole geno) the total no. dead per age per treatment. Forgets about cage.
sapply(dataplot, class) #Shows what variable they are in the console
dataplot$mx=dataplot$Dead/dataplot$At.Risk  #Dividing the no. dead by no. at risk in same row to estimate mortality
dataplot$mx=dataplot$mx/2 #Divide mortality by 2 to get every day mortality (currently we have data for every other day because of scoring)

##-------------------MORTALITY PLOT AGGREGATION ####
dataplot[which(dataplot$mx ==0), 5] <- NA #Make sure second number is the correct number of columns in dataplot (should be 5)
dataplot$mx[is.nan(dataplot$mx)] <- NA #mx is the mortality rate column. Want to say that anytime there is a mortality rate of 0, 
#make it NA so it doesn't plot O

str(dataplot) #Check looks correct in console

sapply(dataplot, class) #Check what variables they are
dataplot$Treatment=as.factor(as.factor(dataplot$Treatment)) #Turns treatment variable into a factor


##-------------------MORTALITY PLOT || Changed to make it look like the survival... why it wasn't already ?? ####
#If you have >13 treatments, delete 'linetype =  Treatment' because it is limiting. Use colours only

bk_dataplot <- dataplot
dataplot <- bk_dataplot

unique(dataplot$Treatment)

dataplot$Treatment <- gsub("2-8", "Switch to rich", dataplot$Treatment)
dataplot$Treatment <- gsub("8-2", "Switch to DR", dataplot$Treatment)
dataplot$Treatment <- gsub("2"    , "Continuous DR", dataplot$Treatment)
dataplot$Treatment <- gsub("8"    , "Continuous rich", dataplot$Treatment)
unique(dataplot$Treatment)
dataplot$Treatment <- as.factor(dataplot$Treatment)
str(dataplot)


## Plot an intitial QC dotplot to show all values align with scoring days. ####
QC_Mortality_Scoring_intervals <- ggplot(dataplot, aes(x=Age, y= log(mx))) +
  geom_point(size = 2, aes(colour = Treatment))+
  scale_color_manual(values=cols) +
  theme_classic() + # Optional extras
  labs(x="Age (Days)") +
  theme(aspect.ratio = 1) +
  theme(axis.title.x = element_text(face="bold", size=22),
        axis.text.x = element_text(face= 'bold', size=17, colour = 'gray30')) +
  theme(axis.title.y = element_text(face="bold", size=22),
        axis.text.y = element_text(face= 'bold', size=17, colour = 'gray30')) +
  theme(axis.line = element_line(colour = 'black', size = 1)) +
  theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  theme(legend.text=element_text(size=13)) +
  theme(legend.key.width = unit(2.5, "line")) +
  theme(legend.position = c(0.82, 0.10)) +     # Change legend position
  theme(legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank()) +
  scale_y_continuous(expand = c(0.002, 0.002))+
  scale_x_continuous(breaks = seq(0, max(plot.data$time)+10, by = 10),   
                     labels = seq(0, max(plot.data$time)+10, by = 10),
                     expand = c(0.005, 0.005),
                     limits = c(0, 90))  +
  geom_vline(xintercept = seq(1,89,2), linetype = "dashed", colour = "dark grey") +
  guides(color=guide_legend(override.aes=list(fill=NA)))

QC_Mortality_Scoring_intervals # ALL points should aline with the abline grid. This is a quick way to confirm that the data was correctly processed prior to use for plotting.

## Plot actual mortality plot ####

MortalityPlot <- ggplot(subset(dataplot, Treatment=="Continuous DR" | Treatment=="Continuous rich"), 
                        aes(x=Age, y= log(mx))) +
  geom_line(show.legend = FALSE, size = 1.4, aes(colour = Treatment, linetype = Treatment))+
  geom_point(data = subset(dataplot, Treatment=="Switch to rich" | Treatment=="Switch to DR"), 
             aes(colour = Treatment), size=2) +
  scale_color_manual(values=cols) +
  scale_linetype_manual(values=linetype) +
  theme_classic() + # Optional extras
  labs(x="Age (Days)") +
  theme(aspect.ratio = 1) +
  theme(axis.title.x = element_text(face="bold", size=22),
        axis.text.x = element_text(face= 'bold', size=17, colour = 'gray30')) +
  theme(axis.title.y = element_text(face="bold", size=22),
        axis.text.y = element_text(face= 'bold', size=17, colour = 'gray30')) +
  theme(axis.line = element_line(colour = 'black', size = 1.5)) +
  theme(axis.ticks = element_line(colour = "black", size = 1.75)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  theme(legend.text=element_text(size=17)) +
  theme(legend.spacing.x = unit(-.2, "cm")) + # fix text spaced from colours..
  theme(legend.title = element_text(size=20)) +
  theme(legend.key.width = unit(2.5, "line")) +
  theme(legend.position = c(0.82, 0.135)) +     # Change legend position
  theme(legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank()) +
  scale_y_continuous(expand = c(0.002, 0.002))+
  scale_x_continuous(breaks = seq(0, max(plot.data$time)+10, by = 10),   
                     labels = seq(0, max(plot.data$time)+10, by = 10),
                     expand = c(0.005, 0.005),
                     limits = c(0, 90))  +
  geom_vline(xintercept = 21, linetype = "dotted", colour = "gray30", size = 1.25) 

MortalityPlot







##--------End of plotting----- Statistical analysis should be conducted on the non-DOB-adjusted csv ####


##------------------ADD IN DIET CONDITION FOR ALL FLIES ####

# The current 'age' of the fly must match with the diet used previously due to the two-day sampling window. i.e. the dietary switch (ON the 16th May) is still COUNTING the demography of the previous diet, therefore the first ACTUAL SWITCH DAY (in the demography) is the 18th May. 
# Split the dataframe apart by cages which were made on the same days and perform filtering to generate accurate diets. Once conducted for all three days (26th, 27th & 28th) merge them back to create the FINAL modelling frame.

# Cages  1:20 = 26th April    ||    Age at switch = 20   ||      Age at initial demography count =  22
# Cages 21:44 = 27th April    ||    Age at switch = 19   ||      Age at initial demography count =  21
# Cages 45:60 = 28th April    ||    Age at switch = 18   ||      Age at initial demography count =  20

# Data wrangling:
Split_Data <- dataset  
str(Split_Data)
## Split_Data$Cage <- as.numeric(Split_Data$Cage)## NOT FINE!!! 
Split_Data$Treatment <- gsub("2-8", "SU", Split_Data$Treatment)
Split_Data$Treatment <- gsub("8-2", "SD", Split_Data$Treatment) # run the hypens first 
Split_Data$Treatment <- gsub("2", "2%", Split_Data$Treatment)
Split_Data$Treatment <- gsub("8", "8%", Split_Data$Treatment)
Split_Data$Treatment <- as.factor(Split_Data$Treatment)


library(tidyverse)
# workaround for the numeric conversion issue. For some reason as.numeric changed some of the data.... 
# possible that as.numeric(as.character(x)) fixes this... something to bring up?? Leave as is for now.  
tmp_col <- as.data.frame(Split_Data[,2]); names(tmp_col) <- "Cage"
tmp_col <- as.data.frame(as.numeric(paste0(tmp_col$Cage)));names(tmp_col) <- "Cage"
str(tmp_col)
Split_Data <- Split_Data[,-2]
Split_Data <- cbind(tmp_col,Split_Data)
str(Split_Data)


## missing my diet treatment x diet df... 

##### Still need to run for Treatment*Diet model - Filter Associated a new diet to each based on fly age (as specified above).  ####
#####    Also associate a 'relative to switch age, i.e. relswitchage' so that I can correctly model.

# I've done this by first further separating each split data section, making changes (according to SwitchAge) then stitch together. I'm sure there could be a better way to do this, but at least this way I'm confident that the diets are correct for the actual changes in demography. Plus coding it is a lot less effor than creating this in excel by hand... 

# Split Cages into relevant sections (DOBs) || Cage number needs to be numberic! 
Split_Start <- subset(Split_Data,Cage<=20&Cage>=1)    # Filtered based on DOB above | dropped cages shouldn't matter
Split_Middle <- subset(Split_Data,Cage<=44&Cage>=21)  # Filtered based on DOB above

# For the flies born on the 26th - i.e. Split_Start
Start_switchAge=22
Split_Start$Diet <- NA # create empty column.
Two_cont <- Split_Start %>% filter(Treatment=="2%") ; Two_cont$Diet <- "2%"
Eight_cont <- Split_Start %>% filter(Treatment=="8%") ; Eight_cont$Diet <- "8%"
SU_start <- Split_Start %>% filter(Age<Start_switchAge, Treatment=="SU") ; SU_start$Diet <- "2%" #BEFORE TREATMENT = DIET
SD_start <- Split_Start %>% filter(Age<Start_switchAge, Treatment=="SD") ; SD_start$Diet <- "8%"
SD_after <- Split_Start %>% filter(Age>=Start_switchAge, Treatment=="SD") ; SD_after$Diet <- "2%" # >= so first day
SU_after <- Split_Start %>% filter(Age>=Start_switchAge, Treatment=="SU") ; SU_after$Diet <- "8%"
Split_Start <- rbind(Two_cont,Eight_cont,SU_start,SD_start,SD_after,SU_after)

# Creat relative switching ages post-switch up until last fly age... 
Split_Start$relSwitchAge <- NA # create empty column in advance
start_age_range <- range(Split_Start$Age) # calculate age range
start_age_ps_seq <- seq(Start_switchAge, start_age_range[2], by=2) # calculate post switch age range
start_relswitchage_ps_seq <- start_age_ps_seq-Start_switchAge+2 # Add 2 as first demo record is day 2..
Split_Start$relSwitchAge[Split_Start$Age<Start_switchAge] <- 0 # Set all 'ages' before the switch as 0 (for group comp)

### while loop to index above so age = relative age from switch to subset for modelling later. 
p=0
while(p<(length(start_relswitchage_ps_seq))) {
  Split_Start$relSwitchAge[Split_Start$Age==start_age_ps_seq[1+p]] <- start_relswitchage_ps_seq[1+p]
  p=p+1
}



# For the flies born on the 27th - i.e. Split_Middle
Middle_switchAge=21
Split_Middle$Diet <- NA # create empty column.
Two_cont <- Split_Middle %>% filter(Treatment=="2%") ; Two_cont$Diet <- "2%"
Eight_cont <- Split_Middle %>% filter(Treatment=="8%") ; Eight_cont$Diet <- "8%"
SU_start <- Split_Middle %>% filter(Age<Middle_switchAge, Treatment=="SU") ; SU_start$Diet <- "2%" 
SD_start <- Split_Middle %>% filter(Age<Middle_switchAge, Treatment=="SD") ; SD_start$Diet <- "8%"
SD_after <- Split_Middle %>% filter(Age>=Middle_switchAge, Treatment=="SD") ; SD_after$Diet <- "2%" 
SU_after <- Split_Middle %>% filter(Age>=Middle_switchAge, Treatment=="SU") ; SU_after$Diet <- "8%"
Split_Middle <- rbind(Two_cont,Eight_cont,SU_start,SD_start,SD_after,SU_after)

# Creat relative switching ages post-switch up until last fly age... 
Split_Middle$relSwitchAge <- NA # create empty column in advance
Middle_age_range <- range(Split_Middle$Age) # calculate age range
Middle_age_ps_seq <- seq(Middle_switchAge, Middle_age_range[2], by=2) # calculate post switch age range
Middle_relswitchage_ps_seq <- Middle_age_ps_seq-Middle_switchAge+2 # Add 2 as first demo record is day 2..
Split_Middle$relSwitchAge[Split_Middle$Age<Middle_switchAge] <- 0 #Set all 'ages' before the switch as 0 (for group comp)

### while loop to index above so age = relative age from switch to subset for modelling later. 
p=0
while(p<(length(Middle_relswitchage_ps_seq))) {
  Split_Middle$relSwitchAge[Split_Middle$Age==Middle_age_ps_seq[1+p]] <- Middle_relswitchage_ps_seq[1+p]
  p=p+1
}


# Final data wrangling parts: 

# Bind final df for modelling (with respect to diet-based on actual fly DOBs)
Modelling_Data <- rbind(Split_Start,Split_Middle)
Simple_Modelling_Data <- Modelling_Data # Pull this aside so treatments are still unique.

# final prep modelling data
Modelling_Data$Treatment <- gsub("2%", "CONTINUOUS", Modelling_Data$Treatment) # fixes dataset for interactive model 
Modelling_Data$Treatment <- gsub("8%", "CONTINUOUS", Modelling_Data$Treatment) # can't be having 4 'groups' think!!
Modelling_Data$Treatment <- gsub("SU", "SWITCH", Modelling_Data$Treatment)
Modelling_Data$Treatment <- gsub("SD", "SWITCH", Modelling_Data$Treatment)
Modelling_Data$Treatment <- as.factor(Modelling_Data$Treatment)
Modelling_Data$Cage <- as.factor(Modelling_Data$Cage)
Modelling_Data$Diet <- as.factor(Modelling_Data$Diet)
str(Modelling_Data)

# final prep simple modelling data 
str(Simple_Modelling_Data)  
Simple_Modelling_Data$Cage <- as.factor(Simple_Modelling_Data$Cage)
Simple_Modelling_Data$Diet <- as.factor(Simple_Modelling_Data$Diet)





## New input: ####


#Mirre: we have two sets of cages:
#---one with transfer = 1 (cages 1-20), these switch on age 22.
  # yes, they are 20 days old at the point they are switched, but age 22 at first demography record (age >=22 as first day new diet)

#---one with tranfer = 2, these switch on age 21 >cages 21 and up
  # yes, they are age 19 at point of switch, so age 21 is first record; must be included >= 21.


#We want to model the long switch group, compared to the 8% continuous group. We can subset the data for ease first before we change the intervals.

#---- Select treatment comparator to set intervals across. 
  # We no not need to make intervals for continuous diets as these keep the same Diet throughout  
    dataset_t=Split_Data[which(Split_Data$Treatment=='SU'),] #----- SU comparator||
    dataset_t=Split_Data[which(Split_Data$Treatment=='SD'),] #----- SD comparator||
    dataset_t=Split_Data[which(Split_Data$Treatment=='8%'),] #----- 8% comparator||
    
# Summary stats 
 table(dataset_t$Age)
 unique(dataset_t$Age)
 aggregate(dataset_t$Dead~dataset_t$Age,FUN=sum)



##---------------first interval for switch is the same as the previous continous (e.g. SU = 2% < Switch) ####
  ## interval 1 - (for ease we just make dataframes that we stitch together later).
  interval1=dataset_t
  interval1$Age1=1 #start of all cages approx.
  interval1$Age2=interval1$Age
     # if flies die in this interval or are censored in this interval keep their original scores.
     # if they have not died or have not been censored in that interval, censor them at end of interval, 
     # Age2 = 20 for T=1, Age2 = 19 for T  =2
  
  #set ages
  interval1$Age2[which(interval1$Age>=20 & interval1$Cage<21)]=20
  interval1$Age2[which(interval1$Age>=19 & interval1$Cage>20)]=19
  table(interval1$Dead)
  #set dead column, to censor individuals stilll alive at end of interval - skip for continuous only
  interval1$Dead[which(interval1$Age>=20 & interval1$Cage<21)]=0
  interval1$Dead[which(interval1$Age>=19 & interval1$Cage>20)]=0
  #set treatment column
              
  interval1$interval="2%" #----- When SU comparator|| just the factor name for period before switch!  
  interval1$interval="2%<switch" #----- for 2% v SU model to see whether 2%=2%.... 
  interval1$interval="8%" #----- When SD comparator|| just the factor name for period before switch! 
  interval1$interval="8%<switch" #----- for 8% v SU model to see whether 8%=8%.... 
  
  table(interval1$Dead)


##----------interval 2 is 2 days after switch - i.e. SWITCHING PERIOD####
   ##interval 2
     interval2=dataset_t
     interval2$Age2=interval2$Age
     interval2$Age1=NA
   #set ages
     interval2$Age1[which(interval2$Cage<21)]=20
     interval2$Age1[which(interval2$Cage>20)]=19
   
     interval2$Age2[which(interval2$Age!=(22) & interval2$Cage<21)]=22 # Mirre adjusted the values here
     interval2$Age2[which(interval2$Age!=(21) & interval2$Cage>20)]=21
   #set dead column
     interval2$Dead[which(interval2$Age!=(22) & interval2$Cage<21)]=0
     interval2$Dead[which(interval2$Age!=(21) & interval2$Cage>20)]=0
   #remove the ones that have already died or been censored prior to this interval
     interval2$Dead[which(interval2$Age<(22) & interval2$Cage<21)]=NA 
     interval2$Dead[which(interval2$Age<(21) & interval2$Cage>20)]=NA
   #set treatment column
     interval2$interval="switch"
     table(interval2$Dead)


##----------interval 3 is 2 days after switch (after the observed response from the switching period / 4 days from vial swap)####
   ##interval 3
    interval3=dataset_t
     interval3$Age2=interval3$Age
     interval3$Age1=NA
   #set ages
     interval3$Age1[which(interval3$Cage<21)]=20+2
     interval3$Age1[which(interval3$Cage>20)]=19+2
     interval3$Age2[which(interval3$Age!=(22+2) & interval3$Cage<21)]=22+2
     interval3$Age2[which(interval3$Age!=(21+2) & interval3$Cage>20)]=21+2
   #set dead column
     interval3$Dead[which(interval3$Age!=(22+2) & interval3$Cage<21)]=0
     interval3$Dead[which(interval3$Age!=(21+2) & interval3$Cage>20)]=0
   #remove the ones that have already died or been censored prior to this interval
     interval3$Dead[which(interval3$Age<(22+2) & interval3$Cage<21)]=NA 
     interval3$Dead[which(interval3$Age<(21+2) & interval3$Cage>20)]=NA
   #set treatment column
     interval3$interval="2 days after switch"
     table(interval3$Dead)


##--------- interval 4 is 4 days after switch ####
   ##interval4
    interval4=dataset_t
     interval4$Age2=interval4$Age
     interval4$Age1=NA
   #set ages
     interval4$Age1[which(interval4$Cage<21)]=20+2+2
     interval4$Age1[which(interval4$Cage>20)]=19+2+2
     interval4$Age2[which(interval4$Age!=(22+2+2) & interval4$Cage<21)]=22+2+2
     interval4$Age2[which(interval4$Age!=(21+2+2) & interval4$Cage>20)]=21+2+2
   #set dead column
     interval4$Dead[which(interval4$Age!=(22+2+2) & interval4$Cage<21)]=0
     interval4$Dead[which(interval4$Age!=(21+2+2) & interval4$Cage>20)]=0
   #remove the ones that have already died or been censored prior to this interval
     interval4$Dead[which(interval4$Age<(22+2+2) & interval4$Cage<21)]=NA 
     interval4$Dead[which(interval4$Age<(21+2+2) & interval4$Cage>20)]=NA
   #set treatment column
     interval4$interval="4 days after switch"
     table(interval4$Dead)


##----------interval 5 is 6 days after switch ####
   ##interval5
     interval5=dataset_t
     interval5$Age2=interval5$Age
     interval5$Age1=NA
   #set ages
     interval5$Age1[which(interval5$Cage<21)]=20+2+2+2
     interval5$Age1[which(interval5$Cage>20)]=19+2+2+2
     interval5$Age2[which(interval5$Age!=(22+2+2+2) & interval5$Cage<21)]=22+2+2+2
     interval5$Age2[which(interval5$Age!=(21+2+2+2) & interval5$Cage>20)]=21+2+2+2
   #set dead column
     interval5$Dead[which(interval5$Age!=(22+2+2+2) & interval5$Cage<21)]=0
     interval5$Dead[which(interval5$Age!=(21+2+2+2) & interval5$Cage>20)]=0
   #remove the ones that have already died or been censored prior to this interval
     interval5$Dead[which(interval5$Age<(22+2+2+2) & interval5$Cage<21)]=NA 
     interval5$Dead[which(interval5$Age<(21+2+2+2) & interval5$Cage>20)]=NA
   #set treatment column
     interval5$interval="6 days after switch"
     table(interval5$Dead)


##----------interval 6 is >=8 days after switch (so >= 10 days from vial swap technically) ####
   #interval6 rest 
     interval6=dataset_t
     interval6$Age2=interval6$Age
     interval6$Age1=NA
   #set ages
     interval6$Age1[which(interval5$Cage<21)]=20+2+2+2+2
     interval6$Age1[which(interval5$Cage>20)]=19+2+2+2+2
   #no need to set Age2 as all will die or be censored now. only remove ones that died or have been censored before this interval.
   #set dead column
   #remove the ones that have already died or been censored prior to this interval
     interval6$Dead[which(interval6$Age<(22+2+2+2+2) & interval6$Cage<21)]=NA 
     interval6$Dead[which(interval6$Age<(21+2+2+2+2) & interval6$Cage>20)]=NA
   #set treatment column
     interval6$interval=">8 days after switch"
     table(interval6$Dead)


##---------- combine datasets ####
   ### < switch | switch | 2 days | 4 days | 6 days | >=8 days ###
     
   #---- SET BASE FOR COMPARISON --
     dataset_BASE=Split_Data[which(Split_Data$Treatment=='8%'),] #------ 8%
     dataset_BASE=Split_Data[which(Split_Data$Treatment=='2%'),] #------ 2%
     
     dataset_BASE$Age1=1
     dataset_BASE$Age2=dataset_BASE$Age
        
       dataset_BASE$interval="8%" #------ 8%
       dataset_BASE$interval="2%" #------ 2%
       
     str(dataset_BASE)
     str(interval1)
     
  # combine final interval frame
     dataset_timedep=rbind(dataset_BASE,interval1,interval2,interval3,interval4,interval5,interval6)
     dataset_timedep=rbind(dataset_BASE,interval1)
     str(dataset_timedep)
   #remove rows with NA in Dead
     dataset_timedep=dataset_timedep[which(is.na(dataset_timedep$Dead)==F),]
     
   #final df adjustment
     dataset_timedep$interval=factor(dataset_timedep$interval)
     
     dataset_timedep$interval=relevel(dataset_timedep$interval,ref="8%")  #------ 8%
     dataset_timedep$interval=relevel(dataset_timedep$interval,ref="2%")  #------ 2%
     
     
     str(dataset_timedep) # chek base levelsl etc.

  # save separate dataframes (prior to re-running above to obtain the other)
    switch_rich_timedep <- dataset_timedep       # 8% v SU (example)
    switch_restricted_timedep <- dataset_timedep # 2% v SD (mine)
    
    switch_restricted_delay <- dataset_timedep       # 8% v SD (does it take longer to starve?)
    switch_rich_rapid <- dataset_timedep # 2% v SU (is it quicker to respond to food when starved?)
    
    cont_model_data <- dataset_timedep
    
    # save all modelling dataframes now that they've been created (can load back in if needed)
    write.csv(switch_rich_timedep, "switch_rich_timedep.csv", row.names = FALSE)
    write.csv(switch_restricted_timedep, "switch_restricted_timedep.csv", row.names = FALSE)
    write.csv(switch_restricted_delay, "switch_restricted_delay.csv", row.names = FALSE)
    write.csv(switch_rich_rapid, "switch_rich_rapid.csv", row.names = FALSE)
    
    write.csv(cont_model_data, "cont_model_data.csv", row.names = FALSE)

##---------- Fit the non-interval based model (simple, no need to censor anything) ####

    # Model 1 for NON-INTERVAL model (4 treatment, across all time. Censored flies are  censors from experiment only).
    
    cont_model_data <- Split_Data
    table(cont_model_data$Dead)
    
    # censor all switch flies - absolutely no different than using the interval method above!  :) 
    cont_model_data$Dead[cont_model_data$Treatment=="SU"]=NA
    cont_model_data$Dead[cont_model_data$Treatment=="SD"]=NA
    cont_model_data$Treatment[cont_model_data$Treatment=="SU"]=NA
    cont_model_data$Treatment[cont_model_data$Treatment=="SD"]=NA
    cont_model_data=droplevels(cont_model_data)
    table(cont_model_data$Dead)
    str(cont_model_data)
    
    write.csv(cont_model_data, "cont_model_data.csv", row.names = FALSE)
    
    # Simple additive model for treatments (base )
    levels(cont_model_data$Treatment)
    cont_model_data$Treatment<-relevel(cont_model_data$Treatment,ref="8%") # make 8% the base?? is this more appropriate?
    fit1<-coxme(Surv(Age, Dead)~Treatment+(1|Cage), data=cont_model_data) # so no Treatment*Diet, just 'interval'
    summary(fit1)
    anova(fit1)       
    
    install.packages("StatWriter")
    
    # table it up
    fit1_table <- summary(fit1)
    fit1_table <- tidy(fit1_table)
    fit1_table$exp.coef <- exp(fit1$coefficients)
    fit1_table$Xsq <- unlist(anova(fit1)[2])[2]
    fit1_table$df <- fit1$df[1]
    fit1_table$n.included <- fit1$n[1]
    fit1_table$n.censored <- fit1$n[2]-fit1$n[1]
    fit1_table$n.tot.mod <- fit1$n[2]
    fit1_table$n.excluded <- length(cont_model_data$Age)-fit1$n[2]
    
    
###------------ Fit the models on the interval-based datasets. 
# https://www.researchgate.net/post/What_to_report_from_a_Cox_Proportional_Hazards_Regression_analysis - reporting output advice
   library(coxme)
  
  # Model 2 for INTERVAL-based 8% v SU 
     
   fit2<-coxme(Surv(Age1,Age2, Dead)~interval+(1|Cage), data=switch_rich_timedep) # so no Treatment*Diet, just 'interval'
   summary(fit2)
   anova(fit2)
   
   # table it up
   fit2_table <- summary(fit2)
   fit2_table <- tidy(fit2_table)
   fit2_table$exp.coef <- exp(fit2$coefficients)
   fit2_table$Xsq <- unlist(anova(fit2)[2])[2]
   fit2_table$df <- fit2$df[1]
   fit2_table$n.included <- fit2$n[1]
   fit2_table$n.censored <- fit2$n[2]-fit2$n[1]
   fit2_table$n.tot.mod <- fit2$n[2]
   fit2_table$n.excluded <- length(switch_rich_timedep$Age)-fit2$n[2]
   
   # 8% v SU normalisation plotting
   plot_fit2 <- fit2_table
   plot_fit2$index <- c("6","3","1","4","5","2")
   plot_fit2$term <- sub(",.*", "", substring(plot_fit2$term, 9)) #sack off first 9 chars
   plot_fit2<- plot_fit2[order(plot_fit2$index), ]
   plot_fit2$term <- factor(plot_fit2$term, levels = plot_fit2$term[order(plot_fit2$index)])
   plot_fit2$term
   str(plot_fit2)
   plot_fit2$ymin <-  plot_fit2$estimate - plot_fit2$std.error
   plot_fit2$ymax <-  plot_fit2$estimate + plot_fit2$std.error
   
   plot(plot_fit2$term, plot_fit2$estimate) # quick check... 
   
   gg_fit2 <- 
     ggplot(plot_fit2) +
     aes(x = term , y = estimate, group=1) +
     geom_line(colour = "cornflowerblue", size = 1) +
     geom_point(shape=21, fill="black", color="cornflowerblue", size=2, stroke=1)+
     geom_errorbar(ymin=plot_fit2$conf.low, ymax=plot_fit2$conf.high, width = 0.2, size=1, colour = "gray20") +
     labs(x=element_blank(), y="ln(Hazard ratio)") +
     scale_y_continuous(limits=c(-3,3), breaks=seq(-3,3,1)) +
     theme_classic() +
     theme(axis.title.x = element_text(face="bold", size=17),
           axis.text.x = element_text(face= 'bold', size=13, colour = 'gray30', angle = 45, hjust = 1)) +
     theme(axis.title.y = element_text(face="bold", size=17),
           axis.text.y = element_text(face= 'bold', size=13, colour = 'gray30')) +
     theme(axis.line = element_line(colour = 'black', size = 1)) +
     theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
     theme(plot.margin = unit(c(0.5, 0, 0, 0), "cm")) +
     theme(aspect.ratio = 1) +
     geom_hline(yintercept = 0, linetype = "dashed", colour = "red", size=1) 


  # Model 3 for INTERVAL-based 2% v SD 
   fit3<-coxme(Surv(Age1,Age2, Dead)~interval+(1|Cage), data=switch_restricted_timedep) # so no Treatment*Diet, just 'interval'
   summary(fit3)
   anova(fit3)   
   
   fit3_table <- summary(fit3)
   fit3_table <- tidy(fit3_table)
   fit3_table$exp.coef <- exp(fit3$coefficients)
   fit3_table$Xsq <- unlist(anova(fit3)[2])[2]
   fit3_table$df <- fit3$df[1]
   fit3_table$n.included <- fit3$n[1]
   fit3_table$n.censored <- fit3$n[2]-fit3$n[1]
   fit3_table$n.tot.mod <- fit3$n[2]
   fit3_table$n.excluded <- length(switch_restricted_timedep$Age)-fit3$n[2]
  
   
   # 2% v SD normalisation plotting
   plot_fit3 <- fit3_table
   plot_fit3$index <- c("1","6","3","4","5","2")
   plot_fit3$term <- sub(",.*", "", substring(plot_fit3$term, 9)) #sack off first 9 chars
   plot_fit3<- plot_fit3[order(plot_fit3$index), ]
   plot_fit3$term <- factor(plot_fit3$term, levels = plot_fit3$term[order(plot_fit3$index)])
   plot_fit3$term
   str(plot_fit3)
   plot_fit3$ymin <-  plot_fit3$estimate - plot_fit3$std.error
   plot_fit3$ymax <-  plot_fit3$estimate + plot_fit3$std.error
   
   plot(plot_fit3$term, plot_fit3$estimate) # quick check... 
   
   gg_fit3 <- 
     ggplot(plot_fit3) +
     aes(x = term , y = estimate, group=1) +
     geom_line(colour = "goldenrod", size = 1) +
     geom_point(shape=21, fill="black", color="goldenrod", size=2, stroke=1)+
     geom_errorbar(ymin=plot_fit3$conf.low, ymax=plot_fit3$conf.high,width = 0.2, size=1, colour = "gray 20") +
     geom_hline(yintercept = 0, linetype = "dashed", colour = "black", size=1) +
     labs(x=element_blank() ,y="ln(Hazard ratio)" ) +
     scale_y_continuous(limits=c(-3,3), breaks=seq(-3,3,1)) +
     theme_classic() +
     theme(axis.title.x = element_text(face="bold", size=17),
           axis.text.x = element_text(face= 'bold', size=13, colour = 'gray30', angle = 45, hjust = 1)) +
     theme(axis.title.y = element_text(face="bold", size=17),
           axis.text.y = element_text(face= 'bold', size=13, colour = 'gray30')) +
     theme(axis.line = element_line(colour = 'black', size = 1)) +
     theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
     theme(plot.margin = unit(c(0.5, 0, 0, 0), "cm")) +
     theme(aspect.ratio = 1) 
     
   
   library(cowplot)
### Construct 'normalisation' figure (3) ---------------
   plot_grid(gg_fit2,gg_fit3, labels = c("A) Switch to rich","B) Switch to DR"), hjust = 0) +
   draw_label(
     label="Interval", hjust=0.6, x = 0.515, y = 0.05, color = "black", size = 18,fontface = "bold")  +
     draw_label(label=expression(paste("Continous \n     rich")), # left line label
       hjust=0.5, x = 0.11, y = 0.6025, color = "red", size = 10, fontface = "italic") +
     draw_label(label=expression(paste("Continous \n     DR")), # right line label
       hjust=0.5, x = 0.61, y = 0.6025, color = "black", size = 10, fontface = "italic") 
   
   
   

# Model 4&5.... simplified?
   ## Use the non-interval data, censor relswitch!=0, then compare.
   ## Set a continouous as base, the corresponding switch should not be sig, whilst both others should.
   
   # 2% base against SU
   cont_switch_comp_2 <- subset(Simple_Modelling_Data, relSwitchAge=="0") # pull relevant data.
   levels(cont_switch_comp_2$Treatment)
   cont_switch_comp_2$Dead[cont_switch_comp_2$Treatment=="8%"]=NA
   cont_switch_comp_2$Dead[cont_switch_comp_2$Treatment=="SD"]=NA
   cont_switch_comp_2$Treatment[cont_switch_comp_2$Treatment=="8%"]=NA
   cont_switch_comp_2$Treatment[cont_switch_comp_2$Treatment=="SD"]=NA
   cont_switch_comp_2=droplevels(cont_switch_comp_2)
   levels(cont_switch_comp_2$Treatment)
   
   fit_csc2 <- coxme(Surv(Age,Dead)~Treatment+(1|Cage), data=cont_switch_comp_2)
   summary(fit_csc2)
   anova(fit_csc2)
   
   # table it up
   fit4_table <- summary(fit_csc2)
   fit4_table <- tidy(fit4_table)
   fit4_table$exp.coef <- exp(fit_csc2$coefficients)
   fit4_table$Xsq <- unlist(anova(fit_csc2)[2])[2]
   fit4_table$df <- fit_csc2$df[1]
   fit4_table$n.included <- fit_csc2$n[1]
   fit4_table$n.censored <- fit_csc2$n[2]-fit_csc2$n[1]
   fit4_table$n.tot.mod <- fit_csc2$n[2]
   fit4_table$n.excluded <- length(cont_switch_comp_2$Age)-fit_csc2$n[2]
 
   # 8% base against SD 
   cont_switch_comp_8 <- subset(Simple_Modelling_Data, relSwitchAge=="0") 
   levels(cont_switch_comp_8$Treatment)
   cont_switch_comp_8$Dead[cont_switch_comp_8$Treatment=="2%"]=NA
   cont_switch_comp_8$Dead[cont_switch_comp_8$Treatment=="SU"]=NA
   cont_switch_comp_8$Treatment[cont_switch_comp_8$Treatment=="2%"]=NA
   cont_switch_comp_8$Treatment[cont_switch_comp_8$Treatment=="SU"]=NA
   cont_switch_comp_8=droplevels(cont_switch_comp_8)
   levels(cont_switch_comp_8$Treatment)
   
   fit_csc8 <- coxme(Surv(Age,Dead)~Treatment+(1|Cage), data=cont_switch_comp_8)
   summary(fit_csc8)
   anova(fit_csc8)  
   
   # table it up
   fit5_table <- summary(fit_csc8)
   fit5_table <- tidy(fit5_table)
   fit5_table$exp.coef <- exp(fit_csc8$coefficients)
   fit5_table$Xsq <- unlist(anova(fit_csc8)[2])[2]
   fit5_table$df <- fit_csc8$df[1]
   fit5_table$n.included <- fit_csc8$n[1]
   fit5_table$n.censored <- fit_csc8$n[2]-fit_csc8$n[1]
   fit5_table$n.tot.mod <- fit_csc8$n[2]
   fit5_table$n.excluded <- length(cont_switch_comp_8$Age)-fit_csc8$n[2]
   
##----- as expected, the switches ARE comparable to the previous treatments, before switching.  
   
   
   # Model 6 -  fit model for non-interval interactive dataset
   levels(Modelling_Data$Treatment);levels(Modelling_Data$Diet) # levels
   Modelling_Data$Diet<-relevel(Modelling_Data$Diet,ref="8%") # Set 8% (rich) as diet base for comparison
   fit6<-coxme(Surv(Age,Dead)~Treatment*Diet+(1|Cage), data=Modelling_Data)
   summary(fit6)
   anova(fit6)  
   
   # table it up
   fit6_table <- summary(fit6)
   fit6_table <- tidy(fit6_table)
   fit6_table$exp.coef <- exp(fit6$coefficients)
   fit6_table$Xsq <- unlist(anova(fit6)[2])[2]
   fit6_table$df <- fit6$df[1]
   fit6_table$n.included <- fit6$n[1]
   fit6_table$n.censored <- fit6$n[2]-fit6$n[1]
   fit6_table$n.tot.mod <- fit6$n[2]
   fit6_table$n.excluded <- length(Modelling_Data$Age)-fit6$n[2]

   

   
# Need[?] to report the P val, chisq value, wald, sample size? hazard ratios, confidence intervals of HR's.
   
   
## Create summary statistics from the kaplan data ####

   # calculate median, maximal lifespan,  and confidence intervals for drosophila age on different diets.
   # think this is very useful.... why hasn't this been used before? 
   
   ## Needs to use the kaplan output, not simply the data long (which even includes the censors). ## 
   
   install.packages("DescTools")
   library(DescTools)
   ?MedianCI
 
   ## run kaplan section - made use of the kaplan stats.... no idea whether it's useful, but was easy enough. 
   kaplan_data <- plot.data
    
   basic_stats_by_treatment <- kaplan_data %>% 
     group_by(Treatment) %>% 
     summarise(
       lowerCI = MedianCI(time, conf.level = 0.95, sides = "two.sided")[2],
       median = MedianCI(time, conf.level = 0.95, sides = "two.sided")[1],
       upperCI = MedianCI(time, conf.level = 0.95, sides = "two.sided")[3],
       max = max(time),
       std_x = sd(time),
       surv= mean(surv),
       std.err_x = mean(std.err[!is.infinite(std.err)]), # R can deal with NA, but must specify Inf omissions!
       upper = mean(upper[!is.na(upper)]),
       lower = mean(lower[!is.na(lower)]),
       n = sum(n.event) # crude but correct; total number of flies per condition.
     );print(basic_stats_by_treatment)
      # switch rich is fractionated as there are an even number of events... fine. 
      # I believe I have ALL the data I could need to do my writing now. Sorted!! Write up tonight!! :D 
   
   
   
   
   ###### Methods plot ##### no too technical - shelved. Made a more 'graphical' one - see what mirre says... 
   
   methods <- matrix(,nrow = 90, ncol = 3); methods <- as.data.frame(methods)
   names(methods) <- c("day","Continuous rich","Continuous restricted")
   methods$day <- seq(1:90)
   
   mysigmoid <- function(n, min, max, maxe = 4, mine = -maxe) {
     y <- 1 / (1 + exp(-seq(mine, maxe, len = max(0, n))))
     min + (max - min) * (y - y[1]) / (y[n] - y[1])
   }
   
   methods$`Continuous rich`[1:50] <-       round(mysigmoid(50, -7, -1, 2, -1), 2) 
   methods$`Continuous restricted`[1:90] <- round(mysigmoid(90, -7, -1, 2, -1), 2) 
  
   # now need to melt into longdata 
   plot.methods <- melt(methods, id = "day")
   names(plot.methods) <- c("day", "treatment", "logmx")
   
   # gg_fit3 <- 
   
     ggplot(plot.methods, aes(day, logmx, colour = treatment)) +
       geom_line()+
     
     geom_vline(xintercept = 19, linetype = "dotted", colour = "dark grey", size=1) +
     scale_y_continuous(limits=c(-7.5,-0.5), breaks=seq(-7,-1,1)) +
     scale_x_continuous(limits=c(0,95), breaks=seq(0,90,10)) +
     labs(x="Days" ,y="log(mx)") +
     theme_classic() +
       theme(axis.title.x = element_text(face="bold", size=17),
             axis.text.x = element_text(face= 'bold', size=13, colour = 'gray30')) +
       theme(axis.title.y = element_text(face="bold", size=17),
             axis.text.y = element_text(face= 'bold', size=13, colour = 'gray30')) +
       theme(axis.line = element_line(colour = 'black', size = 1)) +
       theme(axis.ticks = element_line(colour = "black", size = 1.5)) +
       theme(plot.margin = unit(c(0.5, 0, 0, 0), "cm")) +
       theme(aspect.ratio = 1) 
     
     
     
     
     
     geom_line(colour = "goldenrod", size = 1) +
     geom_point(shape=21, fill="black", color="goldenrod", size=2, stroke=1)+
     geom_errorbar(ymin=plot_fit3$conf.low, ymax=plot_fit3$conf.high,width = 0.2, size=1, colour = "gray 20") +
     geom_hline(yintercept = 0, linetype = "dashed", colour = "black", size=1) +
     


   
  
   