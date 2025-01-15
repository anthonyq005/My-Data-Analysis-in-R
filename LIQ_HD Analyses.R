#LIQ-HD:
# Read and store LIQ-HD Analysis data:
setwd("/Users/annabelle/Downloads/2024/LIQ_HD")
library(readr)
LIQHDpilot <- as.data.frame(read_csv("LIQ HD 7.19.24.csv")) #This contains the LIQ-HD Data from 1 female mouse
View(LIQHDpilot)

LIQHDpilot_mouse1 <- as.data.frame(LIQHDpilot$datetime)
names(LIQHDpilot_mouse1) [1] <- "DateTime"
substr(LIQHDpilot_mouse1$DateTime, 11,16)
LIQHDpilot_mouse1$Time <- substr(LIQHDpilot_mouse1$DateTime, 11,16)
LIQHDpilot_mouse1$ExperimentalSide <- LIQHDpilot$Experimental_side
LIQHDpilot_mouse1$LickNumberLeft <- LIQHDpilot$LickNumber13
LIQHDpilot_mouse1$LickNumberRight <- LIQHDpilot$LickNumber14
LIQHDpilot_mouse1$LickDurationLeft <- LIQHDpilot$LickDuration13
LIQHDpilot_mouse1$LickDurationRight <- LIQHDpilot$LickDuration14
LIQHDpilot_mouse1$BoutNumberLeft <- LIQHDpilot$BoutNumber13
LIQHDpilot_mouse1$BoutNumberRight <- LIQHDpilot$BoutNumber14
LIQHDpilot_mouse1$BoutDurationLeft <- LIQHDpilot$BoutDuration13
LIQHDpilot_mouse1$BoutDurationRight <- LIQHDpilot$BoutDuration14
LIQHDpilot_mouse1$BoutLickNumberLeft <- LIQHDpilot$BoutLickNumber13
LIQHDpilot_mouse1$BoutLickNumberRight <- LIQHDpilot$BoutLickNumber14
LIQHDpilot_mouse1$BoutLickDurationLeft <- LIQHDpilot$BoutLickDuration13
LIQHDpilot_mouse1$BoutLickDurationRight <- LIQHDpilot$BoutLickDuration14

View(LIQHDpilot_mouse1)

# Store new data table in a csv file:
write.csv(LIQHDpilot_mouse1, "~/Downloads/2024/LIQ_HD/LIQHD_Mouse1.csv", row.names=FALSE)


#The "names" show the title of each column
names(CREST60X_data)
head(CREST60X_data)
CREST60X_data_2 <- as.data.frame(CREST60X_data$Filename)
names(CREST60X_data_2) [1] <- "Filename"
View(CREST60X_data_2)
substr(CREST60X_data$Filename, 1,1)
CREST60X_data_2$MouseID <- substr(CREST60X_data$Filename, 1,1)
substr(CREST60X_data$Filename, 3,3)
CREST60X_data_2$Sex <- substr(CREST60X_data$Filename, 3,3)
CREST60X_data_2$ROI <- substr(CREST60X_data$Filename, 37, 40)
nchar(CREST60X_data$Filename)
CREST60X_data_2$Slice<-substr(CREST60X_data$Filename, nchar(CREST60X_data$Filename)-14, nchar(CREST60X_data$Filename)-9)

#CREST60X_data$MouseName<-substr(CREST60X_data$Filename, 3,7)
CREST60X_data_2$CodedMouseID <- CREST60X_data_2$MouseID

myFun <- function(n = 5000) {
  a <- do.call(paste0, replicate(1, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%02d", sample(999, n, TRUE)))}

#paste(CREST60X_data$CodedMouseID, myFun(98), sep = "-")
CREST60X_data_2$CodedMouseID [1:6] <- myFun(1)
CREST60X_data_2$CodedMouseID [7:12] <- myFun(1)
CREST60X_data_2$CodedMouseID [13:18] <- myFun(1)
CREST60X_data_2$CodedMouseID [19:24] <- myFun(1)
CREST60X_data_2$CodedMouseID [25:30] <- myFun(1)
CREST60X_data_2$CodedMouseID [31:36] <- myFun(1)
CREST60X_data_2$CodedMouseID [37:42] <- myFun(1)
CREST60X_data_2$CodedMouseID [43:48] <- myFun(1)
CREST60X_data_2$CodedMouseID [49:54] <- myFun(1)
CREST60X_data_2$CodedMouseID [55:60] <- myFun(1)
CREST60X_data_2$CodedMouseID [61:66] <- myFun(1)
CREST60X_data_2$CodedMouseID [67:72] <- myFun(1)
CREST60X_data_2$CodedMouseID [73:78] <- myFun(1)
CREST60X_data_2$CodedMouseID [79:84] <- myFun(1)
CREST60X_data_2$CodedMouseID [85:90] <- myFun(1)
CREST60X_data_2$CodedMouseID [91:96] <- myFun(1)
CREST60X_data_2$CodedMouseID [97:102] <- myFun(1)
CREST60X_data_2$CodedMouseID [103:108] <- myFun(1)
CREST60X_data_2$CodedMouseID [109:114] <- myFun(1)


#CREST60X_data$Name <- substr(CREST60X_data$Filename, 4, 7)
#CREST60X_data$Name
#Assign Treatment: First, give all data the treatment of "01.Control", then assign those in the Alcohol group to "02.Alcohol"). 
CREST60X_data_2$Treatment [1:114] <- "01.Control" 

CREST60X_data_2$Treatment [1:4] <- "02.Alcohol" 
CREST60X_data_2$Treatment [13:18] <- "02.Alcohol"
CREST60X_data_2$Treatment [61:66] <- "02.Alcohol"
CREST60X_data_2$Treatment [73:102] <- "02.Alcohol"
CREST60X_data_2$Treatment [109:114] <- "02.Alcohol"

#CREST60X_data$MouseID [1:114] <- "Clear"

#Replace the Sex labels "M" and "F" to "01.Male" and "02.Female"
CREST60X_data_2$Sex <- as.character(CREST60X_data_2$Sex)
CREST60X_data_2$Sex[CREST60X_data_2$Sex == "M"] <- "01.Male"
CREST60X_data_2$Sex[CREST60X_data_2$Sex == "F"] <- "02.Female"
CREST60X_data_2$Sex <- as.factor(CREST60X_data_2$Sex)
names(CREST60X_data_2)
names(CREST60X_data)

#Add data to Crest60X_data_2
CREST60X_data_2$AllCRFcObjectCount <- CREST60X_data$AllCRFcObjectCount
CREST60X_data_2$AllGFPcObjectCount <- CREST60X_data$AllGFPcObjectCount
CREST60X_data_2$AllDAPIObjectCount <- CREST60X_data$AllDAPIObjectCount
CREST60X_data_2$CRFcellcount3D <- CREST60X_data$CRFcellcount3D
CREST60X_data_2$GFPcellcount3D <- CREST60X_data$GFPcellcount3D
CREST60X_data_2$Hoechstcellcount3D <- CREST60X_data$Hoechstcellcount3D
CREST60X_data_2$AllCRFcObjectCount <- CREST60X_data$AllCRFcObjectCount
CREST60X_data_2$AllCRFcObjectCount <- CREST60X_data$AllCRFcObjectCount

#Calculate the percent difference between the two measurements (2D vs 3D):

Diff <- function(measurement2, measurement1) {(measurement1 - measurement2)}

CREST60X_data_2$Diff_CRF <- Diff(CREST60X_data_2$AllCRFcObjectCount, CREST60X_data_2$CRFcellcount3D)
CREST60X_data_2$Diff_GFP <- Diff(CREST60X_data_2$AllGFPcObjectCount, CREST60X_data_2$GFPcellcount3D)
CREST60X_data_2$Diff_Hoechst <-Diff(CREST60X_data_2$AllDAPIObjectCount, CREST60X_data_2$Hoechstcellcount3D)

# Store new data table in a csv file:
write.csv(CREST60X_data_2, "~/Downloads/2024/Modified Drinking in the Dark/CREST60X_data_2_CellsonlyPractice_ROI4.csv", row.names=FALSE)


#Create a function that calculates the ratio (also means fraction):
Ratio <- function(FractionValue, TotalValue) {FractionValue/TotalValue}
Percent <- function(Ratio) {Ratio * 100}
TotalSum <- function(Total1, Total2) {Total1 + Total2}

#Define the constant value of the section thickness (35 micrometers):
image_60X_px <- 1608
image_60X_mmperpx <- 188.4 /1000
image_60X_mm <- image_60X_px / image_60X_mmperpx
image_60X_mm3 <- image_60X_mm * image_60X_mm 


sectionthickness_um <- 3203340

# XYZ for density measurements: Size of image is 302.52 um x 302.52 um x 18 um

sectionthickness_mm3 <- (302.52 * 302.52 * 18)/(1000*1000*1000)
sectionthickness_mm2 <-(302.52 * 302.52)/(1000*1000)

302.52 *302.52 * 35
3203142/1000000000
495 / 0.003203142

CREST60X_data$celldensity3D <- CREST60X_data$Hoechstcellcount3D / sectionthickness_mm3
CREST60X_data$celldensity2D <- CREST60X_data$AllDAPIObjectCount / sectionthickness_mm2

CREST60X_data$celldensity2D <- CREST60X_data$AllDAPIObjectCount / sectionthickness_mm3

#Apply functions:

# mGFPf = membrane-bound GFP fibers
# CRFp = CRF peptide not co-locolized with DAPI 
# CRFcp = CRF cell-peptide (CRF peptide co-locolized with DAPI)
# GFPxCRF = GFP touching CRF
# CRFxGFP = CRF touching GFP
# CRFANDGFP = overlap between CRF and GFP

#A. Volume measurement to % of Area covered by X

  #This is Hoechst expression measured as percent of area covered:
CREST60X_data_2$RatioAreaTotalHoechst <- Ratio(CREST60X_data$TotalHoechstVolume, sectionthickness_um)
CREST60X_data_2$PercentAreaTotalHoechst <- Percent(CREST60X_data_2$RatioAreaTotalHoechst) 

  #This is CRF and GFP expression measured as percent of area covered:
CREST60X_data_2$RatioAreaTotalCRF <- Ratio(CREST60X_data$TotalCRFVolume, sectionthickness_um)
CREST60X_data_2$PercentAreaTotalCRF <- Percent(CREST60X_data_2$RatioAreaTotalCRF) 

CREST60X_data_2$RatioAreaTotalGFP <- Ratio(CREST60X_data$TotalGFPVolume, sectionthickness_um)
CREST60X_data_2$PercentAreaTotalGFP <- Percent(CREST60X_data_2$RatioAreaTotalGFP) 

CREST60X_data_2$RatioAreaCRFp <- Ratio(CREST60X_data$CRFpVolume, sectionthickness_um)
CREST60X_data_2$PercentAreaCRFp <- Percent(CREST60X_data_2$RatioAreaCRFp) 

CREST60X_data_2$RatioAreamGFPf<- Ratio(CREST60X_data$mGFPfVolume, sectionthickness_um)
CREST60X_data_2$PercentAreamGFPf <- Percent(CREST60X_data_2$RatioAreamGFPf) 

  #This is CRFxGFP or GFPxCRF or CRFHavingGFP or GFPHavingCRF or CRFANDGFP expression measured as percent of area covered:
CREST60X_data_2$RatioAreaAllCRFxAllGFP <- Ratio(CREST60X_data$AllCRFxAllGFPVolume, CREST60X_data$ImageVolume)
CREST60X_data_2$PercentAreaAllCRFxAllGFP <- Percent(CREST60X_data_2$RatioAreaAllCRFxAllGFP) 

CREST60X_data_2$RatioAreaAllGFPxAllCRF <- Ratio(CREST60X_data$AllGFPxAllCRFVolume, CREST60X_data$ImageVolume)
CREST60X_data_2$PercentAreaAllGFPxAllCRFP <- Percent(CREST60X_data_2$RatioAreaAllGFPxAllCRF) 

CREST60X_data_2$RatioAreaCRFpHavingmGFPf <- Ratio(CREST60X_data$CRFpHavingmGFPfVolume, CREST60X_data$ImageVolume)
CREST60X_data_2$PercentAreaCRFpHavingmGFPf <- Percent(CREST60X_data_2$RatioAreaCRFpHavingmGFPf) 

CREST60X_data_2$RatioAreamGFPfHavingCRFp<- Ratio(CREST60X_data$mGFPfHavingCRFpVolume, CREST60X_data$ImageVolume)
CREST60X_data_2$PercentAreamGFPfHavingCRFp <- Percent(CREST60X_data_2$RatioAreamGFPfHavingCRFp) 

CREST60X_data_2$RatioAreaAllCRFANDAllGFP<- Ratio(CREST60X_data$AllCRFANDAllGFPVolume, CREST60X_data$ImageVolume)
CREST60X_data_2$PercentAreaAllCRFANDAllGFP <- Percent(CREST60X_data_2$RatioAreaAllCRFANDAllGFP) 

CREST60X_data_2$RatioAreaCRFpANDmGFPf<- Ratio(CREST60X_data$CRFpANDmGFPfVolume, CREST60X_data$ImageVolume)
CREST60X_data_2$PercentAreaCRFpANDmGFPf <- Percent(CREST60X_data_2$RatioAreaCRFpANDmGFPf) 

#B. Cell Density: number of cells per um3 or mm3

  #First, calculate the Image Volume from um3 to mm3: conversion is 1 cubic micrometer = 1x10-9 (ten to the minus 9) cubic milimeters
CREST60X_data_2$ImageVolumemm3 <- CREST60X_data$ImageVolume / 1000000000
  
  #This is number of cells per mm3:
CREST60X_data_2$Density_mm3_TotalHoechstCells <- Ratio(CREST60X_data$HoechstCellCount, CREST60X_data_2$ImageVolumemm3)
CREST60X_data_2$Density_mm3_CRFCells <- Ratio(CREST60X_data$CRFcellsCount, CREST60X_data_2$ImageVolumemm3)
CREST60X_data_2$Density_mm3_GFPCells <- Ratio(CREST60X_data$GFPcellsCount, CREST60X_data_2$ImageVolumemm3)

#C. Percentage of cells:

CREST60X_data_2$RatioCRFCells<- Ratio(CREST60X_data_2$Density_mm3_CRFCells, CREST60X_data_2$Density_mm3_TotalHoechstCells)
CREST60X_data_2$PercentCRFCells<-Percent(CREST60X_data_2$RatioCRFCells)

CREST60X_data_2$RatioGFPCells<-Ratio(CREST60X_data_2$Density_mm3_GFPCells, CREST60X_data_2$Density_mm3_TotalHoechstCells)
CREST60X_data_2$PercentGFPCells<-Percent(CREST60X_data_2$RatioGFPCells)

CREST60X_data_2$PercentofLabeledCells <- TotalSum(CREST60X_data_2$PercentCRFCells, CREST60X_data_2$PercentGFPCells)
CREST60X_data_2$PercentofUnlabeledCells <- 100 - CREST60X_data_2$PercentofLabeledCells 

#D. Count measured as Density of signal per mm3:

CREST60X_data_2$Density_mm3_TotalCRF <- Ratio(CREST60X_data$TotalCRFCount, CREST60X_data_2$ImageVolumemm3)
CREST60X_data_2$Density_mm3_TotalGFP <- Ratio(CREST60X_data$TotalGFPCount, CREST60X_data_2$ImageVolumemm3)
CREST60X_data_2$Density_mm3_CRFp <- Ratio(CREST60X_data$CRFpCount, CREST60X_data_2$ImageVolumemm3)
CREST60X_data_2$Density_mm3_mGFPf <- Ratio(CREST60X_data$mGFPfCount, CREST60X_data_2$ImageVolumemm3)
CREST60X_data_2$Density_mm3_AllCRFxAllGFP <- Ratio(CREST60X_data$AllCRFxAllGFPCount, CREST60X_data_2$ImageVolumemm3)
CREST60X_data_2$Density_mm3_AllGFPxAllCRF <- Ratio(CREST60X_data$AllGFPxAllCRFCount, CREST60X_data_2$ImageVolumemm3)
CREST60X_data_2$Density_mm3_CRFpHavingmGFPf <- Ratio(CREST60X_data$CRFpHavingmGFPfCount, CREST60X_data_2$ImageVolumemm3)
CREST60X_data_2$Density_mm3_mGFPfHavingCRFp <- Ratio(CREST60X_data$mGFPfHavingCRFpCount, CREST60X_data_2$ImageVolumemm3)
CREST60X_data_2$Density_mm3_AllCRFANDAllGFP <- Ratio(CREST60X_data$AllCRFANDAllGFPCount, CREST60X_data_2$ImageVolumemm3)
CREST60X_data_2$Density_mm3_CRFpANDmGFPf <- Ratio(CREST60X_data$CRFpANDmGFPfCount, CREST60X_data_2$ImageVolumemm3)

# Store new data table in a csv file:
write.csv(CREST60X_data_2, "~/Downloads/2024/Modified Drinking in the Dark/CREST60X_data_2_ForGraphs.csv", row.names=FALSE)

########################################################################################################################
# In this section, I am calculating the average value for each treatment and for each sex.

# Select only a specific ROI data:
CREST60X_data_ROI1 <- CREST60X_data[CREST60X_data$ROI == "ROI1", ]

CREST60X_data_ROI2 <- CREST60X_data[CREST60X_data$ROI == "ROI2", ]

CREST60X_data_ROI3 <- CREST60X_data[CREST60X_data$ROI == "ROI3", ]

CREST60X_data_ROI4 <- CREST60X_data[CREST60X_data$ROI == "ROI4", ]

CREST60X_data_ROI5 <- CREST60X_data[CREST60X_data$ROI == "ROI5", ]

CREST60X_data_ROI6 <- CREST60X_data[CREST60X_data$ROI == "ROI6", ]

CREST60X_data_ROI7 <- CREST60X_data[CREST60X_data$ROI == "ROI7", ]

#Use ggpubr to generate graphs of the means after aggregating the data:
library(ggpubr)
personal_theme = theme(plot.title = element_text(hjust = 0.5)) #This makes sure the Title is Center

View(CREST60X_data_ROI1)
# identify the outlier in all the ROIs and all the Slices, for one measurement 
CREST60X_data_ROI1_mean_val <- aggregate (mGFPFfVolume ~ Treatment + Sex, data = CREST60X_data_ROI1, FUN = mean)
CREST60X_data_ROI1_mean_val

CREST60X_data_ROI1_sd_val <- aggregate (mGFPFfVolume ~ Treatment + Sex, data = CREST60X_data_ROI1, FUN = sd)
CREST60X_data_ROI1_sd_val

Outlier <- function(standarddev) {standarddev * 1.5}
ThresholdforOutlier <- aggregate (mGFPFfVolume ~ Treatment + Sex, data = CREST60X_data_ROI1_sd_val, FUN = Outlier)

ApplyThreshold <- function(Thresholdvalue, meanValue, Measurement) {abs(Measurement - meanValue) > Thresholdvalue}

CREST60X_data_ROI1$Outliers <- abs(CREST60X_data_ROI1$mGFPFfVolume - CREST60X_data_ROI1_mean_val) > ThresholdforOutlier
View(CREST60X_data_ROI1)

CREST60X_data_ROI1_Outliers <- CREST60X_data_ROI1 [CREST60X_data_ROI1$Outliers==TRUE, ]
View(CREST60X_data_ROI1_Outliers)

CREST60X_data_ROI1$Outliers <- as.data.frame(CREST60X_data_ROI1$mGFPFfVolume[abs(CREST60X_data_ROI1$mGFPFfVolume - CREST60X_data_ROI1_mean_val) > Thresholdforexclusion])
View(Outliers)

# Aggregate data into a new dataframe.
# This allows the calculation of the mean of a specific measurement per mouse while collapsing across the "Slice" column or Bregma.
mGFPFfVolume.meanpermouseROI1 <- aggregate(mGFPFfVolume ~ CodedMouseID + Sex + Treatment, data = CREST60X_data_ROI1, FUN = mean)

ROI1.PVT <- ggbarplot(mGFPFfVolume.meanpermouseROI1, 
          x = "Treatment", 
          y = "mGFPFfVolume", 
          ylim = c(0, 20000),
          title = "ROI 1: PVT",
          ylab= "mGFP fibers volume", 
          xlab = "Group", 
          fill = "Treatment", 
          facet.by = "Sex",
          error.plot = "upper_errorbar", 
          label = FALSE,
          add = c("mean_se", "jitter"), 
          position = position_dodge(0.8)) + personal_theme + scale_fill_manual(values=c("#ffffff", "#cc0000"))

png("ROI1_PVT.png")
print(ROI1.PVT)
dev.off()

# Make a linear model for statistical analyses:
mGFPFfVolume.Treatment.Sex_model <- lm(mGFPFfVolume ~ Treatment + Sex, data = mGFPFfVolume.meanpermouseROI1)
summary(mGFPFfVolume.Treatment.Sex_model)
anova(mGFPFfVolume.Treatment.Sex_model)   # This is a Two-way ANOVA for ROI1: the PVT
anova(mGFPFfVolume.Treatment.Sex_model) %>% as.data.frame() %>% write.csv(file = "ROI1.STATS_anova_table.csv")


########################################################################################################################

# Aggregate data into a new dataframe[this allows the calculation of the mean of a specific measurement for per mouse collapsing across the "Slice" column which indicates the Bregma]
mGFPFfVolume.meanpermouseROI2 <- aggregate(mGFPFfVolume ~ CodedMouseID + Sex + Treatment, data = CREST60X_data_ROI2, FUN = mean)

ROI2.PVN <- ggbarplot(mGFPFfVolume.meanpermouseROI2, 
                                 x = "Treatment", 
                                 y = "mGFPFfVolume", 
                                 ylim = c(0, 50000),
                                 title = "ROI 2: PVN",
                                 ylab= "mGFP fibers volume", 
                                 xlab = "Group", 
                                 fill = "Treatment", 
                                 facet.by = "Sex",
                                 error.plot = "upper_errorbar", 
                                 label = FALSE,
                                 add = c("mean_se", "jitter"), 
                                 position = position_dodge(0.8)) + personal_theme + scale_fill_manual(values=c("#ffffff", "#cc0000"))

png("ROI2_PVN.png")
print(ROI2.PVN)
dev.off()

# Make a linear model for statistical analyses:
mGFPFfVolume.Treatment.Sex_model <- lm(mGFPFfVolume ~ Treatment + Sex, data = mGFPFfVolume.meanpermouseROI2)
summary(mGFPFfVolume.Treatment.Sex_model)
anova(mGFPFfVolume.Treatment.Sex_model)   # This is a Two-way ANOVA
anova(mGFPFfVolume.Treatment.Sex_model) %>% as.data.frame() %>% write.csv(file = "ROI2.STATS_anova_table.csv")

########################################################################################################################
# Select only data from Slice 3-6, exclude Slice 1 and 2
#Because of the shape of the CeA, we need to make sure we are selecting for the right Bregmas for the CeA, so we excluded slice 1 and 2 for better accuracy in CeA structure.
CREST60X_data_ROI3$Slice != "Slice1"
CREST60X_data_ROI3 <- CREST60X_data_ROI3[CREST60X_data_ROI3$Slice != "Slice1", ]
CREST60X_data_ROI3 <- CREST60X_data_ROI3[CREST60X_data_ROI3$Slice != "Slice2", ]
View(CREST60X_data_ROI3)
# Aggregate data into a new dataframe[this allows the calculation of the mean of a specific measurement for per mouse collapsing across the "Slice" column which indicates the Bregma]
mGFPFfVolume.meanpermouseROI3 <- aggregate(mGFPFfVolume ~ CodedMouseID + Sex + Treatment, data = CREST60X_data_ROI3, FUN = mean)

ROI3.CeA <- ggbarplot(mGFPFfVolume.meanpermouseROI3, 
                      x = "Treatment", 
                      y = "mGFPFfVolume", 
                      ylim = c(0, 30000),
                      title = "ROI 3: CeA",
                      ylab= "mGFP fibers volume", 
                      xlab = "Group", 
                      fill = "Treatment", 
                      facet.by = "Sex",
                      error.plot = "upper_errorbar", 
                      label = FALSE,
                      add = c("mean_se", "jitter"), 
                      position = position_dodge(0.8)) + personal_theme + scale_fill_manual(values=c("#ffffff", "#cc0000"))

png("ROI3_CeA.png")
print(ROI3.CeA)
dev.off()

# Make a linear model for statistical analyses:
mGFPFfVolume.Treatment.Sex_model <- lm(mGFPFfVolume ~ Treatment + Sex, data = mGFPFfVolume.meanpermouseROI3)
summary(mGFPFfVolume.Treatment.Sex_model)
anova(mGFPFfVolume.Treatment.Sex_model)   # This is a Two-way ANOVA
anova(mGFPFfVolume.Treatment.Sex_model) %>% as.data.frame() %>% write.csv(file = "ROI3.STATS_anova_table.csv")

# Store new data table in a csv file:
write.csv(mGFPFfVolume.meanpermouseROI3, "~/Downloads/Richardson Lab/Modified Drinking in the Dark/CREST60_RData_June.11.2024_ForGraphs_ROI3_mGFPVolume.csv", row.names=FALSE)

########################################################################################################################

# Aggregate data into a new dataframe[this allows the calculation of the mean of a specific measurement for per mouse collapsing across the "Slice" column which indicates the Bregma]
mGFPFfVolume.meanpermouseROI4 <- aggregate(mGFPFfVolume ~ CodedMouseID + Sex + Treatment, data = CREST60X_data_ROI4, FUN = mean)

ROI4.Fimbria <- ggbarplot(mGFPFfVolume.meanpermouseROI4, 
                      x = "Treatment", 
                      y = "mGFPFfVolume", 
                      ylim = c(0, 70000),
                      title = "ROI 4: Fimbria",
                      ylab= "mGFP fibers volume", 
                      xlab = "Group", 
                      fill = "Treatment", 
                      facet.by = "Sex",
                      error.plot = "upper_errorbar", 
                      label = FALSE,
                      add = c("mean_se", "jitter"), 
                      position = position_dodge(0.8)) + personal_theme + scale_fill_manual(values=c("#ffffff", "#cc0000"))

png("ROI4_Fimbria.png")
print(ROI4.Fimbria)
dev.off()

# Make a linear model for statistical analyses:
mGFPFfVolume.Treatment.Sex_model <- lm(mGFPFfVolume ~ Treatment + Sex, data = mGFPFfVolume.meanpermouseROI4)
summary(mGFPFfVolume.Treatment.Sex_model)
anova(mGFPFfVolume.Treatment.Sex_model)   # This is a Two-way ANOVA
anova(mGFPFfVolume.Treatment.Sex_model) %>% as.data.frame() %>% write.csv(file = "ROI4.STATS_anova_table.csv")

########################################################################################################################

# Aggregate data into a new dataframe[this allows the calculation of the mean of a specific measurement for per mouse collapsing across the "Slice" column which indicates the Bregma]
mGFPFfVolume.meanpermouseROI5 <- aggregate(mGFPFfVolume ~ CodedMouseID + Sex + Treatment, data = CREST60X_data_ROI5, FUN = mean)

ROI5.MedialHabenula <- ggbarplot(mGFPFfVolume.meanpermouseROI5, 
                          x = "Treatment", 
                          y = "mGFPFfVolume", 
                          ylim = c(0, 50000),
                          title = "ROI 5: MedialHabenula",
                          ylab= "mGFP fibers volume", 
                          xlab = "Group", 
                          fill = "Treatment", 
                          facet.by = "Sex",
                          error.plot = "upper_errorbar", 
                          label = FALSE,
                          add = c("mean_se", "jitter"), 
                          position = position_dodge(0.8)) + personal_theme + scale_fill_manual(values=c("#ffffff", "#cc0000"))

png("ROI5_MedialHabenula.png")
print(ROI5.MedialHabenula)
dev.off()

# Make a linear model for statistical analyses:
mGFPFfVolume.Treatment.Sex_model <- lm(mGFPFfVolume ~ Treatment + Sex, data = mGFPFfVolume.meanpermouseROI5)
summary(mGFPFfVolume.Treatment.Sex_model)
anova(mGFPFfVolume.Treatment.Sex_model)   # This is a Two-way ANOVA
anova(mGFPFfVolume.Treatment.Sex_model) %>% as.data.frame() %>% write.csv(file = "ROI5.STATS_anova_table.csv")

########################################################################################################################
# Select only data from Slice 1 and 2, exclude Slice 3-6
#Because of the shape of the CeA, we need to make sure we are selecting for the right Bregmas for the CeA, so we excluded slice 1 and 2 for better accuracy in CeA structure.
CREST60X_data_ROI6$Slice != "Slice3"
CREST60X_data_ROI6 <- CREST60X_data_ROI6[CREST60X_data_ROI6$Slice != "Slice3", ]
CREST60X_data_ROI6 <- CREST60X_data_ROI6[CREST60X_data_ROI6$Slice != "Slice4", ]
CREST60X_data_ROI6 <- CREST60X_data_ROI6[CREST60X_data_ROI6$Slice != "Slice5", ]
CREST60X_data_ROI6 <- CREST60X_data_ROI6[CREST60X_data_ROI6$Slice != "Slice6", ]
View(CREST60X_data_ROI6)

# Aggregate data into a new dataframe[this allows the calculation of the mean of a specific measurement for per mouse collapsing across the "Slice" column which indicates the Bregma]
mGFPFfVolume.meanpermouseROI6 <- aggregate(mGFPFfVolume ~ CodedMouseID + Sex + Treatment, data = CREST60X_data_ROI6, FUN = mean)
View(mGFPFfVolume.meanpermouseROI6)

ROI6.CorpusCallosum <- ggbarplot(mGFPFfVolume.meanpermouseROI6, 
                                 x = "Treatment", 
                                 y = "mGFPFfVolume", 
                                 ylim = c(0, 50000),
                                 title = "ROI 6: Corpus Callosum",
                                 ylab= "mGFP fibers volume", 
                                 xlab = "Group", 
                                 fill = "Treatment", 
                                 facet.by = "Sex",
                                 error.plot = "upper_errorbar", 
                                 label = FALSE,
                                 add = c("mean_se", "jitter"), 
                                 position = position_dodge(0.8)) + personal_theme + scale_fill_manual(values=c("#ffffff", "#cc0000"))

png("ROI6_CorpusCallosum.png")
print(ROI6.CorpusCallosum)
dev.off()

# Make a linear model for statistical analyses:
mGFPFfVolume.Treatment.Sex_model <- lm(mGFPFfVolume ~ Treatment + Sex, data = mGFPFfVolume.meanpermouseROI6)
summary(mGFPFfVolume.Treatment.Sex_model)
anova(mGFPFfVolume.Treatment.Sex_model)   # This is a Two-way ANOVA
anova(mGFPFfVolume.Treatment.Sex_model) %>% as.data.frame() %>% write.csv(file = "ROI6.STATS_anova_table.csv")


# Store new data table in a csv file:
write.csv(mGFPFfVolume.meanpermouseROI6, "~/Downloads/Richardson Lab/Modified Drinking in the Dark/CREST60_RData_June.11.2024_ForGraphs_ROI6_mGFPVolume.csv", row.names=FALSE)

########################################################################################################################
# Select only data from Slice 1 and 2, exclude Slice 3-6
#Because of the shape of the CeA, we need to make sure we are selecting for the right Bregmas for the CeA, so we excluded slice 1 and 2 for better accuracy in CeA structure.
CREST60X_data_ROI7$Slice != "Slice3"
CREST60X_data_ROI7 <- CREST60X_data_ROI7[CREST60X_data_ROI7$Slice != "Slice3", ]
CREST60X_data_ROI7 <- CREST60X_data_ROI7[CREST60X_data_ROI7$Slice != "Slice4", ]
CREST60X_data_ROI7 <- CREST60X_data_ROI7[CREST60X_data_ROI7$Slice != "Slice5", ]
CREST60X_data_ROI7 <- CREST60X_data_ROI7[CREST60X_data_ROI7$Slice != "Slice6", ]
View(CREST60X_data_ROI7)
# Aggregate data into a new dataframe[this allows the calculation of the mean of a specific measurement for per mouse collapsing across the "Slice" column which indicates the Bregma]
mGFPFfVolume.meanpermouseROI7 <- aggregate(mGFPFfVolume ~ CodedMouseID + Sex + Treatment, data = CREST60X_data_ROI7, FUN = mean)
View(mGFPFfVolume.meanpermouseROI7)

ROI7.CingulumBundle <- ggbarplot(mGFPFfVolume.meanpermouseROI7, 
                      x = "Treatment", 
                      y = "mGFPFfVolume", 
                      ylim = c(0, 50000),
                      title = "ROI 7: Cingulum Bundle",
                      ylab= "mGFP fibers volume", 
                      xlab = "Group", 
                      fill = "Treatment", 
                      facet.by = "Sex",
                      error.plot = "upper_errorbar", 
                      label = FALSE,
                      add = c("mean_se", "jitter"), 
                      position = position_dodge(0.8)) + personal_theme + scale_fill_manual(values=c("#ffffff", "#cc0000"))

png("ROI7_Cingulum Bundle.png")
print(ROI7.CingulumBundle)
dev.off()

# Make a linear model for statistical analyses:
mGFPFfVolume.Treatment.Sex_model <- lm(mGFPFfVolume ~ Treatment + Sex, data = mGFPFfVolume.meanpermouseROI7)
summary(mGFPFfVolume.Treatment.Sex_model)
anova(mGFPFfVolume.Treatment.Sex_model)   # This is a Two-way ANOVA
anova(mGFPFfVolume.Treatment.Sex_model) %>% as.data.frame() %>% write.csv(file = "ROI7.STATS_anova_table.csv")

# Store new data table in a csv file:
write.csv(mGFPFfVolume.meanpermouseROI7, "~/Downloads/Richardson Lab/Modified Drinking in the Dark/CREST60_RData_June.11.2024_ForGraphs_ROI7_mGFPVolume.csv", row.names=FALSE)

#############################################################################################################

# Transforming data: [In this case, data does not have to be transformed from wide to long]
#require(reshape2) # data for brain measurements is wide, I will transform it to be long: 
# Use the "melt(dataframename)" function to transform the data if needed.
#CREST60X_data_ROI1_long <- melt(CREST60X_data_ROI1)        

####################################################################################################################
