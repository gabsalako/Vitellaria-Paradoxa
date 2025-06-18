#load required libraries (ensure the packages had been pre-installed)
library(raster)#Work with raster data
library(terra) #Work with raster data
library(rgdal) #Export GeoTIFFs and other core GIS functions
library (randomForest)#for random forest
library(rpart)
library(RStoolbox)
library(caret)#for confusion matrix
library(sp)	
library(spThin)#for data thinning
library(sf)   #GIS vector data analysis
library(corrplot)	#Correlation plot
library(cluster)
library(biomod2)
library(caret)
library(usdm)#for vif
library(data.table) #data management
library(ggplot2)    #plot and graph
library(ggpubr)# for ggpar,ggarrange
library(dismo)#for GBM dismo
library(neuralnet)#for artificial neural network (ANN)
library(XLS)#write to excel
library(xlsx)#write to excel

#=========================================================================
##Set working directories (WDr)
#=========================================================================
setwd("C:/Users/gsalako/Documents/Vit_Model/Vit_draft_Pub/Vitellaria_P Package")
#Check files to ensure your data is in the WDr
list.files()
#=========================================================================
##Load downloaded csv/excel file and check the data summary
#=========================================================================
Vit_data <- read.csv("Vitellaria_P-allGBIF.csv")
summary(Vit_data)
head(Vit_data, 5)
str(Vit_data)
#show coordinates with duplicate records
dupsl <- duplicated(Vit_data[, c('decimalLongitude', 'decimalLatitude')])#duplicate coordinates
sum(dupsl)
#retain only the non duplicated
ClVitdata <- Vit_data[!dupsl, ]
summary(ClVitdata)
#=========================================================================
##convert retained non duplicated data to point data
#=========================================================================
proj4Str <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
Vit_ThinPoints <- SpatialPointsDataFrame(coords = ClVitdata[,c("decimalLongitude","decimalLatitude")], data = ClVitdata, proj4string = CRS(proj4Str))
plot(Vit_ThinPoints,  pch = 16, col ="red")
#=========================================================================
##load study area shapefile and crop the occurence point to the extent
#=========================================================================
#Vitellaria paradoxa study area extent shapefile
VitExt <- readOGR("VitSelVect.shp")
#intersect or crop to the study area extent
Vitextc <- intersect(Vit_ThinPoints, VitExt)
plot(Vitextc, col="green", pch=16, add=TRUE)
plot(VitExt, add=TRUE)
head(Vitextc, 10)
summary(Vitextc)
#remove any unwanted column e.g.
Vitextc$GID_2 <- NULL
Vitextc$GID_0 <- NULL
Vitextc$COUNTRY  <- NULL
Vitextc$GID_1   <- NULL
Vitextc$NAME_1 <- NULL
Vitextc$NL_NAME_1 <- NULL
Vitextc$NAME_2 <- NULL
Vitextc$VARNAME_2 <- NULL
Vitextc$ NL_NAME_2 <- NULL
Vitextc$TYPE_2 <- NULL
Vitextc$ENGTYPE_2 <- NULL
Vitextc$CC_2<- NULL
Vitextc$HASC_2<- NULL
Vitextc$optional<- NULL
head(Vitextc, 10)
#...................
#...................
#=========================================================================
##load the raster files (bands & climate)
#=========================================================================

#WGS raster load ready for use
MB3 <- raster("GrBand3.tif")
MB4 <- raster("ReBand4.tif")
MB5 <- raster("NIRBand5.tif")
MB7 <- raster("SWBand7.tif")
MBwet <- raster("Wet.tif")/10 #precipitation data are divided by 10
MBdry <- raster("dry.tif")/10
MBrain <- raster("Precip2.tif")
#plot sample raster
plot(MBwet)
#=========================================================================
##Arithmetic operation of vegetation index
#=========================================================================

##plot NDVI, green cover/biomass
NDVI <- (MB5 - MB4)/(MB5 +MB4)
#RatioVI similar to NDVI but less sensitive to soil effect
RVI <- (MB5 - MB4)/sqrt(MB5 +MB4)#remove for multicollineraity problem
#GreeNDVI, chlorophyll
GNDVI <- (MB5 - MB3)/(MB5 +MB3)
#two band Enhance Vegetation Index2, G=2.5
EVI2 <- 2.5 * (MB5 - MB4)/(MB5 + 2.4 * MB4+ 1)
#SAVI, L=0.5
SAVI <- (MB5 - MB4) * (1 + 0.5)/(MB5 + MB4 + 0.5)
#specific leaf area index (SLAVI)
SLAVI <- MB5/(MB4 + MB7)
#Normalized water index(NDWI), water stress
NDWI <-  (MB3 - MB5)/(MB3 + MB5)
#Modified Normalized water index(NDWI), water stress
MNDWI <-  (MB3 - MB7)/(MB3 + MB7)

#=========================================================================
##stack all vegetation index and climate data
#=========================================================================

#all, veg climate and SDM
StackIndex=stack(NDVI,GNDVI,  EVI2, SAVI,  NDWI, MNDWI, MBwet,MBdry,MBrain)
#rename the raster
names(StackIndex) <-c("NDVI","GNDVI","SAVI", "NDWI", "MNDWI","EVI", "Wet", "Dry", "PPt")

#=========================================================================
#sample random background data & create data frame
#=========================================================================
#sample random background data & data frame
presvals <- extract(StackIndex, Vitextc)#stack raster & the occ.converted to Spatialpoint as your presence point
backgr <- randomPoints(StackIndex, 1000) #obtain random background point, load dismo for randompoint
absvals <- extract(StackIndex, backgr)#now your absence points
Vit_pb1 <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals))) #line up presence and absence points
vitDF <- data.frame(cbind(Vit_pb1, rbind(presvals, absvals)))#make a dataframe by row and column bind
vitDF <- na.omit(vitDF)#Remove NA completely
summary(vitDF)
head(vitDF, 10)
tail(vitDF, 10)

#=========================================================================
#Model fitting/calibration
#=========================================================================
#Model fitting BRT
GBM_Vit <- gbm.step(data=vitDF, gbm.x = 2:10, gbm.y = 1, #x= environmental variable column, y=species or response varible column
                    family = "gaussian", tree.complexity = 5,
                    learning.rate = 0.001, bag.fraction = 0.5) #trees adding n.trees makes difference;look good use for reference values
#response curve
par(mfrow = c(1,1))
gbm.plot(GBM_Vit, write.title = TRUE)#group plot
summary(GBM_Vit)#variable importance
#prediction
Pred_Vit <- predict(StackIndex, GBM_Vit)
plot(Pred_Vit)
writeRaster(Pred_Vit, filename = "Vit_SDM.tif", format='GTiff', overwrite=TRUE)#convert to ASCII or GTiff to tif or CDF for netcdf
#Model fitting RF
RF_vit <- randomForest(x= vitDF[,2:10], y= vitDF[,1], ntree = 1000,nodesize = 10, importance = T)
varImpPlot(RF_vit)
PredRF_vit<- predict(StackIndex, RF_vit)
plot(PredRF_vit)
#partial plot
par(mfrow = c(3,3))
partialPlot(RF_vit, vitDF, "PPt", plot = TRUE)
partialPlot(RF_vit, vitDF, "Wet", plot = TRUE)
partialPlot(RF_vit, vitDF, "Dry", plot = TRUE)
partialPlot(RF_vit, vitDF, "MNDWI", plot = TRUE)
partialPlot(RF_vit, vitDF, "GNDVI", plot = TRUE)
partialPlot(RF_vit, vitDF, "NDWI", plot = TRUE)
partialPlot(RF_vit, vitDF, "SAVI", plot = TRUE)
partialPlot(RF_vit, vitDF, "NDVI", plot = TRUE)
partialPlot(RF_vit, vitDF, "EVI", plot = TRUE)

##Model fitting ANN
#Nueralnet may have to be scaled or normalised
scaleddata<-scale(vitDF)
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
NNVit = neuralnet(
  Vit_pb1~.,
  data=scaleddata,
  hidden=c(4,2),
  rep=10,
  linear.output = FALSE
)
#plot model fit
plot(NNVit,rep = 10)
#prediction
Pred_NNVit <- predict(StackIndex, NNVit)
plot(Pred_NNVit)# poor prediction

#=========================================================================
#Model evaluation/assessment
#=========================================================================
#Split data
set.seed(245)
data_rows <- floor(0.80 * nrow(vitDF))
train_indices <- sample(c(1:nrow(vitDF)), data_rows)
train_data <- vitDF[train_indices,]
test_data <- vitDF[-train_indices,]
#Evaluation
#Confusion matrix, calculating Kappa, accuracy (RF)
Pred_test_vit <- predict(RF_vit, newdata =  test_data)#train model and test data to predict using test data
predicted_class_RF <- as.factor(ifelse(Pred_test_vit >= 0.5, 
                                       "1", "0")) #presence(1) absence(0)

confusionMatrix(predicted_class_RF, as.factor(test_data$Vit_pb1))#predicted classes, train test data$species to evaluate as in dataframe column

#AUC_RF
AUC_RFvit<- Metrics::auc(test_data$Vit_pb1, Pred_test_vit)#test data, and pred of train model and test data 
AUC_RFvit
0.9


#Confusion matrix, calculating Kappa, accuracy (GBM)
Pred_test_vit2 <- predict(GBM_Vit, newdata =  test_data)#train model, and new train_test data to predict using test data
predicted_class_GBM <- as.factor(ifelse(Pred_test_vit2 >= 0.5, 
                                        "1", "0")) #presence(1) absence(0)

confusionMatrix(predicted_class_GBM, as.factor(test_data$Vit_pb1))#predicted classes, train test data$species to evaluate as in dataframe column

##Confusion matrix, calculating Kappa, accuracy (ANN)
Pred_test_vit3 <- predict(NNVit, newdata =  test_data)#train model, and new train_test data to predict using test data
predicted_class_ANN <- as.factor(ifelse(Pred_test_vit3 >= 0.5, 
                                        "1", "0")) #presence(1) absence(0)

confusionMatrix(predicted_class_ANN, as.factor(test_data$Vit_pb1))#predicted classes, train test data$species to evaluate as in dataframe column

#AUC_ANN
AUC_ANNvit<- Metrics::auc(test_data$Vit_pb1, Pred_test_vit3)#test data, and pred of train model and test data 
AUC_ANNvit
#=========================================================================
#Averaging the successful models
#=========================================================================
#Mean BRT RF
EnsVit <- mean(PredRF_vit, Pred_Vit)
plot(EnsVit)
#mapping uncertainty or coefficent of variation
stackvit=stack(Pred_Vit, PredRF_vit)
sdVit <- calc(stackvit, fun=sd)#standard deviation
plot(sdVit)
Vit_cv2 <- (sdVit/EnsVit)*10
plot(Vit_cv2)


##=========================================================================
#Thresholding the vegetation index
#=========================================================================
#Calculating modal or peak of NDVI
# create model function first
mode <- function(x, na.rm = FALSE) {
  
  if(na.rm){ #if na.rm is TRUE, remove NA values from input x
    x = x[!is.na(x)]
  }
  
  val <- unique(x)
  return(val[which.max(tabulate(match(x, val)))])
}
#apply it to data frame
#Mode NDVI
ModeNDVI <- mode(vitDF2$NDVI)
ModeNDVI
modeDry <-mode(vitDF2$Dry)
modeDry
ModeGNDVI <- mode(vitDF2$GNDVI)
ModeGNDVI
ModeSAVI <- mode(vitDF2$SAVI)
ModeSAVI
ModeMNDWI <- mode(vitDF2$MNDWI)
ModeMNDWI
ModeNDWI <- mode(vitDF2$NDWI)
ModeNDWI
ModeEVI <-mode(vitDF2$EVI)
ModeEVI

