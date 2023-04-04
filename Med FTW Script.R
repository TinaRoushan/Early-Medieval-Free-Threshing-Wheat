library('Momocs')
library(geomorph)
#Repeatability and measurement error
MTFTWCoords <- list.files("D:\\Coordinates\\Free Threshing\\FTW _Method_Test", full.names = TRUE)
MTFTWFrame <- read.csv("D:\\Coordinate_Frames\\MethodTest_FTW.csv", header = TRUE)
MTFTWTxt <- import_txt(MTFTWCoords, fileEncoding="UTF-8-BOM")#fileencoding gets rid of weird symbols
MTFTWOut <- Out(MTFTWTxt, fac=MTFTWFrame) #creates an Out object from a specified list- you can specify landmarks in "ldk"
MTFTWOut.l <- filter(MTFTWOut, View == "l")#creates subset of lateral views only
MTFTWOut.d <- filter(MTFTWOut, View == "d")#creates subset of lateral views only
MTFTWOut.l.efour <- efourier(MTFTWOut.l, nb.h=8, norm = FALSE, start = FALSE)
MTFTWOut.d.efour <- efourier(MTFTWOut.d, nb.h=8, norm = FALSE, start = FALSE)
MTFTWOut.d.l <- combine (MTFTWOut.d, MTFTWOut.l)
#Procrustes ANOVA
export(MTFTWOut.l.efour)
export(MTFTWOut.d.efour)
MTFTW_l <- read.csv("D:\\Coordinate_Frames\\MTFTW_l.csv")
MTFTW_d <- read.csv("D:\\Coordinate_Frames\\MTFTW_d.csv")
#Calculate PCA scores for above
PCAObservation_d<-prcomp(MTFTW_d)$x 
PCAObservation_l<-prcomp(MTFTW_l)$x 
#Creates 5 level session factor (consecutive)
sessionfactor<-as.factor(gl(5, 5))
sessionfactor
#Create 5 level individual factor (repeated)
individualfactor<-as.factor(rep(1:5, 5))
individualfactor
#Geomorph data frames for PCA computations of coefficients
gdf_d<-geomorph.data.frame(shape=PCAObservation_d) #ID: a vector containing the specimens ID 
gdf_l<-geomorph.data.frame(shape=PCAObservation_l) #ID: a vector containing the specimens ID 
#Procrustes ANOVA using session factor and then individual factor as dependent variables
summary(procD.lm(shape~sessionfactor, data=gdf_d))
summary(procD.lm(shape~sessionfactor, data=gdf_l))
#=session difference not significant for dorsal or lateral view
mod<-summary(procD.lm(shape~individualfactor, data= gdf_d))
mod2<-summary(procD.lm(shape~individualfactor, data= gdf_l))
mod
mod2
#=indiv difference significant for dorsal and lateral view
#dorsal measurement error (See Claude 2008)
s2within<-mswithin<-mod[[1]][2,3]
mod[[1]][2, 3]
MSamong<-mod[[1]][1, 3]
MSamong
s2among<-(MSamong-mswithin)/5
s2within/(s2within+s2among)*100#measurement error
#lateral measurement error
s2within<-mswithin<-mod2[[1]][2,3]
mod2[[1]][2, 3]
MSamong<-mod2[[1]][1, 3]
MSamong
s2among<-(MSamong-mswithin)/5
s2within/(s2within+s2among)*100#measurement error


#Stage 1 Uncharred
FTWRefCoords <- list.files("D:\\Coordinates\\Free Threshing\\Uncharred.01", full.names = TRUE)
FTWRefCoords
FTWRefFrame <- read.csv("D:\\Coordinate_Frames\\Uncharred_Frame_FT.csv", header = TRUE)
FTWRefTxt <- import_txt(FTWRefCoords, fileEncoding="UTF-8-BOM")#fileencoding gets rid of weird symbols
FTWRefOut <- Out(FTWRefTxt, fac=FTWRefFrame) #creates an Out object from a specified list- you can specify landmarks in "ldk"
FTWRefOut2 <-coo_scale (FTWRefOut)
FTWRefOut3 <-coo_center (FTWRefOut2)
###Filter and Fourier transform
FTWRefOut.l <- filter(FTWRefOut3, View == "l")#creates subset of lateral views only
FTWRefOut.d <- filter(FTWRefOut3, View == "d")#creates subset of dorsal views only
FTWRefOut.p <- filter(FTWRefOut3, View == "p")#creates subset of polar views only
calibrate_harmonicpower_efourier(FTWRefOut.l,nb.h=12)#tells you how many harmonics will be needed to gather x% of harmonic power- here 8= 99%
FTWRefOut.l.efour <- efourier(FTWRefOut.l, nb.h=8, norm = FALSE, start = TRUE)
FTWRefOut.d.efour <- efourier(FTWRefOut.d, nb.h=8, norm = FALSE, start = TRUE)
FTWRefOut.p.efour <- efourier(FTWRefOut.p, nb.h=8, norm = FALSE, start = TRUE)
FTWRefOut.d.l <- combine (FTWRefOut.d, FTWRefOut.l)#dataset of just dorsal and lateral views
FTWRefOut.d.l.p <- combine (FTWRefOut.d, FTWRefOut.l, FTWRefOut.p)#dataset of just dorsal and lateral views
#LDA
#LDA 
FTWRefOut.d %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
FTWRefOut.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
FTWRefOut.p %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
FTWRefOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
FTWRefOut.d.l.p %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 


#Charred reference material 
FTWRefCoords <- list.files("D:\\Coordinates\\Free Threshing\\Charred.02", full.names = TRUE)
FTWRefCoords
FTWRefFrame <- read.csv("D:\\Coordinate_Frames\\Charred_Frame_FT.csv", header = TRUE)
FTWRefTxt <- import_txt(FTWRefCoords, fileEncoding="UTF-8-BOM")#fileencoding gets rid of weird symbols
FTWRefOut <- Out(FTWRefTxt, fac=FTWRefFrame) #creates an Out object from a specified list- you can specify landmarks in "ldk"
FTWRefOut2 <-coo_scale (FTWRefOut)
FTWRefOut3 <-coo_center (FTWRefOut2)
###Filter and Fourier transform
FTWRefOut.l <- filter(FTWRefOut3, View == "l")#creates subset of lateral views only
FTWRefOut.d <- filter(FTWRefOut3, View == "d")#creates subset of dorsal views only
FTWRefOut.p <- filter(FTWRefOut3, View == "p")#creates subset of polar views only
FTWRefOut.l.efour <- efourier(FTWRefOut.l, nb.h=8, norm = FALSE, start = TRUE)
FTWRefOut.d.efour <- efourier(FTWRefOut.d, nb.h=8, norm = FALSE, start = TRUE)
FTWRefOut.p.efour <- efourier(FTWRefOut.p, nb.h=8, norm = FALSE, start = TRUE)
FTWRefOut.d.l <- combine (FTWRefOut.d, FTWRefOut.l)#dataset of just dorsal and lateral views
FTWRefOut.d.l.p <- combine (FTWRefOut.d, FTWRefOut.l, FTWRefOut.p)#dataset of just dorsal and lateral views
#LDA
#LDA 
FTWRefOut.d %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
FTWRefOut.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
FTWRefOut.p %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
FTWRefOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
FTWRefOut.d.l.p %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 

##Reference with compactum
FTWRefCoords <- list.files("D:\\Coordinates\\Free Threshing\\Charred_Compactum", full.names = TRUE)
FTWRefFrame <- read.csv("D:\\Coordinate_Frames\\Charred_Frame_FT2.csv", header = TRUE)
FTWRefTxt <- import_txt(FTWRefCoords, fileEncoding="UTF-8-BOM")#fileencoding gets rid of weird symbols
FTWRefOut <- Out(FTWRefTxt, fac=FTWRefFrame) #creates an Out object from a specified list- you can specify landmarks in "ldk"
FTWRefOut2 <-coo_scale (FTWRefOut)
FTWRefOut3 <-coo_center (FTWRefOut2)
###Filter and Fourier transform
FTWRefOut.l <- filter(FTWRefOut3, View == "l")#creates subset of lateral views only
FTWRefOut.d <- filter(FTWRefOut3, View == "d")#creates subset of dorsal views only
FTWRefOut.l.efour <- efourier(FTWRefOut.l, nb.h=8, norm = FALSE, start = TRUE)
FTWRefOut.d.efour <- efourier(FTWRefOut.d, nb.h=8, norm = FALSE, start = TRUE)
FTWRefOut.p.efour <- efourier(FTWRefOut.p, nb.h=8, norm = FALSE, start = TRUE)
FTWRefOut.d.l <- combine (FTWRefOut.d, FTWRefOut.l)#dataset of just dorsal and lateral views
#LDA 
FTWRefOut.d %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
FTWRefOut.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
FTWRefOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 

##Final dataset
FTWRefCoords <- list.files("D:\\Coordinates\\Free Threshing\\FeedSax_Ref", full.names = TRUE)
FTWRefFrame <- read.csv("D:\\Coordinate_Frames\\Charred_Frame_FeedSax_Ref.csv", header = TRUE)
FTWRefTxt <- import_txt(FTWRefCoords, fileEncoding="UTF-8-BOM")#fileencoding gets rid of weird symbols
FTWRefOut <- Out(FTWRefTxt, fac=FTWRefFrame) #creates an Out object from a specified list- you can specify landmarks in "ldk"
FTWRefOut2 <-coo_scale (FTWRefOut)
FTWRefOut3 <-coo_center (FTWRefOut2)
###Filter and Fourier transform
FTWRefOut.l <- filter(FTWRefOut3, View == "l")#creates subset of lateral views only
FTWRefOut.d <- filter(FTWRefOut3, View == "d")#creates subset of dorsal views only
FTWRefOut.l.efour <- efourier(FTWRefOut.l, nb.h=8, norm = FALSE, start = TRUE)
FTWRefOut.d.efour <- efourier(FTWRefOut.d, nb.h=8, norm = FALSE, start = TRUE)
FTWRefOut.d.l <- combine (FTWRefOut.d, FTWRefOut.l)#dataset of just dorsal and lateral views
#LDA 
FTWRefOut.d %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
FTWRefOut.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
FTWRefOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
#By accession
FTWRefOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('Accession', scale=FALSE, center= TRUE) 
#Plot
FTWRefOut.d.l.lda <- FTWRefOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE, center= TRUE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
plot(FTWRefOut.d.l.lda, points= TRUE , labelspoints= FALSE, zoom= 1.4, cex.labelsgroups= 0.8, rect.labelsgroups= TRUE, cex = 0.1)


###As above but with germinated grains 
FTWRefCoordsG <- list.files("D:\\Coordinates\\Free Threshing\\FeedSax_Ref_Germ", full.names = TRUE)
FTWRefFrameG <- read.csv("D:\\Coordinate_Frames\\Charred_Frame_FeedSax_Ref_Germ.csv", header = TRUE)
FTWRefTxtG <- import_txt(FTWRefCoordsG, fileEncoding="UTF-8-BOM")#fileencoding gets rid of weird symbols
FTWRefOutG <- Out(FTWRefTxtG, fac=FTWRefFrameG) #creates an Out object from a specified list- you can specify landmarks in "ldk"
FTWRefOutG2 <-coo_scale (FTWRefOutG)
FTWRefOutG.l <- filter(FTWRefOutG2, View == "l")#creates subset of lateral views only
FTWRefOutG.d <- filter(FTWRefOutG2, View == "d")#creates subset of dorsal views only
FTWRefOutG.l.efour <- efourier(FTWRefOutG.l, nb.h=8, norm = FALSE, start = TRUE)
FTWRefOutG.d.efour <- efourier(FTWRefOutG.d, nb.h=8, norm = FALSE, start = TRUE)
FTWRefOutG.d.l <- combine (FTWRefOutG.d, FTWRefOutG.l)#dataset of just dorsal and lateral views
FTWRefOutG.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
FTWRefOutG.d.l.lda <- FTWRefOutG.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE, center= TRUE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE, prior= c(1,1,1,1,1)/5)
plot(FTWRefOutG.d.l.lda, points= TRUE , labelspoints= FALSE, zoom= 1.2, cex.labelsgroups= 0.8, rect.labelsgroups= TRUE, cex = 0.1)


#3d plot for Part 1
library(rgl)
gcolours<-c("#d90602","#55059c","#20a6b3","#62b320")
open3d()
par3d(windowRect = c(100, 100, 612, 612))
FTWRefOut.d.l.lda.colour<-gcolours[as.factor(FTWRefOut.d.l.lda$fac)]
LD<-as.data.frame(FTWRefOut.d.l.lda[["mod.pred"]][["x"]])
LD
plot3d(LD$LD1, LD$LD2, LD$LD3,col=FTWRefOut.d.l.lda.colour, xlab= "LD1", ylab="LD2", zlab="LD3")
rgl.snapshot("modernaccessions_FeedSax.png","png")
?rgl.snapshot
rgl.postscript("modernaccessions_FeedSax.eps","eps")

#Comparison with arch
#Unknown LDA reclass- select dataset for reclassification from FeedSax_All- example uses Houghton
UnclassCoords <- list.files("D:\\Coordinates\\Free Threshing\\Unclassified", full.names = TRUE)#select data for reclassification from 'FeedSax_All'
UnclassFrame <- read.csv("D:\\Coordinate_Frames\\Unclassified3.csv", header = TRUE)#select data for reclassification from 'FeedSax_All'
UnclassTxt <- import_txt(UnclassCoords, fileEncoding="UTF-8-BOM")#fileencoding gets rid of weird symbols
UnclassOut <- Out(UnclassTxt, UnclassFrame) #creates an Out object from a specified list- you can specify landmarks in "ldk"
UnclassOut %>% coo_center %>% coo_scale %>% stack # centres --> scales (normalises) ---> stacks your Out objects by centroid size (?check this is right)
UnclassOut1 <-coo_scale(UnclassOut)
UnclassOut2 <-coo_scale(UnclassOut1)
#reclass
UnclassOut.d <- filter(UnclassOut2, View == "d")
UnclassOut.l <- filter(UnclassOut2, View == "l")
UnclassOut.d.l <- combine (UnclassOut.d, UnclassOut.l)
UnclassOut.d.l.efour <- UnclassOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine
reLDA(UnclassOut.d.l.efour, FTWRefOut.d.l.lda)

#with dataset including experimentally germinated grains 
reLDA(UnclassOut.d.l.efour, FTWRefOutG.d.l.lda)
#with priors modified for even probability between groups 
mod<-FTWRefOutG.d.l.lda$mod
newdata<-UnclassOut.d.l.efour
predicted<-predict(mod, newdata$coe, prior= c(1,1,1,1,1)/5)
predicted

######## Arch comparisons e.g. Houghton 
FTWRefCoordsG <- list.files("D:\\Coordinates\\Free Threshing\\FeedSax_Ref_Comparison", full.names = TRUE)
FTWRefFrameG <- read.csv("D:\\Coordinate_Frames\\Charred_Frame_FeedSax_Comparison.csv", header = TRUE)
FTWRefTxtG <- import_txt(FTWRefCoordsG, fileEncoding="UTF-8-BOM")#fileencoding gets rid of weird symbols
FTWRefOutG <- Out(FTWRefTxtG, fac=FTWRefFrameG) #creates an Out object from a specified list- you can specify landmarks in "ldk"
FTWRefOutG2 <-coo_scale (FTWRefOutG)
FTWRefOutG3 <-coo_scale (FTWRefOutG2)
FTWRefOutG.l <- filter(FTWRefOutG3, View == "l")#creates subset of lateral views only
FTWRefOutG.d <- filter(FTWRefOutG3, View == "d")#creates subset of dorsal views only
FTWRefOutG.l.efour <- efourier(FTWRefOutG.l, nb.h=8, norm = FALSE, start = TRUE)
FTWRefOutG.d.efour <- efourier(FTWRefOutG.d, nb.h=8, norm = FALSE, start = TRUE)
FTWRefOutG.d.l <- combine (FTWRefOutG.d, FTWRefOutG.l)#dataset of just dorsal and lateral views
FTWRefOutG.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
FTWRefOutG.d.l.lda <- FTWRefOutG.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE, center= TRUE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE, prior= c(1,1,1,1,1)/5)
plot(FTWRefOutG.d.l.lda, points= TRUE , labelspoints= FALSE, zoom= 1.2, cex.labelsgroups= 0.8, rect.labelsgroups= TRUE, cex = 0.1)



#visualise e.g. Rampton Rivet stacked and centered
#dorsal
FTWRefOut.d.site <-filter (FTWRefOut.d, Accession == "W0508")
#centre, scale by centroid size and stack outlines
FTWRefOut.d.site %>% coo_center %>% coo_scale %>% stack
#lateral
FTWRefOut.l.site <-filter (FTWRefOut.l, Accession == "W0508")
#centre, scale by centroid size and stack outlines
FTWRefOut.l.site %>% coo_center %>% coo_scale %>% stack

#Mean shapes
#MANOVA/Mean shapes
MeanDorsal<-MSHAPES(FTWRefOut.l.efour, fac = 'taxon.code', FUN = mean, nb.pts = 120)
MeanDorsal2<-MSHAPES(FTWRefOut.d.efour, fac = 'taxon.code', FUN = mean, nb.pts = 120)
FTW.ms <- MSHAPES(FTWRefOut.d.efour, 'taxon.code')
FTW.ms2<- MSHAPES(FTWRefOut.l.efour, 'taxon.code')
#pairwise comparison
plot_mshapes (FTW.ms, size= 3/4)
plot_mshapes (FTW.ms2, size= 3/4)



#All arch data
FeedSaxCoords <- list.files("D:\\Coordinates\\Free Threshing\\FeedSax_All", full.names = TRUE)
FeedSaxFrame <- read.csv("D:\\Coordinate_Frames\\FeedSax_All.csv", header = TRUE)
FeedSaxtxt <- import_txt(FeedSaxCoords, fileEncoding="UTF-8-BOM")#fileencoding gets rid of weird symbols
FeedSaxOut <- Out(FeedSaxtxt, fac=FeedSaxFrame) #creates an Out object from a specified list- you can specify landmarks in "ldk"
FeedSaxOut2 <-coo_scale (FeedSaxOut)
FeedSaxOut3 <-coo_center (FeedSaxOut2)
FeedSaxOut.l <- filter(FeedSaxOut2, View == "l")#creates subset of lateral views only
FeedSaxOut.d <- filter(FeedSaxOut2, View == "d")#creates subset of lateral views only
FeedSaxOut.l.efour <- efourier(FeedSaxOut.l, nb.h=8, norm = FALSE, start = FALSE)
FeedSaxOut.d.efour <- efourier(FeedSaxOut.d, nb.h=8, norm = FALSE, start = FALSE)
FeedSaxOut.d.l <- combine (FeedSaxOut.d, FeedSaxOut.l)

#PCA
calibrate_harmonicpower_efourier(FTWRefOut.d.l,nb.h=12)#tells you how many harmonics will be needed to gather x% of harmonic power- here 8= 99%
FTWRefOut.d.l.efour <- efourier(FTWRefOut.d.l, nb.h=8, norm = FALSE, start = TRUE)# for why norm= false see note on ?efourier regarding roughly circular objects
FTWRefOut.d.l.pca <- PCA(FTWRefOut.d.l.efour, scale= FALSE, center= TRUE)
FTWRefOut.l.pca <- PCA(FTWRefOut.l.efour, scale= FALSE, center= TRUE)
plot(FTWRefOut.l.pca, cex= 1, zoom = 0.8, points= TRUE, 'taxon.code', labelspoints= FALSE)
FTWRefOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% PCA (scale=FALSE, center= TRUE) %>% plot ( cex= 0.8, zoom = 0.8, labelspoints=FALSE, points= TRUE,'taxon.code')
#and polar
FTWRefOut.d.l.p.pca <- PCA(FTWRefOut.d.l.efour, scale= FALSE, center= TRUE)
FTWRefOut.d.l.p %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% PCA (scale=FALSE, center= TRUE) %>% plot ( cex= 0.8, zoom = 0.8, labelspoints=FALSE, points= TRUE,'taxon.code')

#(Procrustes) MANOVA of Modern Reference material (between taxa vs between individuals vs random grouping)
export(FTWRefOut.l.efour)
export(FTWRefOut.d.efour)
FTWRef_l <- read.csv("D:\\Coordinate_Frames\\FTWRef_l.csv")
FTWRef_d <- read.csv("D:\\Coordinate_Frames\\FTWRef_d.csv")
#Calculate PCA scores for above
PCAObservation2_d<-prcomp(FTWRef_d)$x 
PCAObservation2_l<-prcomp(FTWRef_l)$x 
#Creates 4 level session factor i.e. taxon distinction (coefficients reordered as Taes, Tdur, Tturg, Tcom)
taxonfactor1<-as.factor(gl(4, 81))
taxonfactor1
#Create random factor- grains evenly split between 4 groups
randomfactor1<-as.factor(rep(1:4,81))
randomfactor1
#Create 324 level individual factor 
individualfactor<-as.factor(gl(324,1))
individualfactor
#Geomorph data frames for PCA computations of coefficients
gdf2_d<-geomorph.data.frame(shape=PCAObservation2_d) #ID: a vector containing the specimens ID 
gdf2_l<-geomorph.data.frame(shape=PCAObservation2_l) #ID: a vector containing the specimens ID 
#taxon factor
summary(procD.lm(shape~taxonfactor1, data=gdf2_d))
summary(procD.lm(shape~taxonfactor1, data=gdf2_l))
#individual factor
mod<-summary(procD.lm(shape~individualfactor1, data= gdf2_d))
mod2<-summary(procD.lm(shape~individualfactor1, data= gdf2_l))
mod
mod2
#random factor 
summary(procD.lm(shape~randomfactor1, data=gdf2_d))
summary(procD.lm(shape~randomfactor1, data=gdf2_l))

#(Procrustes) MANOVA of Archaeological material (between taxa vs between individuals)
export(FeedSaxOut.l.efour)
export(FeedSaxOut.d.efour)
FeedSaxRef_l <- read.csv("D:\\Coordinate_Frames\\FeedSax_l.csv")
FeedSaxRef_d <- read.csv("D:\\Coordinate_Frames\\FeedSax_d.csv")
#Calculate PCA scores for above
PCAObservation3_d<-prcomp(FeedSaxRef_d)$x 
PCAObservation3_l<-prcomp(FeedSaxRef_l)$x 
#Creates 4 level session factor (consecutive) i.e. taxon distinction (coefficients reordered as Taes, Tdur, Tturg, Tcom)
taxonfactor2<-as.factor(c(rep(1:1,33),rep(2:2,195), rep(3:3, 147), rep(4:4, 88)))
taxonfactor2
#Create random factor (repeated)
randomfactor2<-as.factor(c(rep(1:4,115),rep(1:3,1)))
randomfactor2
#Create 463 level individual factor (repeated)
indivfactor2<-as.factor(gl(463,1))
indivfactor2
#Geomorph data frames for PCA computations of coefficients
gdf3_d<-geomorph.data.frame(shape=PCAObservation3_d) #ID: a vector containing the specimens ID 
gdf3_l<-geomorph.data.frame(shape=PCAObservation3_l) #ID: a vector containing the specimens ID 
#taxon factor
summary(procD.lm(shape~taxonfactor2, data=gdf3_d))
summary(procD.lm(shape~taxonfactor2, data=gdf3_l))
#individual factor
mod<-summary(procD.lm(shape~indivfactor2, data= gdf3_d))
mod2<-summary(procD.lm(shape~indivfactor2, data= gdf3_l))
mod
mod2
#random factor
summary(procD.lm(shape~randomfactor2, data=gdf3_d))
summary(procD.lm(shape~randomfactor2, data=gdf3_l))

#(Procrustes) MANOVA of Archaeological material- by site 
FeedSaxRefS_l <- read.csv("D:\\Coordinate_Frames\\FeedSax_l_Site.csv")
FeedSaxRefS_d <- read.csv("D:\\Coordinate_Frames\\FeedSax_d_Site.csv")
#Calculate PCA scores for above
PCAObservation4_d<-prcomp(FeedSaxRefS_d)$x 
PCAObservation4_l<-prcomp(FeedSaxRefS_l)$x 
#Creates 4 level session factor (consecutive) i.e. taxon distinction (coefficients reordered as Taes, Tdur, Tturg, Tcom)
sitefactor<-as.factor(c(rep(1:1,34),rep(2:2,30),rep(3:3,26),rep(4:4,58),rep(5:5,31),rep(6:6,28),rep(7:7,46),rep(8:8,20),rep(9:9,40),rep(10:10,86),rep(11:11,27),rep(12:12,37)))
sitefactor
#Create random factor (repeated)
randomfactor3<-as.factor(c(rep(1:12,38),rep(1:7,1)))
randomfactor3
#Geomorph data frames for PCA computations of coefficients
gdf4_d<-geomorph.data.frame(shape=PCAObservation4_d) #ID: a vector containing the specimens ID 
gdf4_l<-geomorph.data.frame(shape=PCAObservation4_l) #ID: a vector containing the specimens ID 
#site factor
summary(procD.lm(shape~sitefactor, data=gdf3_d))
summary(procD.lm(shape~sitefactor, data=gdf3_l))
#random factor
summary(procD.lm(shape~randomfactor3, data=gdf3_d))
summary(procD.lm(shape~randomfactor3, data=gdf3_l))
