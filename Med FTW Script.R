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
  ###Filter and Fourier transform
  FTWRefOut.l <- filter(FTWRefOut2, View == "l")#creates subset of lateral views only
  FTWRefOut.d <- filter(FTWRefOut2, View == "d")#creates subset of dorsal views only
  FTWRefOut.p <- filter(FTWRefOut2, View == "p")#creates subset of dorsal views only
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
  ###Filter and Fourier transform
  FTWRefOut.l <- filter(FTWRefOut2, View == "l")#creates subset of lateral views only
  FTWRefOut.d <- filter(FTWRefOut2, View == "d")#creates subset of dorsal views only
  FTWRefOut.p <- filter(FTWRefOut2, View == "p")#creates subset of dorsal views only
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
  FTWRefCoords
  FTWRefFrame <- read.csv("D:\\Coordinate_Frames\\Charred_Frame_FT2.csv", header = TRUE)
  FTWRefTxt <- import_txt(FTWRefCoords, fileEncoding="UTF-8-BOM")#fileencoding gets rid of weird symbols
  FTWRefOut <- Out(FTWRefTxt, fac=FTWRefFrame) #creates an Out object from a specified list- you can specify landmarks in "ldk"
  FTWRefOut2 <-coo_scale (FTWRefOut)
  ###Filter and Fourier transform
  FTWRefOut.l <- filter(FTWRefOut2, View == "l")#creates subset of lateral views only
  FTWRefOut.d <- filter(FTWRefOut2, View == "d")#creates subset of dorsal views only
  FTWRefOut.p <- filter(FTWRefOut2, View == "p")#creates subset of dorsal views only
  FTWRefOut.l.efour <- efourier(FTWRefOut.l, nb.h=8, norm = FALSE, start = TRUE)
  FTWRefOut.d.efour <- efourier(FTWRefOut.d, nb.h=8, norm = FALSE, start = TRUE)
  FTWRefOut.p.efour <- efourier(FTWRefOut.p, nb.h=8, norm = FALSE, start = TRUE)
  FTWRefOut.d.l <- combine (FTWRefOut.d, FTWRefOut.l)#dataset of just dorsal and lateral views
  FTWRefOut.d.l.p <- combine (FTWRefOut.d, FTWRefOut.l, FTWRefOut.p)#dataset of just dorsal and lateral views
  #LDA
  #LDA 
  FTWRefOut.d %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
  FTWRefOut.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
  FTWRefOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 

  ##Final dataset
  FTWRefCoords <- list.files("D:\\Coordinates\\Free Threshing\\FeedSax_Ref", full.names = TRUE)
  FTWRefCoords
  FTWRefFrame <- read.csv("D:\\Coordinate_Frames\\Charred_Frame_FeedSax_Ref.csv", header = TRUE)
  FTWRefTxt <- import_txt(FTWRefCoords, fileEncoding="UTF-8-BOM")#fileencoding gets rid of weird symbols
  FTWRefOut <- Out(FTWRefTxt, fac=FTWRefFrame) #creates an Out object from a specified list- you can specify landmarks in "ldk"
  FTWRefOut2 <-coo_scale (FTWRefOut)
  ###Filter and Fourier transform
  FTWRefOut.l <- filter(FTWRefOut2, View == "l")#creates subset of lateral views only
  FTWRefOut.d <- filter(FTWRefOut2, View == "d")#creates subset of dorsal views only
  FTWRefOut.p <- filter(FTWRefOut2, View == "p")#creates subset of dorsal views only
  FTWRefOut.l.efour <- efourier(FTWRefOut.l, nb.h=8, norm = FALSE, start = TRUE)
  FTWRefOut.d.efour <- efourier(FTWRefOut.d, nb.h=8, norm = FALSE, start = TRUE)
  FTWRefOut.p.efour <- efourier(FTWRefOut.p, nb.h=8, norm = FALSE, start = TRUE)
  FTWRefOut.d.l <- combine (FTWRefOut.d, FTWRefOut.l)#dataset of just dorsal and lateral views
  FTWRefOut.d.l.p <- combine (FTWRefOut.d, FTWRefOut.l, FTWRefOut.p)#dataset of just dorsal and lateral views
  #LDA
  #LDA 
  FTWRefOut.d %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
  FTWRefOut.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
  FTWRefOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
  FTWRefOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('Accession', scale=FALSE, center= TRUE) 
   FTWRefOut.d.l.lda <- FTWRefOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE, center= TRUE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
  FTWRefOut.d.l.lda
  plot(FTWRefOut.d.l.lda, points= TRUE , labelspoints= FALSE, zoom= 1.4, cex.labelsgroups= 0.8, rect.labelsgroups= TRUE, cex = 0.1)
  
  
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

  #Comparison with arch
 #Unknown LDA reclass
  UnclassCoords <- list.files("D:\\Coordinates\\Free Threshing\\Unclassified", full.names = TRUE)
  UnclassFrame <- read.csv("D:\\Coordinate_Frames\\Unclassified3.csv", header = TRUE)
  UnclassTxt <- import_txt(UnclassCoords, fileEncoding="UTF-8-BOM")#fileencoding gets rid of weird symbols
  UnclassOut <- Out(UnclassTxt, UnclassFrame) #creates an Out object from a specified list- you can specify landmarks in "ldk"
  UnclassOut %>% coo_center %>% coo_scale %>% stack # centres --> scales (normalises) ---> stacks your Out objects by centroid size (?check this is right)
  UnclassOut1 <-coo_scale(UnclassOut)
#reclass
  UnclassOut.d <- filter(UnclassOut1, View == "d")
  UnclassOut.l <- filter(UnclassOut1, View == "l")
  UnclassOut.d.l <- combine (UnclassOut.d, UnclassOut.l)
  UnclassOut.d.l.efour <- UnclassOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine
  reLDA(UnclassOut.d.l.efour, FTWRefOut.d.l.lda)
  #with germinated 
  reLDA(UnclassOut.d.l.efour, FTWRefOutG.d.l.lda)
  #with priors modified for even probability between groups 
  mod<-FTWRefOutG.d.l.lda$mod
  newdata<-UnclassOut.d.l.efour
  predicted<-predict(mod, newdata$coe, prior= c(1,1,1,1,1)/5)
  predicted
  
  ######## 
  
#visualise
 #dorsal
    FTWRefOut.d.site <-filter (FTWRefOut.d, Accession == "Germinated")
  #centre, scale by centroid size and stack outlines
  FTWRefOut.d.site %>% coo_center %>% coo_scale %>% stack
#lateral
   FTWRefOut.l.site <-filter (FTWRefOut.l, Accession == "High09")
  #centre, scale by centroid size and stack outlines
  FTWRefOut.l.site %>% coo_center %>% coo_scale %>% stack
#Mean shapes
 #MANOVA/Mean shapes
  MeanDorsal<-MSHAPES(FTWRefOut.l.efour, fac = 'taxon.code', FUN = mean, nb.pts = 120)
  MeanDorsal2<-MSHAPES(FTWRefOut.d.efour, fac = 'taxon.code', FUN = mean, nb.pts = 120)
  MeanDorsal.pca <- PCA(MeanDorsal, scale= FALSE, center= TRUE)
  FTW.ms <- MSHAPES(FTWRefOut.d.efour, 'taxon.code')
  FTW.ms2<- MSHAPES(FTWRefOut.l.efour, 'taxon.code')
  Out(FTW.ms$shp) %>% panel(names=TRUE) %>% plot_mshapes()
  Taes_com <- FTW.ms$shp$Taes_com    %>% coo_plot(border="blue")
  Timo <- FTW.ms$shp$Timo    %>% coo_plot(border="red")
  #pairwise comparison
  plot_mshapes (FTW.ms, size= 3/4)
  plot_mshapes (FTW.ms2, size= 3/4)
  
  #no durum
  Meandorsal2 <- filter(FTWRefOut.d.efour, Helper == "x")#
  Meanlateral2 <- filter(FTWRefOut.l.efour, Helper == "x")#
  FTW.ms3 <- MSHAPES(Meandorsal2, 'taxon.code2')
  FTW.ms4<- MSHAPES(Meanlateral2, 'taxon.code2')
  plot_mshapes (FTW.ms3, size= 3/4)
  plot_mshapes (FTW.ms4, size= 3/4)
  
  #pca
  MethodTest<- FTWRefOut.d.l%>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% PCA (scale=FALSE, center= TRUE) 
  MethodTestPlot<-plot(MethodTest, cex= 1, zoom = 1, points= TRUE, labelspoints= FALSE,'taxon.code')  
  
  
  