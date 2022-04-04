  library('Momocs')
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
  FTWRefOut.d.l.lda <- FTWRefOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE, center= TRUE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
  plot(FTWRefOut.d.l.lda, points= TRUE , labelspoints= FALSE, zoom= 1.4, cex.labelsgroups= 0.8, rect.labelsgroups= TRUE, cex = 0.1)
  
   ###
  
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
 
  
  
  
#visualise
 #dorsal
    FTWRefOut.d.site <-filter (FTWRefOut.d, taxon.code == "Adur")
  #centre, scale by centroid size and stack outlines
  FTWRefOut.d.site %>% coo_center %>% coo_scale %>% stack
#lateral
   FTWRefOut.l.site <-filter (FTWRefOut.l, taxon.code == "Adur")
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
  