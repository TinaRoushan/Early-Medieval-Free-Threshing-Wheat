library('Momocs')
#Adding modern reference outlines 
FTWRefCoords <- list.files("D:\\Coordinates\\Free Threshing\\FeedSax_Ref", full.names = TRUE)
FTWRefFrame <- read.csv("D:\\Coordinate_Frames\\Charred_Frame_FeedSax_Ref.csv", header = TRUE)
FTWRefTxt <- import_txt(FTWRefCoords, fileEncoding="UTF-8-BOM")#fileencoding gets rid of weird symbols
FTWRefOut <- Out(FTWRefTxt, fac=FTWRefFrame) #creates an Out object from a specified list- you can specify landmarks in "ldk"
FTWRefOut2 <-coo_scale (FTWRefOut)
###Filter and Fourier transform
FTWRefOut.l <- filter(FTWRefOut2, View == "l")#creates subset of lateral views only
FTWRefOut.d <- filter(FTWRefOut2, View == "d")#creates subset of dorsal views only
calibrate_harmonicpower_efourier(FTWRefOut.l,nb.h=12)#tells you how many harmonics will be needed to gather x% of harmonic power- here 8= 99%
calibrate_harmonicpower_efourier(FTWRefOut.d,nb.h=12)#tells you how many harmonics will be needed to gather x% of harmonic power- here 8= 99%
FTWRefOut.l.efour <- efourier(FTWRefOut.l, nb.h=8, norm = FALSE, start = TRUE)
FTWRefOut.d.efour <- efourier(FTWRefOut.d, nb.h=8, norm = FALSE, start = TRUE)
FTWRefOut.d.l <- combine (FTWRefOut.d, FTWRefOut.l)#dataset of just dorsal and lateral views
#LDA 
FTWRefOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE) 
FTWRefOut.d.l.lda <- FTWRefOut.d.l %>% chop(~View) %>% lapply(efourier,nb.h=8, norm = FALSE, start = FALSE, center= TRUE) %>% combine %>% LDA ('taxon.code', scale=FALSE, center= TRUE)
plot(FTWRefOut.d.l.lda, points= TRUE , labelspoints= TRUE, zoom= 1.25, cex.labelsgroups= 0.95, rect.labelsgroups= TRUE, cex = 0.1)
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
