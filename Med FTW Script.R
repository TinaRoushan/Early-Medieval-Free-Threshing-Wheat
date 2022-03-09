packageurl <- "https://cran.r-project.org/src/contrib/Archive/Momocs//Momocs_1.2.9.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
packageVersion('Momocs')
library('Momocs')
#Adding modern reference outlines 
FTWRefCoords <- list.files("D:\\Coordinates\\Free Threshing\\FeedSax_Ref", full.names = TRUE)
FTWRefFrame <- read.csv("D:\\Coordinate_Frames\\Charred_Frame_FeedSax_Ref.csv", header = TRUE)
FTWRefTxt <- import_txt(FTWRefCoords, fileEncoding="UTF-8-BOM")#fileencoding gets rid of weird symbols
FTWRefOut <- Out(FTWRefTxt, fac=FTWRefFrame) #creates an Out object from a specified list- you can specify landmarks in "ldk"
FTWRefOut2 <-coo_scale (FTWRefOut)
###
FTWRefOut.l <- filter(FTWRefOut2, View == "l")#creates subset of lateral views only
FTWRefOut.d <- filter(FTWRefOut2, View == "d")#creates subset of dorsal views only
calibrate_harmonicpower_efourier(FTWRefOut.l,nb.h=12)#tells you how many harmonics will be needed to gather x% of harmonic power- here 8= 99%
calibrate_harmonicpower_efourier(FTWRefOut.d,nb.h=12)#tells you how many harmonics will be needed to gather x% of harmonic power- here 8= 99%
FTWRefOut.l.efour <- efourier(FTWRefOut.l, nb.h=8, norm = FALSE, start = TRUE)
