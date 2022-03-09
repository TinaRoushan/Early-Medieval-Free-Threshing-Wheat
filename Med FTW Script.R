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
#comment2222