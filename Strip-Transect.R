###exemple complet avec carte et fichier quebec
#FBolduc
#18 dec 2014

require(GeoAviR)
#test

data(quebec)
d<-filterECSAS(quebec)


### faire un shapefile à partir des début d'observations dans les données
require(rgdal)
transect <- data.frame(lat=d$LatStart,lon=d$LongStart)
coordinates(transect) <- ~lon + lat
transect<-SpatialPointsDataFrame(transect,data=d[,"Count",drop=FALSE])
proj4string(transect)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")


dist_res <- "U:\\Oiseaux souillés\\Inventaires pélagiques\\GeoAviR\\InterfaceGraphique0.1\\Output\\Exemples"
setwd(dist_res)

setwd(paste(getwd(), "/Canada",sep=""))
setwd("..") #remonte d'un cran ds l'arborescence  (si nécéssaire)

#récupération des données dans un shp local - 
can1 <- readOGR(dsn = paste(getwd(),"/Canada",sep = ""),layer = "CAN_adm1")
can1@data[,"NAME_1"] <- as.character(can1@data[,"NAME_1"])    
can1@data[11,"NAME_1"] <- "Quebec"
can1@data <- can1@data[,c(1,2,3,4,5)]
colnames(can1@data) <- c("ID_0","Pays","ID_1","Province","Type")

posm <- which(can1@data[,"Type"] == "Quebec" |
                can1@data[,"Type"] == "Newfoundland and Labrador" |
                can1@data[,"Type"] == "Prince Edward Island" |
                can1@data[,"Type"] == "Nova Scotia"  |
                can1@data[,"Type"] == "New Brunswick" ) 
shpm <- can1[posm,]

###mettre les données et le fond de carte ensemble - assurer la même projection au cas où
prj <- proj4string(transect)
shpm<-spTransform(shpm,CRS(prj)) 


### Build a 50000 x 50000 meters grid - lat et long à partir du jeu de données? on voir plus tard avec plot(transect) que la grille sert aussi à restreinfre l'étendue spatiale
size<-50000
new.grid<-create.grid(Latitude=c(44,52),Longitude=c(-70,-56),Grid.size=c(size,size),Clip=FALSE,clip.shape=canshp,projection=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
new.grid$ID<-paste("parc",new.grid$ID,sep="")

#couper la grille pour ne garder que les cellules de la grille visitées
new.grid<-new.grid[apply(gIntersects(transect,new.grid,byid=TRUE),1,any),] #limitation de taille de vecteur?

### Overlay transects and grid and attribute squares to observations
x<-over(transect,new.grid)
d$square<-x$ID
d<-d[!is.na(d$square),]
d$square_area<-(size/1000)^2 #in kilometers
d$SMP_LABEL<-paste(d$square,d$Date,sep="_")

### Construct sample labels considering that day transects can overlap with multiple squares
#demander à l'utilisateur, quelle variable sera le niveau de l'échantillon (WatchID, Date, autre?)
temp<-aggregate(WatchLenKm~SMP_LABEL,data=unique(d[,c("SMP_LABEL","WatchID","WatchLenKm")]),sum)
names(temp)[2]<-"SMP_EFFORT"
d<-merge(d,temp,sort=FALSE)
d<-d[,c("square","square_area","Date","SMP_LABEL","SMP_EFFORT","Distance","Count","Alpha")]
dd<-ddply(d,.(SMP_LABEL),function(i){sum(i$Count,na.rm=TRUE)}) #eliminate duplicate lines for transect without observations
dd<-dd[dd$V1==0,] #get the label name for transect without observations
d<-d[(d$SMP_LABEL%in%dd$SMP_LABEL & !duplicated(d$SMP_LABEL)) | (!d$SMP_LABEL%in%dd$SMP_LABEL & !(d$Alpha=="")),] #keep only lines for empty transects or non-empty lines for non-empty transects
d<-d[order(d$square),]

###distance sampling
path<-"C:/temp/distance"
pathMCDS<-"U:\\Oiseaux souillés\\Inventaires pélagiques\\GeoAviR"
SMP_LABEL<-"SMP_LABEL"
SMP_EFFORT<-"SMP_EFFORT"
SIZE<-"Count"
STR_LABEL<-"square"
STR_AREA<-"square_area"
lsub<-list(Alpha=c("NOFU","NOGA")) 
split<-TRUE
stratum<-"STR_LABEL"
empty<-NULL
detection <- "All"
breaks<-c(0,300)

#detection<-"All"
#DISTANCE <- NULL

#DISTANCE<-"Distance"
d$Distance <- NULL #remove distance var from dataset


##strip estimates

x<-strip.wrap(d,stratum=stratum,empty=empty,detection=detection,lsub=lsub,split=split,
              path=path,pathMCDS=pathMCDS,breaks=breaks,STR_LABEL=STR_LABEL,STR_AREA=STR_AREA,
              SMP_LABEL=SMP_LABEL,SMP_EFFORT=SMP_EFFORT,SIZE=SIZE,
              units = list(Type = "Line", Distance = "Perp", Length_units = "Kilometers", 
                           Distance_units ="Meters", Area_units = "Square kilometers"), verbose=FALSE)


#choose best model for one species of the lsub object
#all.sp.best <- keep.best.model(x$NOFU)
#if lsub=NULL, all.sp.best <- keep.best.model(x)
# mod.selected <- which.min(sapply(1:6, function(i)x$NOFU[[i]]$AIC[3])) 
# global.summary(model=x, species="NOFU", file="alcidae_global", directory="C:/temp/distance")





#x<-strip.wrap(d,stratum=stratum,empty=empty,lsub=lsub,split=split,path=path,pathMCDS=pathMCDS,STR_LABEL=STR_LABEL,STR_AREA=STR_AREA,SMP_LABEL=SMP_LABEL,SMP_EFFORT=SMP_EFFORT,SIZE=SIZE,verbose=FALSE)

global.summary(model=x[["NOFU"]], species=c("fulmars"), file="DetModel-NOFU2", directory="C:/temp/distance")



###################################
######fabrication de la carte######


#aller chercher les résultats dans l'objet x
#une espèce en particulier :
tmp <- x$NOFU$density_estimate$Stratum 
######si aucun sous-groupe = x$density_estimate$Stratum#####

densities <- tmp[tmp$Parameters == "D",c("Stratum","Estimates","% of var.")]

####save shp - grille + données assocées aux cellules
densities$Estimates <- as.numeric(densities$Estimates)
names(densities)[3] <- "CoefVar" #probleme avec % dans nom
densities$CoefVar <- as.numeric(densities$CoefVar)
names(densities)[names(densities) == "Stratum"] <- "ID"

##change extremes with 99% percentile
p99 <- quantile(densities$Estimates, c(.995)) #ok
densities$Estimates <- ifelse(densities$Estimates > p99,p99,densities$Estimates)

brks <- quantile(densities$Estimates, c(0,.5,.75,.95))
brks[length(brks)+1] <- max(densities$Estimates)+0.001
brks <- round(brks,3)


##associate data w/quantile intervals
for (i in 1:length(brks)) {
  if (i ==1 ) {
    densities$class  <- ifelse(densities$Estimates == brks[i], as.character(brks[i]),"-")
  }
  if(i>1) {
    densities$class  <- ifelse(densities$Estimates > brks[i-1] & densities$Estimates <= brks[i], paste(i-1,": > ",brks[i-1]," - ",brks[i],sep=""), densities$class)
  }
}

#données sauvées dans un fichier txt - pour join avec le shp dans une autre application, par exemple
# filename <- "densities"
# storecsv <- dist_res
# name <- paste(storecsv, "//",filename, ".txt", sep="")
# write.csv(densities,name)

#join les estimations aux shp des cellules
require(maptools)
o <- match(new.grid$ID, densities$ID)
temp <- densities[o,]
row.names(temp) <- row.names(new.grid)
new.grid2 <- spCbind(new.grid,temp)
new.grid2 <- new.grid2[!is.na(new.grid2$Estimates),] #enlever les cellules non visitées dans le subset


#colors for legend
library(RColorBrewer) 
br.palette <- colorRampPalette(c("green","yellow", "red"), space = "rgb")
nb<- length(brks)-1 #must be below 9 for now - class labeling order
br.palette(nb)
pal <- br.palette(n=nb)
pal <- c("lightgray",pal)#ajout class 0


##mapping in R
new.grid2 <- new.grid2[order(new.grid2$class),]
tags <- unique(new.grid2$class)

titre <- "Observed densities"
setwd("C:\\temp\\distance")
legendtitle <- "birds/km2"

#l'emplacement de la legende pourrait être problématique...
#png(paste(titre,".png",sep=""),width=960,height=960)
plot(new.grid2,axes=T)
plot(new.grid2[new.grid2$class==tags[1],],col=pal[1],add=T)
for (i in 2:length(tags)){
  plot(new.grid2[new.grid2$class==tags[i],],col=pal[i],add=T)
}
plot(shpm,add=T,col="darkkhaki")
title(main =titre, cex.main=1.5)
legend("bottomright", bty = "n",  legend = tags,fill=pal, title=legendtitle)
#dev.off()


##save corrected estimates 
nom <- "testshp"
writeOGR(new.grid2, ".", nom, driver="ESRI Shapefile",overwrite_layer=T)

