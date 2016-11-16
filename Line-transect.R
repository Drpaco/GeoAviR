###Complete example for line-transect estimates from GeoAviR
#FBolduc
#18 dec 2014

require(GeoAviR)

data(quebec)
d<-filterECSAS(quebec)


### spatial object from watches lat-long
require(rgdal)
transect <- data.frame(lat=d$LatStart,lon=d$LongStart)
coordinates(transect) <- ~lon + lat
transect<-SpatialPointsDataFrame(transect,data=d[,"Count",drop=FALSE])
proj4string(transect)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")


dist_res <- "U:\\Oiseaux souillés\\Inventaires pélagiques\\GeoAviR\\InterfaceGraphique0.1\\Output\\Exemples"
setwd(dist_res)



#### earth background
setwd("U:\\GIS\\Terreshp")
shpm <- readOGR(".","ne_10m_land")


###ensure projection agreement
prj <- proj4string(transect)
shpm<-spTransform(shpm,CRS(prj)) 


### Build a 50000 x 50000 meters grid - 
size<-100000
new.grid<-create.grid(Latitude=c(44,52),Longitude=c(-70,-56),Grid.size=c(size,size),Clip=F,clip.shape=shpm,projection=CRS(prj))
new.grid$ID<-paste("parc",new.grid$ID,sep="")
#proj4string(new.grid)


#keep only visited cells
new.grid<-new.grid[apply(gIntersects(transect,new.grid,byid=TRUE),1,any),] 

###as an option?
##clip new.grid w/earth, as an option, better for coastal cells
# test <- gDifference(new.grid, shpm, byid = T)# note: may take several minutes
# plot(test)
test <- new.grid

##compute cell/zone area (km2) - may need to clip cells with earth for exact area along the coast
#needs lambert projection for distance
prjm <- "+proj=lcc +lat_1=46 +lat_2=60 +lat_0=44 +lon_0=-68.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
test<-spTransform(test,CRS(prjm))
area <- data.frame(km2=gArea(test,byid=T)/1000000)
area$ID <- 1:nrow(area)
area$ID <- paste("parc",area$ID,sep="")
test <- SpatialPolygonsDataFrame(test,data=area)
new.grid<-spTransform(test,CRS(prj)) 
#plot(test[test$km2<500,],add=T,col="red")
# plot(new.grid)
# plot(transect,add=T)

### Overlay transects and grid and attribute squares to observations
x<-over(transect,new.grid)
d$square<-x$ID
d$square_area <- x$km2
d<-d[!is.na(d$square),]


##sampling unit = choose var-
#names(d)
sample <- names(d)[6] #"Date" #Select variable WatchID, or Date or else
sample
#sample labels
d$SMP_LABEL<-paste(d$square,d[,c(sample)],sep="_")

### Construct sample labels considering that day transects can overlap with multiple squares
temp<-aggregate(WatchLenKm~SMP_LABEL,data=unique(d[,c("SMP_LABEL","WatchID","WatchLenKm")]),sum)
names(temp)[2]<-"SMP_EFFORT"
d<-merge(d,temp,sort=FALSE)
#stats to merge with shp later
# kms <- aggregate(WatchLenKm~square,data=unique(d[,c("square","WatchID","WatchLenKm")]),sum)#effort total par zone
# samples <- aggregate(.~square,data=unique(d[,c("square",sample)]),length) #nb échantillons par zone
# names(samples)[2]<-"samples"
# seastate <- aggregate(seastate~SMP_LABEL,data=unique(d[,c("SMP_LABEL","WatchID","seastate")]),mean)
# birds <- aggregate(Count~square,data=unique(d[,c("square","WatchID","Count")]),sum)
# IDlevels <- aggregate(.~square,data=unique(d[,c("square","Alpha")]),length) #includes unknowns for now
# cruises <- aggregate(.~square,data=unique(d[,c("square","CruiseID")]),length)
# watchid <- aggregate(.~square,data=unique(d[,c("square","WatchID")]),length)
# dates <- aggregate(.~square,data=unique(d[,c("square","Date")]),length)

#d<-d[,c("square","square_area","Date","SMP_LABEL","SMP_EFFORT","Distance","Count","Alpha")]
dd<-ddply(d,.(SMP_LABEL),function(i){sum(i$Count,na.rm=TRUE)}) #eliminate duplicate lines for transect without observations
dd<-dd[dd$V1==0,] #get the label name for transect without observations
d<-d[(d$SMP_LABEL%in%dd$SMP_LABEL & !duplicated(d$SMP_LABEL)) | (!d$SMP_LABEL%in%dd$SMP_LABEL & !(d$Alpha=="")),] #keep only lines for empty transects or non-empty lines for non-empty transects
d<-d[order(d$square),]

###distance sampling
path<-"C:/temp/distance"
pathMCDS<-"U:\\Oiseaux souillés\\Inventaires pélagiques\\GeoAviR"
breaks<-c(0,50,100,200,300)#un peu redondant avec filtre/données
SMP_LABEL<-"SMP_LABEL"
SMP_EFFORT<-"SMP_EFFORT"
DISTANCE<-"Distance"
SIZE<-"Count"
STR_LABEL<-"square"
STR_AREA<-"square_area"
lsub<-list(Alpha=c("HERG","NOFU","BLKI")) 
split<-TRUE
stratum<-"STR_LABEL"
detection<-"All"
empty<-NULL

#estimators
x<-distance.wrap(d,stratum=stratum,empty=empty,detection=detection,lsub=lsub,split=split,
                 path=path,pathMCDS=pathMCDS,breaks=breaks,STR_LABEL=STR_LABEL,STR_AREA=STR_AREA,
                 SMP_LABEL=SMP_LABEL,SMP_EFFORT=SMP_EFFORT,DISTANCE=DISTANCE,SIZE=SIZE,
                 units = list(Type = "Line", Distance = "Perp", Length_units = "Kilometers", 
                              Distance_units ="Meters", Area_units = "Square kilometers"), verbose=FALSE)


#choose best model for one species of the lsub object
all.sp.best <- keep.best.model(x$NOFU)
#if lsub=NULL, all.sp.best <- keep.best.model(x)
mod.selected <- which.min(sapply(1:6, function(i)x$NOFU[[i]]$AIC[3])) 
global.summary(model=all.sp.best, species="NOFU", file="alcidae_global", directory="C:/temp/distance")


###################################
######maps######

####Extract the probability of detection


###extract predictions
tmp <- all.sp.best$ density_estimate$Stratum

densities <- tmp[tmp$Parameters == "D",c("Stratum","Estimates","% of var.")]
####save shp - grille + données assocées aux cellules
#    densities$Estimates <- as.numeric(densities$Estimates)
names(densities)[3] <- "CoefVar" #probleme avec % dans nom
#    densities$CoefVar <- as.numeric(densities$CoefVar)
names(densities)[names(densities) == "Stratum"] <- "ID"

#     ###merge stats with densities
#     densities <- merge(densities,samples,by.x="ID",by.y="square")
#     densities <- merge(densities,kms,by.x="ID",by.y="square")
#     densities <- merge(densities,dates,by.x="ID",by.y="square")
#     densities <- merge(densities,cruises,by.x="ID",by.y="square")
#     densities <- merge(densities,watchid,by.x="ID",by.y="square")
#     densities <- merge(densities,IDlevels,by.x="ID",by.y="square")
#     densities <- merge(densities,birds,by.x="ID",by.y="square")
#     densities$CoefVar <- ifelse(is.na(densities$CoefVar),0,densities$CoefVar)#avoid NAs


#join les estimations aux shp des cellules
#changement dans l'ordre seulement, mais a des impacts dans la suite des choses
require(maptools)
o <- match(new.grid$ID, densities$ID) #IDs for new.grid before??
temp <- densities[o,]
row.names(temp) <- row.names(new.grid)
new.grid2 <- spCbind(new.grid,temp)
new.grid2 <- new.grid2[!is.na(new.grid2$Estimates),] #enlever les cellules non visitées dans le subset
#test <- new.grid2@data


##change extremes with 99% percentile
#put as a checkbox as an option
#           p99 <- quantile(new.grid2@data$Estimates, c(.995)) #ok
#         new.grid2@data$Estimates <- ifelse(new.grid2@data$Estimates > p99,p99,new.grid2@data$Estimates)


#select variable to map
#names(new.grid2@data)
#changement dans le nom de l'objet pour avoir plus de flexibilité à l'avenir dans la variable à cartographier
map <- names(new.grid2@data)[4] #or else
map

##change extremes with 99% percentile
# p99 <- quantile(densities$Estimates, c(.995)) #ok
# densities$Estimates <- ifelse(densities$Estimates > p99,p99,densities$Estimates)

brks <- quantile(new.grid2@data[,c(map)], c(0,.5,.75,.95))
brks[length(brks)+1] <- max(new.grid2@data[,c(map)])+0.01
brks <- round(brks,2)
#    brks
#select classification scheme
#     classify <- "even" #or "kmeans", or "quantile", qui envoi à l'Ancienne boucle qui créer les quantiles
#     nbclasses <- 10 #comme il y a déjà une entrée dansgeoavirweb

# brks <- quantile(new.grid2@data[,c(map)], c(0,.5,.75,.95))
# brks
# ?quantile


#if classify = even
#         brks <- round(seq(0, max(new.grid2@data[,c(map)],na.rm = T)+0.01, length.out = nbclasses), 2)
#         brks

#if classify = kmeans
#         library(classInt) 
#         classes_km <- classIntervals(new.grid2@data[,c(map)], n=nbclasses, style = "kmeans", rtimes = 1)
#         classes_km$brks[1] <- ifelse(classes_km$brks[1]!=0,0,classes_km$brks[1])
#         brks <- round(classes_km$brks,2)
#     #brks
#     if(brks[length(brks)] < max(new.grid2@data[,c(map)])) {brks[length(brks)] <- brks[length(brks)] +0.01}                   

##associate data w/brks intervals
tags<-vector()
for (i in 1:length(brks)) {
  if (i ==1 ) {
    new.grid2@data$class  <- ifelse(new.grid2@data[,c(map)] == brks[i], as.character(brks[i]),"-")
    new.grid2@data$classno  <- ifelse(new.grid2@data[,c(map)] == brks[i], as.numeric(i),"-")
    tags[1]<- as.character(brks[i])
    
  }
  if(i>1) {
    new.grid2@data$class  <- ifelse(new.grid2@data[,c(map)] > brks[i-1] & new.grid2@data[,c(map)] <= brks[i], paste(" > ",brks[i-1]," - ",brks[i],sep=""), new.grid2@data$class)
    new.grid2@data$classno  <- ifelse(new.grid2@data[,c(map)] > brks[i-1] & new.grid2@data[,c(map)] <= brks[i], as.numeric(i), new.grid2@data$classno)
    tags[i]<- paste(" > ",brks[i-1]," - ",brks[i],sep="")
  }
}


#colors for legend
library(RColorBrewer) 
br.palette <- colorRampPalette(c("green","yellow", "red"), space = "rgb")
nb<- length(brks)-1 
br.palette(nb)
pal <- br.palette(n=nb)
pal <- c("lightgray",pal)#ajout class 0
new.grid2$color <- pal[as.numeric(new.grid2$classno)]


#titre <- "Corrected densities"
SousTitre <- "BLKI" #nom du sous groupe sélectionné 
setwd("C:\\temp\\distance")
legendtitle <- "birds/km2"
new.grid2$classno <- as.numeric(new.grid2$classno)
new.grid2 <- new.grid2[order(new.grid2$classno),]

#classes <- unique(new.grid2$class) 
classes <- unique(new.grid2$class[!is.na(new.grid2$class)])
classes

#l'emplacement de la legende pourrait être problématique...
#png(paste(titre,".png",sep=""),width=960,height=960)
plot(new.grid2,bg="lightblue",border="lightblue",axes=T)
for (i in 1:length(classes)){
  plot(new.grid2[new.grid2$class[!is.na(new.grid2$class)]==classes[i],],col=new.grid2@data[new.grid2$class[!is.na(new.grid2$class)]==classes[i],"color"],border="darkgray",add=T)
}
plot(shpm,add=T,col="darkkhaki",border="darkkhaki",axes=T)
#title(main =titre, cex.main=1.5)
mtext(SousTitre,cex=1)
legend("bottomright", bty = "n",  legend = tags,fill=pal, title=legendtitle,cex=0.6)
#a chaque fois que l'utilisateur visualize une stat, sa classification est sauvée en changeant son nom:
#     names(new.grid2)[names(new.grid2) == "class"] <- paste("class",map,sep="_")
#     names(new.grid2)[names(new.grid2) == "classno"] <- paste("classno",map,sep="_")
#dev.off()
#test <- new.grid2@data

#données sauvées dans un fichier txt - pour join avec le shp dans une autre application, par exemple
filename <- "densities"
storecsv <- dist_res
name <- paste(storecsv, "//",filename, ".txt", sep="")
write.csv(new.grid2@data,name)

##save corrected estimates 
nom <- "testshp"
writeOGR(new.grid2, ".", nom, driver="ESRI Shapefile",overwrite_layer=T)
