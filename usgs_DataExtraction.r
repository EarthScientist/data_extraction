library(raster)
library(sp)
library(maptools)

setwd("/workspace/UA/malindgren/projects/DataExtract_USGS/working")

# lets set some variables needed:
# this is the input folder that has the data to be listed
in_data <- "/Data/Base_Data/Climate/AK_CAN_2km/AK_CAN_2km/cru_TS31/historical/pr/"
out_data <- "/workspace/UA/malindgren/projects/DataExtract_USGS/working/"

# this is the input shapefile
pointSPDF <- readShapePoints("/workspace/UA/malindgren/projects/DataExtract_USGS/Data_given/orig_all_pts_shape.shp")

# here are a couple of things I need from that shapefile
# the xy coordinates
xy <- cbind(pointSPDF$POINT_X,pointSPDF$POINT_Y)

# the data.frame
site_names <- pointSPDF$Site_Name
# lets bring them together
xy_site_names <- cbind(xy, as.character(site_names))
#set some new column names
colnames(xy_site_names) <- c("POINT_X","POINT_Y","SITES")

# stuff which will be iterators through the years/months
years <- 1901:2009
months <- c("01","02","03","04","05","06","07","08","09","10","11","12")

# naming stuff to get everything in the list chronologically
variable <- "pr"
metric <- "total_mm"
model <- "cru_TS31"


# this is where we list the filenames we want to put into the raster stack
# a way to get the files in chronological order
list.2km <- character()
for(y in years){
	for(m in months){
		# summation of the naming vars into a single convention
		name_convention <- paste(in_data,variable,"_",metric,"_",model,"_",m,"_",y,".tif", sep="")
		list.2km <- append(list.2km, name_convention)
	}
}

# here we take the list we just created and turn it into a raster stack
s <- stack(list.2km)

# get the extent-class of an object of the stack
e <- extent(raster(s, layer=1))

# create a new blank raster with what we know about the layers in the stack to use in adjacent query
blankRaster <- raster(e, nrows=nrow(raster(s, layer=1)), ncols=ncol(raster(s, layer=1)), crs=projection(s))

# get the cell numbers from the xy coordinates of the input points
in_cells <- data.frame(cellFromXY(blankRaster, xy))

# now we need to unlist that var into something we can iterate through
in_cells <- unlist(in_cells)

# a counter for the below loop
count=0

for(i in in_cells){
	
	count=count+1

	# ask it which 8 cells are adjacent to the input cell
	adj_cells <- adjacent(blankRaster, i, directions=8, pairs=F, target=NULL, sorted=T, include=F, id=T)
	print(adj_cells)
	site_name <- xy_site_names[count,3]
	print(paste(count, site_name))
	
	if(count == 1){
		adj_cells_sites <- data.frame(adj_cells, site_name)	
	}else{
		adj_cells_sites_tmp <- data.frame(adj_cells, site_name)
		adj_cells_sites <- rbind(adj_cells_sites, adj_cells_sites_tmp)
	}
}

# create a new table of the xy values and a zero column since all of these are the adjacent cells not the focal 
adj_cells_extractor <- data.frame(adj_cells_sites, as.data.frame("0"))

# give some names to those new columns
colnames(adj_cells_extractor) <- c("CellNum","SiteName","FocalCell")

# here we are doing the same as above but with the focal cells
in_cells_extractor <- data.frame(cellFromXY(blankRaster, xy), xy_site_names[,3], as.data.frame("1"))

# give the same column names as above with the adjacent cells
colnames(in_cells_extractor) <- c("CellNum","SiteName","FocalCell")

# here is the combination of the input cells and the adjacent ones
combined_extractor <- rbind(in_cells_extractor, adj_cells_extractor)

# here we grab the xy values of the centroids of the input and adjacent pixels
out_xy <- xyFromCell(blankRaster, combined_extractor[,1], spatial=FALSE)

# this is the new extractor that will be used to get the desired data
combined_extractor_xy <- data.frame(combined_extractor[,1], out_xy, combined_extractor[,2:3])

# set their column names
colnames(combined_extractor_xy) <- c("CellNum","POINT_X","POINT_Y","SiteName","FocalCell")

# extract the values using a vector of the cell numbers
extraction <- extract(s, combined_extractor_xy[,1])

# now we want to create a new output table
output_extraction <- cbind(combined_extractor_xy[,2:ncol(combined_extractor_xy)], extraction)

# here I sort the entire data.frame based on the site_name variable
output_extraction <- output_extraction[order(output_extraction[,3]),]

# write out the final csv file with the extraction
write.csv(output_extraction, file=paste(out_data, variable,"_",metric,"_",model,"_","SNAP_extraction.csv", sep=""), row.names=F)

rm(list=ls())
# END