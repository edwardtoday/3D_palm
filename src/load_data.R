source("morphology.R")

dilate <- function (mat, strel)
{
	width <- dim(mat)[1]
	height <- dim(mat)[2]
	sw <- max(strel$xrange)
	sh <- max(strel$yrange)
	result <- matrix(0,nrow=width+2*sw,ncol=height+2*sh)
	for (k in strel$xrange)
		for(l in strel$yrange)
			if(!is.na(strel$value(k,l)) && strel$value(k,l))
			{
				for(w in k-sw:k+sw)
					for(h in l-sh:l+sh)
						result[w,h] <- sup(result[w,h],strel$value(k,l))
			}
	result <- result[sw:(sw+width),sh:(sh+height)]
	result
}

# datfile <- file("/Users/edwardtoday/polyu_3d_palmprint_classification/ROI/Sub3D_float/Sub3D_I_13_1.dat","rb")
# bindata <- readBin(datfile, double(), n = 128*128, size = 4, endian='little')
datfile <- file("/Users/edwardtoday/polyu_3d_palmprint_classification/3DPalm_zonly/3D_I_1_0.zonly","rb")
bindata <- readBin(datfile, double(), n = 768*576, size = 4, endian='little')
close(datfile)
matdata <- matrix(bindata,nrow=768)

roi_size <- 400
# Extract 400*400 Region of Interest
roi <- matdata[235:(234+roi_size),69:(68+roi_size)]

refplane <- roi[1:(roi_size*0.2),(roi_size*0.3):(roi_size*0.7)]
maxregion <- roi[(roi_size*0.4):(roi_size*0.95),(roi_size*0.3):(roi_size*0.7)]

# Calculte max depth (MD)
d_r <- max(refplane)
d_max <- min(maxregion)
MD <- d_r - d_max

# Show preview

# library(rgl)
# persp3d(roi,col='lightyellow')
# image(roi)
# image(refplane)
# image(maxregion)
# filled.contour(roi)

num_levels <- 8
step <- MD/num_levels

level <- matrix(NA, nrow=roi_size, ncol=roi_size)
# ptime <- system.time({
for(row in 1:roi_size)
{
	for(col in 1:roi_size)
	{
		d_pixel <- roi[row,col]
		if(d_pixel > d_r | d_pixel < d_max)
		{
			level[row,col] <- -1
		}		
		else if(d_pixel == d_max)
		{
			if( row > roi_size*0.4 & row < roi_size*0.95 &
				col > roi_size*0.3 & col < roi_size*0.7)
			{
				level[row,col] <- 0
			}
			else
			{
				level[row,col] <- -1
			}
		}
		else
		{
			pixel_layer <- 1 + as.integer((d_pixel-d_max)/step)
			level[row,col] <- pixel_layer
		}
	}
}
# })[3]
# ptime

# library(rgl)
# persp3d(level,col='lightyellow')
# image(level)

# disc <- structuring.element.disc(35-3*1)
# L_1 <- level==1
# DL_1 <- dilatation(L_1,disc)

HCA <- rep(0, num_levels)
for(l in 1:num_levels)
{
	HCA[l] = sum(level>=l)
}
HCA

# RLL
# From deepest point, go through the pixels on the rasterized line, stop when reaching -1, then calc the distance.
