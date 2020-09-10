rm(list = ls()) # clear workspace/environment (all variables and functions)
cat("\014") # clear the console

library(stats)




validation <- 'zurich' #to choose from zurich, vumc, pmh, design, bd2decide

#radiomics data
rad <- read.csv("hpv_validation_zurich.csv",sep=",")


#clinical data
clinical <- read.csv("hpv_validation_zurich_clinical.csv",sep=",")

#selected features
namesfile <- read.csv("selected_features_hpv_zurich.csv",sep=",",header = FALSE, na = '')


# remove patients with less than 64 voxels 
ind_vol <- which(rad$voxels>64) 
rad <- rad[ind_vol,]
clinical <- clinical[ind_vol,]


names <- na.omit(t(namesfile[2,]))
names<- t(names)

fit_data <- subset(rad, select = names)
fit_data["c"]<-1
fit_data["HPV"]<-clinical$ï..HPV
fit_data["center"]<-clinical$center
fit_data<- fit_data[complete.cases(fit_data),]
#scaling
fit_data[,1:(dim(fit_data)[2]-3)]<- scale(fit_data[,1:(dim(fit_data)[2]-3)])


#split to centers
centers <- fit_data$center
centers_levels <- levels(centers)
#centers_levels <- centers_levels[centers_levels != 'design']

#fit_data <- fit_data[which(fit_data$center != 'design'),]

# remove a column
#fit_data = fit_data[,-c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)]
#fit_data = fit_data[,-c(2,3,4,5,6,7,8,9)]


y <- fit_data$HPV
X <- fit_data[,-dim(fit_data)[2]] #remove the OS and center column
X <- X[,-dim(X)[2]]


ind <- which(centers == centers_levels[1])
X1 <-X[ind,]
y1 <- y[ind]
ind <- which(centers == centers_levels[2])
X2 <-X[ind,]
y2 <- y[ind]
ind <- which(centers == centers_levels[3])
X3 <-X[ind,]
y3 <- y[ind]
ind <- which(centers == centers_levels[4])
X4 <-X[ind,]
y4 <- y[ind]

fit_data$c <- NULL
fit_data$center <-NULL