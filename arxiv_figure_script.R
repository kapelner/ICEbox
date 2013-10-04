#######################################################################
# This R script generates the figures in 
#   "Peeking Inside the Black Box: Visualizing Statistical Learning
#    with Plots of Individual Conditional Expectation,"
# available on arXiv at:
# <arxiv link here>. 
#
# The stable version of the ICEbox package is available on CRAN, 
# and all examples below will work with this version.  Beta
# versions of the package with additional features are available at
# <github link here>.
#
# Note that examples pertaining to the depression dataset are omitted, 
# as this data cannot currently be made public.
#######################################################################

### load R packages
library(ICEbox)
library(randomForest)
library(gam)
library(gbm)
library(nnet)
library(missForest)
library(MASS) #has Boston Housing data, Pima


#######################################################################
############ Section: The ICE Toolbox #################################
#######################################################################
# load the BH data:
data(Boston)
X = Boston
y = X$medv
X$medv = NULL

# build an RF with default settings:
rf_mod = randomForest(X, y)


#################### ICE Procedure ####################################
#make frac_to_build < 1 for fewer curves but faster computation
rf.ice = ice(rf_mod, X, y, predictor = "age", frac_to_build = 1) 

# Note we can plot Friedman's PDP from Section 2
# by making the individual curves in the ICE white:
plot(rf.ice, x_quantile = TRUE, plot_pdp = TRUE, 
		colorvec = rep("white", nrow(X)), plot_orig_pts_preds = FALSE)

### plot ICE:
# Make frac_to_plot < 1 to plot a fraction of the curves built.  
# This an make the plot less cluttered sometimes.
plot(rf.ice, x_quantile = TRUE, plot_pdp = TRUE, frac_to_plot = 1) 

################### c-ICE #############################################
#just a matter of saying 'centered=TRUE'
plot(rf.ice, x_quantile = TRUE, plot_pdp = TRUE, centered = TRUE)

################### d-ICE #############################################
#create the object:
rf.dice = dice(rf.ice)

# plot the d-ICE:
plot(rf.dice, x_quantile = T)

################### Visualizing a second feature ######################
# We investigate the c-ICE  by coloring it by the "rm" variable. 
# First create an indicator variable:
rf.ice$Xice$I_rm = ifelse(rf.ice$Xice$rm > 6.2, 1, 0)  

# then plot using 'color_by'.
plot(rf.ice, frac_to_plot = 1, centered = TRUE, prop_range_y = TRUE,  
		x_quantile = T, plot_orig_pts_preds = T, color_by = "I_rm")
#######################################################################




#######################################################################
############ Section: Simulations #####################################
#######################################################################

#################### Additivity Assessment ############################
#function that generates simulated data:
additivity_ex_sim = function(n,seednum=NULL){
	if(!is.null(seednum)){
		set.seed(seednum)
	}
	p = 2
	X = as.data.frame(matrix(runif(n * p, -1, 1), ncol = p))	
	colnames(X) = paste("x_", 1 : p, sep = "")
	bbeta = c(1,1)

	y = bbeta[1] * X[,1]^2  + bbeta[2] * X[,2]
	y = y + rnorm(n)
	Xy = as.data.frame(cbind(X, y))
	return(list(Xy=Xy,X=X,y=y))
}

# generate data:
additivity_ex_data = additivity_ex_sim(1000)
Xy = additivity_ex_data$Xy
X  = additivity_ex_data$X
y  = additivity_ex_data$y

# build gam with possible interactions:
gam_mod = gam(y~s(x_1)+s(x_2)+s(x_1*x_2),data=Xy)   

# build ICE and d-ICE:
gam.ice = ice(gam_mod, X, predictor = 1, frac_to_build = 1) 
gam.dice = dice(gam.ice)

# plot the ICE plot with pdp, and d-ICE with dpdp
plot(gam.ice, x_quantile = F, plot_pdp = T, frac_to_plot = 1)  
plot(gam.dice, x_quantile = F, plot_dpdp = T, frac_to_plot = 1) 

#################### Finding Interactions #############################
# function that generates simulated data:
interaction_ex_sim = function(n,seednum=NULL){
	if(!is.null(seednum)){
		set.seed(seednum)
	}

	p = 3
	X = matrix(runif(n * p, -1, 1), ncol = p)		
	colnames(X) = paste("x_", 1 : p , sep = "")

	#coefficients dependent on level of x_3:
	bbeta1 = as.matrix(c(0.2, 5, 0))
	bbeta2 = as.matrix(c(0.2, -5, 0))

	#now generate y
	y = array(NA, n)
	for (i in 1 : n){
		if (X[i, 3] >= 0 ){
			y[i] = X[i, ] %*% bbeta1
		} else {
			y[i] = X[i, ] %*% bbeta2
		}
	}
	Xy = as.data.frame(cbind(X, y))
	X = as.data.frame(X)
	return(list(Xy=Xy,X=X,y=y))
}

# generate data:
interaction_ex_data = interaction_ex_sim(1000)
Xy = interaction_ex_data$Xy
X  = interaction_ex_data$X
y  = interaction_ex_data$y

# build stochastic gradient boosting with cross validated number of trees and interaction depth=3
gbm_mod = gbm(y ~ ., data = Xy, n.tree = 500, 
				interaction.depth = 3, shrinkage = 0.1, cv.folds = 5)
ntree = gbm.perf(gbm_mod, method = "cv")

# Create and plot ICE object:
# To predict with the cross-validated number of trees, we pass a custom predict function
gbm.ice = ice(gbm_mod, X, predictor = "x_3", 
			predictfcn = function(object, newdata){predict(object, newdata, n.tree = ntree)},
			frac_to_build = 1)

#plot only 5% of curves with quantiles, actual pdp, and original points. 
plot(gbm.ice, x_quantile = F, plot_pdp = T, frac_to_plot = 0.05)

# Create and plot d-ICE object:
gbm.dice = dice(gbm.ice)
plot(gbm.dice, plot_orig_pts_deriv = FALSE)

#################### Extrapolation Detection ##########################
#function that generates simulated data:
extrap_ex_sim = function(N,seednum=NULL){
	
	if(!is.null(seednum)){
		set.seed(seednum)
	}

	#helper that simulates one observation:
	sim1 = function(){
		s = runif(1)
		if(s <  (1/3)){
			x1 = runif(1,min=-1,max=0)
			x2 = runif(1,min=-1,max=0)
		}
		else{
			if( s < (2/3)){
				x1 = runif(1,min=-1,max=0)
				x2 = runif(1,min=0,max=1)
			}
			else{
				x1 = runif(1,min=0,max=1)
				x2 = runif(1,min=-1,max=0)
			}
		}
		return(cbind(x1,x2))
	} #end of helper fcn.

	#parameters:	
	b1=10; b2=1; sd=.1
	X = matrix(NA,ncol=2,nrow=N)

	for(i in 1:N){X[i,] = sim1()}

	Y = b1*X[,1]^2 + b2*(X[,2]>0)+rnorm(N,sd=sd)
	df = as.data.frame(cbind(Y,X))
	names(df) = c("Y","x1","x2")
	df
}

# generate data:
extrap_ex_data = extrap_ex_sim(1000)
X = extrap_ex_data[,2:3] #predictors are second and third columns
# fit a rf:
rf_mod = randomForest(Y~.,extrap_ex_data)

# Create an ICE object:
rf.ice_extrap = ice(rf_mod, X = X, predictor="x1", frac_to_build=1)

# Set up 'color_by' to visualize extrapoloations in X-space.
# x2_indic = 1 <==> no extrapolation 
rf.ice_extrap$Xice$x2_indic = ifelse(rf.ice_extrap$Xice$x2>0,0,1) 

# full ICE plot:
plot(rf.ice_extrap, plot_pdp=FALSE, color_by="x2_indic") 
#######################################################################



#######################################################################
############ Section: Real Data #######################################
#######################################################################

##################### Depression Clinical Data ########################
# Examples pertaining to the depression dataset are omitted, 
# as this data cannot currently be made public.


##################### White Wine ######################################
#load the data:
data(WhiteWine) #part of the ICEbox package.
X = WhiteWine[,-1]
y = WhiteWine$quality

predictors = names(X)
N = nrow(X)
frac_to_plot = 300/N; frac_to_plot = min(frac_to_plot, 1)   #suitable frac_to_plot

# neural nets require standardization of the X matrix:
X_std = scale(x = X, center = T, scale = T)
X_center = attributes(X_std)$`scaled:center`  
X_scale = attributes(X_std)$`scaled:scale`  

## fit the nnet:
nnet_mod = nnet(x = X_std, y = as.matrix(y), size = 3, 
			maxit = 500, decay = 5e-4, linout = ifelse(is.factor(y), F, T))


## ICE for the pH predictor:
# We need a custom predict function because of the standardization.  We standardize
# "new" X data by the training data. Note all ice_curves/ plots are returned to the scale of 
# the original data.
nn.ice = ice(nnet_mod, X=X, predictor = "pH", 
             predictfcn = function(object, newdata){
                            newdata_std = scale(newdata, center = X_center, scale = X_scale)
                            predict(object, newdata_std)
                           }, 
			y=y)

## ICE 
# create an indicator for alcohol content to pass to 'color_by':
nn.ice$Xice$al_ind = ifelse(nn.ice$Xice$alcohol > 10, 1, 0)

# Make a c-ICE plot, with frac_to_plot set so that ~300 observations are plotted.
frac_to_plot = 300 / nrow(WhiteWine)
plot(nn.ice, centered=TRUE, centered_percentile=0.01, frac_to_plot=frac_to_plot, color_by = "al_ind", x_quantile=TRUE)

## d-ICE
nn.dice = dice(nn.ice)
plot(nn.dice, x_quantile=TRUE, frac_to_plot = frac_to_plot, 
				plot_sd = TRUE, plot_dpdp = TRUE, color_by = "al_ind")


################### Diabetes Classification in Pima Indians ###########
data(Pima.te)
X = Pima.te[,-8]
y = Pima.te[, 8]

# build rf:
pima_rf = randomForest(x = X, y = y)

## Create an ICE object for the predictor "skin":
# For classification we plot the centered log-odds. If we pass a predict
# function that returns fitted probabilities, setting logodds = TRUE instructs
# the function to set each ice curve to the centered log-odds of the fitted 
# probability.
pima.ice = ice(pima_rf, X = X, predictor = "skin", logodds = TRUE,
                    predictfcn = function(object, newdata){
							predict(object, newdata, type = "prob")[, 2]
					}
			 )

# make a c-ICE plot:
plot(pima.ice, x_quantile = TRUE, centered = TRUE)

## make a d-ICE object and plot it.
pima.dice = dice(pima.ice)
plot(pima.dice, x_quantile = TRUE)


#######################################################################
############ Section: Additivity Lineup ###############################
#######################################################################

### 1) source the functions:
source("additivityLineup.R"); source("backfitter.R")

### 2) load the data:
data(WhiteWine) 
X = WhiteWine[,-1]
y = WhiteWine$quality
predictors = names(X)
pred_name = "pH" 

### 3) preliminary functions:
# The standardization means we need to be careful about standardizing
# with and without pH.
# all:
X_std = scale(x = X, center = T, scale = T)
X_center = attributes(X_std)$`scaled:center`  
X_scale = attributes(X_std)$`scaled:scale`

# without pH (pH is 9th col):
X_std_noPH = scale(x = X[,-9], center = T, scale = T)
X_center_noPH = attributes(X_std_noPH)$`scaled:center`  
X_scale_noPH = attributes(X_std_noPH)$`scaled:scale`  

# Simlarly, we need a predict function for backfitting as per null (no pH):
nn_predictfcn = function(object, newdata){ #for backfitting, one col removed from std
  X_std = scale(newdata, center = .GlobalEnv$X_center_noPH, scale = .GlobalEnv$X_scale_noPH)
  predict(object, X_std)
}

# ... and a fit function for backfitting as per null (no pH):
nn_Fit = function(X,y){ #for backfitting.
  X_std = scale(X, center = .GlobalEnv$X_center_noPH, scale = .GlobalEnv$X_scale_noPH)
  nnet(x = X_std, y = as.matrix(y), size = 3, maxit = 500, decay = 5e-4, linout = ifelse(is.factor(y), F, T))
}

# ... and finally a fit using all X.
nn_FitFull= function(X,y){ #used for fitting full X to y whose cond. exp. is additive in g(pH) and h(other predictors).
  X_std = scale(X, center = .GlobalEnv$X_center, scale = .GlobalEnv$X_scale)
  nnet(x = X_std, y = as.matrix(y), size = 3, maxit = 500, decay = 5e-4, linout = ifelse(is.factor(y), F, T))
}


### 4) make the real ICE object.
# fit the nnet: 
nnet_mod = nnet(x = X_std, y = as.matrix(y), size = 3, 
			maxit = 500, decay = 5e-4, linout = ifelse(is.factor(y), F, T))

# create the ice object, with the predict function using the full scaling.
nn_wine_ice = ice(nnet_mod, X=X, y=y, predictor = pred_name, 
                    predictfcn = function(object, newdata){
                          newdata_std = scale(newdata, center = .GlobalEnv$X_center, scale = .GlobalEnv$X_scale)
                          predict(object, newdata_std)
                    } 
                 )

### 5) backfitting.
bf_pH = backfitter(X=X,y=y, predictor="pH", eps=.001, fitMethod=nn_Fit,
                           predictfcn = nn_predictfcn, iter.max=30)

### 6) create the lineup.
# this function generates the coloring according to alcohol content:
colorFcn = function(ice_obj){
	ifelse(ice_obj$Xice$alcohol > 10, "RED", "GREEN")
}

# Lineup:
# In the paper, the real plot is in the bottom right position, but this is randomly
# chosen on each call.
alu_pH = additivityLineup(bf_pH, fitMethod=nn_FitFull, figs=12, realICE=nn_wine_ice, 
                          centered=TRUE, x_quantile=TRUE, frac_to_plot=.1,
						  colorvecfcn=colorFcn, usecolorvecfcn_inreal=TRUE, plot_orig_pts=FALSE)
