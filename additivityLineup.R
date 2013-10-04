#nonparametric parametric bootstrap.

additivityLineup = function(backfit_obj, fitMethod, realICE, figs=10, colorvecfcn, usecolorvecfcn_inreal=F,
null_predictfcn,...){
# colorvec fcn is there to allow you to color null plots by the levels of a variable without
# introducing the indicator variable into the X matrix.  introducing it into the X matrix
# will let the algos use the indicator as a predictor -- prob not the behavior we want.
  
  ######## check inputs
  # backfit_obj
  if(class(backfit_obj)!="backfitter"){
    stop("'backfit_obj' is not of class 'backfitter'")
  }
  
  # fitMethod 
  if(missing(fitMethod)){
    stop("Must pass fitMethod that accepts arguments 'X' and 'y'.")
  }else{
    fcn_args = names(formals(fitMethod))
    if(!("X" %in% fcn_args) || !("y" %in% fcn_args)){
      stop("fitMethod must accept arguments X and y.")
    }
  }
  
  # predictfcn for nulls:
  if(missing(null_predictfcn)){
	if(is.null(realICE$predictfcn)){
		null_predictfcn = NULL
	}
	else{
		null_predictfcn = realICE$predictfcn	
	}
  }
  
  predictor = backfit_obj$predictor
  #some pre-processing regarding centering.
  arg_list = list(...)
  
  #frac_to_build is not allowed
  if(!is.null(arg_list$frac_to_build)){
    cat("Cannot specify  frac_to_build. Can specify frac_to_plot, which applies to the realICE, and frac_to_build for the null ICEs is then inferred to plot the same number of curves.")
    cat("\n")
  }
  
  centered = arg_list$centered
  if(is.null(centered)){
    centered = FALSE
  }
  centered_percentile = arg_list$centered_percentile
  if(is.null(centered_percentile) && centered==TRUE){ 
    centered_percentile = .01  #default in plot.ice
  }

  
  #and some more for frac_to_plot
  frac_to_build_null = 1

  if(!is.null(arg_list$frac_to_plot)){
    frac_to_plot = arg_list$frac_to_plot
    warning_msg = paste("'frac_to_plot' only applies to plotting 'realICE'.",
    "'frac_to_build' is set in null ICEs to ensure the same number of curves are plotted for null and real plots.",sep="\n")  
    warning(warning_msg)
    frac_to_build_null = nrow(realICE$ice_curves)*arg_list$frac_to_plot / nrow(backfit_obj$X)

	#fix indices to plot so that the ylim's can be constrained to only those
    #curves actually plotted.
    plot_points_indices = which(as.logical(rbinom(nrow(realICE$ice_curves), 1, frac_to_plot)))
    realICE$ice_curves = realICE$ice_curves[plot_points_indices, ]
	realICE$grid =  realICE$grid[plot_points_indices]
	realICE$Xice =  realICE$Xice[plot_points_indices,]
	realICE$xj =  realICE$xj[plot_points_indices]
	frac_to_plot = 1
  }

  #figure out min and max of real ice object -- depends on centering  
  if(centered){
    centering_vector = realICE$ice_curves[, ceiling(ncol(realICE$ice_curves) * centered_percentile + 0.00001)]
    rg = range(realICE$ice_curves - centering_vector) 
  }else{  #just the min and max.
    rg = range(realICE$ice_curves)
  }
  icecurve_min = rg[1]
  icecurve_max = rg[2]
  
	additive_fit = backfit_obj$g1_of_Xs+backfit_obj$g2_of_Xc
	additive_res = backfit_obj$y - additive_fit
	
	null_additive_fits = list()
	null_ices = list()
  
	for(i in 1:(figs-1)){
		response = additive_fit + sample(additive_res, size=length(additive_res), replace = F)
		new_fit = fitMethod(X=backfit_obj$X, y=response)
		null_additive_fits[[i]] = new_fit

		if(is.null(null_predictfcn)){ #no predictfcn found, use generic
	  		null_ices[[i]] = ice(new_fit, X=backfit_obj$X, predictor=predictor, y = backfit_obj$y,
								frac_to_build=frac_to_build_null)
		}else{
		  null_ices[[i]] = ice(new_fit, X=backfit_obj$X, predictor=predictor, y = backfit_obj$y, 
		                         frac_to_build=frac_to_build_null, predictfcn = null_predictfcn)
		}
    
		### keep track of min and max. 
		if(!is.null(centered) && centered==TRUE){  #keep track of range after centered

		  centering_vector = null_ices[[i]]$ice_curves[, ceiling(ncol(null_ices[[i]]$ice_curves) * centered_percentile + 0.00001)]
		  rg = range(null_ices[[i]]$ice_curves - centering_vector) #range for centered plot
		}
		else{  #regular pre-centered range
		  rg = range(null_ices[[i]]$ice_curves)      
		}

		#update min range and max range
		if(rg[1] < icecurve_min){
		    icecurve_min = rg[1]
		}
		if(rg[2] > icecurve_max){
		    icecurve_max = rg[2]
		}
		cat("Finished null ice ",i,"\n")
	} #end loop through null ices
  
	#graphics
	num_plot_cols = ceiling(figs/4)
	num_plot_rows = ceiling(figs/num_plot_cols)
	par(mfrow=c(num_plot_rows, num_plot_cols))
	par(cex=.3)
	par(mar=c(0.13,0.13,0.13,0.13))
    ylim = c(icecurve_min,icecurve_max)
  
   #argument list for the null plots.
   null_arg_list = arg_list
   null_arg_list$ylim = ylim
   null_arg_list$frac_to_plot = 1
  
	#randomly place the real plot somewhere...
	where_to_place = sample(1:figs,1)
	plotted_truth = FALSE
	for(i in 1:figs){

		if(plotted_truth){
			idx_to_plot = i-1
		}else{
			idx_to_plot = i
		}

		if((!plotted_truth) && (i==where_to_place)){
			if(usecolorvecfcn_inreal){
				colors = colorvecfcn(realICE)	
				plot(realICE,ylim = ylim,colorvec=colors,...)
			}else{
				plot(realICE,ylim = ylim,...)
			}
			plotted_truth = TRUE
		}
		else{
      		null_arg_list$x = null_ices[[idx_to_plot]] #generic plot has argument 'x'
			if(!missing(colorvecfcn)){
				colors = colorvecfcn(null_ices[[idx_to_plot]])
				null_arg_list$colorvec = colors
			}
      		do.call(plot, null_arg_list) 
		}
	}
  al_obj = list(location=where_to_place, null_additive_fits = null_additive_fits, 
                        null_ices = null_ices, frac_to_build_null = frac_to_build_null)
  class(al_obj) = "additivityLineup"  
	invisible(al_obj)
}
