#fits 
# \hat{f}(x) = \hat{g}_{1}(x_S)+\hat{g}_{2}(x_C)


backfitter = function(X, y, predictor, fitMethod, predictfcn, eps = .01, iter.max=10, verbose=TRUE,...){

  ######## check inputs
  # fitMethod 
  if(missing(fitMethod)){
    stop("Must pass fitMethod that accepts arguments 'X' and 'y'.")
  }else{
    fcn_args = names(formals(fitMethod))
    if(!("X" %in% fcn_args) || !("y" %in% fcn_args)){
      stop("fitMethod must accept arguments X and y.")
    }
  }
  
  # predictfcn checking
  # check for valid prediction routine...
  if(missing(predictfcn)){
    stop("Must pass predictfcn.")
  }
  else {
    fcn_args = names(formals(predictfcn))
    if (!("object" %in% fcn_args)){
      stop("predictfcn must have an 'object' argument. New X data is optionally passed using 'newdata'.")
    } 
  }

  # some useful constants.
	N = nrow(X)

	#order by the predictor
	xorder = order(X[, predictor])
	X = X[xorder, ]
	y = y[xorder]

  if(!is.numeric(predictor)){
    which_col = which(names(X)==predictor)
    Xc = X[, -which_col]
  }
	else{
	  Xc = X[, -predictor]
	}
	Xs = as.vector(X[, predictor])	

	#initialize
	g1_of_Xs = rep(0,N)
	g2_of_Xc = rep(0,N)
	current_g2 = NULL

  # supsmu  will condense ties in x, so length(new_g1) = length(unqiue(Xs))
  # needless to say, this is not ideal.
  times_to_repeat = as.numeric(table(Xs)) #remember Xs is sorted
  
	OneStep = function(){
		#do g2 first
		new_g2_mod = fitMethod(X=Xc, y=(y-g1_of_Xs))
		new_g2 = predictfcn(object=new_g2_mod, Xc) 
		new_g1 = supsmu(x=Xs, y=(y-new_g2))$y
    	new_g1 = rep(new_g1, times_to_repeat) #matches length of new_g2 now.
		return(list(new_g1=new_g1,new_g2=new_g2,new_g2_mod=new_g2_mod))
	}

	delta = Inf
	iter = 0
	while( delta > eps && iter < iter.max){
		#one iteration
		nextStep = OneStep()
		
		#compute delta
		delta = sum((nextStep$new_g1 - g1_of_Xs)^2) / sum(g1_of_Xs^2)
    	delta = delta + sum((nextStep$new_g2 - g2_of_Xc)^2) / sum(g2_of_Xc^2)  

		#update
		current_g2 = nextStep$new_g2_mod
		g1_of_Xs = nextStep$new_g1
		g2_of_Xc = nextStep$new_g2
		iter = iter + 1
		
		#print message
		if(verbose){
			cat(paste("iter",iter," delta: ", round(delta,3),sep=""))
			cat("\n")
		}
	}
	#leaving us with...
	bf_obj = list(g1_of_Xs=g1_of_Xs, g2_of_Xc = g2_of_Xc, g2_mod=current_g2,X=X,y=y,
				 predictor = predictor, iter=iter,delta=delta, 
         fitMethod=fitMethod, predictfcn=predictfcn)
	class(bf_obj) = "backfitter"
	return(bf_obj)
}
