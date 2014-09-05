#################################################################################
##
##   R package rmgarch by Alexios Ghalanos Copyright (C) 2008-2013.
##   This file is part of the R package rmgarch.
##
##   The R package rmgarch is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rmgarch is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################

ewmacov = function(X, lambda.v = 0.96, lambda.m = 0.96, diag.v = FALSE, vec.m = FALSE,
		use.m = TRUE, cor = FALSE)
{
	m = ncol(X)
	n = nrow(X)
	if(!use.m){
		if(diag.v){
			if(length(lamba.v)==1) lambda.v = diag(rep(lambda.v, m))
			if(length(lambda.v)<m) stop("\nwrong length for lambda.v (must be equal to no.cols X)")
			lambda.v = diag(lambda.v[1:m])
			ans = .Call("dewmacov2", X = coredata(X), vlambda = lambda.v, PACKAGE="rmgarch")
		} else{
			lambda.v = lambda.v[1]
			ans = .Call("sewmacov2", X = coredata(X), vlambda = lambda.v, PACKAGE="rmgarch")
		}
	} else{
		if(diag.v){
			if(length(lamba.v)==1) lambda.v = diag(rep(lambda.v, m))
			if(length(lambda.v)<m) stop("\nwrong length for lambda.v (must be equal to no.cols X)")
			lambda.v = diag(lambda.v[1:m])
		} else{
			lambda.v = lambda.v[1]
		}
		if(vec.m){
			if(length(lambda.m)==1) lambda.m = rep(lambda.m, m)
			if(length(lambda.m)<m) stop("\nwrong length for lambda.m (must be equal to no.cols X)")
			lambda.m = lambda.m[1:m]
		} else{
			lambda.m = rep(lambda.m, m)
		}
		if(diag.v){
			ans = .Call("dewmacov1", X = coredata(X), vlambda = lambda.v, mlambda = lambda.m, PACKAGE="rmgarch")
		} else{
			ans = .Call("sewmacov1", X = coredata(X), vlambda = lambda.v, mlambda = lambda.m, PACKAGE="rmgarch")
		}
	}
	if(cor) C = array(apply(ans$V, 3, function(x) cov2cor(x)), dim=c(m,m,n)) else C = NULL
	return(list(cov = ans$V, cor = R, mu = ans$M))
}

