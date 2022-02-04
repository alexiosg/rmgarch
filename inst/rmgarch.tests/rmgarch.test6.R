#################################################################################
##
##   R package rmgarch by Alexios Ghalanos Copyright (C) 2008-2022
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

# DCC Beta Model

rmgarch.test6a = function(cluster = NULL)
{
  # DCC under different specifications
  tic = Sys.time()
  require(xts)
  data(dji30retw)
  dji30retw = xts(coredata(dji30retw), as.Date(rownames(dji30retw)))
  benchmark = xts(rowMeans(dji30retw), index(dji30retw))
  Z = na.omit(cbind(benchmark, lag(dji30retw[, 1:15])))
  uspec = replicate(16, ugarchspec(mean.model=list(armaOrder=c(0,0)), variance.model=list(model = "gjrGARCH",variance.targeting=FALSE)))
  B1<-betaDCCFilter(y=Z[,1], x=Z[,-1], uspec, use.n = 1, fit=NULL, only.betas=FALSE)
  gspec<-	gogarchspec(
    mean.model = list(model = c("constant", "AR", "VAR")[1]),
    variance.model = list(model = "gjrGARCH", garchOrder = c(1,1), submodel = NULL, variance.targeting = FALSE),
    distribution.model = c("mvnorm", "manig", "magh")[1], ica = c("fastica", "radical")[2])
  B2<-betaGOGARCHFilter(y=Z[,1], x=Z[,-1], gspec, use.n = 1, fit=NULL, only.betas=FALSE, maxiter1=15000)

  plot(B1$betas[,1])
  lines(B2$betas[,1], col=2)

  plot(B1$betas[,2])
  lines(B2$betas[,2], col=2)

  plot(B1$betas[,3])
  lines(B2$betas[,3], col=2)

  plot(B1$betas[,4])
  lines(B2$betas[,4], col=2)

  caret:::RMSE(B1$fitted, Z[,1])
  caret:::RMSE(B2$fitted, Z[,1])

  toc = Sys.time()-tic
  cat("Elapsed:", toc, "\n")
  return(toc)
}
