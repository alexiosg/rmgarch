/*################################################################################
##
##   R package rmgarch by Alexios Ghalanos Copyright (C) 2009-2013
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
#################################################################################*/
#ifndef _rdcc_H
#define _rdcc_H
#include <RcppArmadillo.h>
#include "rmdist.h"
RcppExport SEXP dccnormC1(SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP );
RcppExport SEXP dccnormC2(SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP );
RcppExport SEXP dccstudentC1(SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP );
RcppExport SEXP dccstudentC2(SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP );
RcppExport SEXP dcclaplaceC1(SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP );
RcppExport SEXP dcclaplaceC2(SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP );
RcppExport SEXP dccfilter(SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP );
RcppExport SEXP dccsimmvn(SEXP , SEXP , SEXP, SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP );
RcppExport SEXP dccsimmvt(SEXP , SEXP , SEXP, SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP );
RcppExport SEXP dccsimmvl(SEXP , SEXP , SEXP, SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP );
RcppExport SEXP fdccnormC1(SEXP , SEXP , SEXP , SEXP , SEXP , SEXP );
RcppExport SEXP fdccnormC2(SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP );
RcppExport SEXP fdccsimmvn(SEXP , SEXP , SEXP , SEXP , SEXP , SEXP, SEXP , SEXP , SEXP );
RcppExport SEXP dewmacov1(SEXP , SEXP , SEXP );
RcppExport SEXP dewmacov2(SEXP , SEXP );
RcppExport SEXP sewmacov1(SEXP , SEXP , SEXP );
RcppExport SEXP sewmacov2(SEXP , SEXP );
#endif
