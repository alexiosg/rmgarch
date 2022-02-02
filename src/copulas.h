/*################################################################################
##
##   R package rmgarch by Alexios Galanos Copyright (C) 2009-2013
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
#ifndef _copulas_H
#define _copulas_H
#include "rmdist.h"
#include <RcppArmadillo.h>
RcppExport SEXP copulaNormalC1(SEXP , SEXP);
RcppExport SEXP copulaNormalC2(SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP);
RcppExport SEXP copulaStudentC1(SEXP , SEXP , SEXP , SEXP , SEXP);
RcppExport SEXP copulaStudentC2(SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP);
RcppExport SEXP copuladccsimmvn(SEXP , SEXP , SEXP , SEXP, SEXP , SEXP , SEXP , SEXP , SEXP , SEXP );
RcppExport SEXP copuladccsimmvt(SEXP , SEXP , SEXP , SEXP, SEXP , SEXP , SEXP , SEXP , SEXP , SEXP );
#endif
