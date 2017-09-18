#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP ArrayDiag(SEXP, SEXP);
extern SEXP copuladccsimmvn(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP copuladccsimmvt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP copulaNormalC1(SEXP, SEXP);
extern SEXP copulaNormalC2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP copulaStudentC1(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP copulaStudentC2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Cov2Cor(SEXP, SEXP);
extern SEXP Cov2Res(SEXP, SEXP, SEXP);
extern SEXP dcclaplaceC1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dcclaplaceC2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dccnormC1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dccnormC2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dccsimmvl(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dccsimmvn(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dccsimmvt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dccstudentC1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dccstudentC2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fdccnormC1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fdccnormC2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fdccsimmvn(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gogarchcksigma(SEXP, SEXP);
extern SEXP gogarchCor(SEXP, SEXP);
extern SEXP gogarchCov(SEXP, SEXP);
extern SEXP gogarchcssigma(SEXP, SEXP);
extern SEXP gogarchSigma(SEXP, SEXP);
extern SEXP radicalrot(SEXP, SEXP);
extern SEXP tvbetacokurt(SEXP, SEXP, SEXP);
extern SEXP tvbetacoskew(SEXP, SEXP, SEXP);
extern SEXP tvbetacovar(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"ArrayDiag",        (DL_FUNC) &ArrayDiag,         2},
    {"copuladccsimmvn",  (DL_FUNC) &copuladccsimmvn,  10},
    {"copuladccsimmvt",  (DL_FUNC) &copuladccsimmvt,  10},
    {"copulaNormalC1",   (DL_FUNC) &copulaNormalC1,    2},
    {"copulaNormalC2",   (DL_FUNC) &copulaNormalC2,    8},
    {"copulaStudentC1",  (DL_FUNC) &copulaStudentC1,   5},
    {"copulaStudentC2",  (DL_FUNC) &copulaStudentC2,   9},
    {"Cov2Cor",          (DL_FUNC) &Cov2Cor,           2},
    {"Cov2Res",          (DL_FUNC) &Cov2Res,           3},
    {"dcclaplaceC1",     (DL_FUNC) &dcclaplaceC1,      8},
    {"dcclaplaceC2",     (DL_FUNC) &dcclaplaceC2,      9},
    {"dccnormC1",        (DL_FUNC) &dccnormC1,         8},
    {"dccnormC2",        (DL_FUNC) &dccnormC2,         9},
    {"dccsimmvl",        (DL_FUNC) &dccsimmvl,        10},
    {"dccsimmvn",        (DL_FUNC) &dccsimmvn,        10},
    {"dccsimmvt",        (DL_FUNC) &dccsimmvt,        10},
    {"dccstudentC1",     (DL_FUNC) &dccstudentC1,      8},
    {"dccstudentC2",     (DL_FUNC) &dccstudentC2,      9},
    {"fdccnormC1",       (DL_FUNC) &fdccnormC1,        6},
    {"fdccnormC2",       (DL_FUNC) &fdccnormC2,        7},
    {"fdccsimmvn",       (DL_FUNC) &fdccsimmvn,        9},
    {"gogarchcksigma",   (DL_FUNC) &gogarchcksigma,    2},
    {"gogarchCor",       (DL_FUNC) &gogarchCor,        2},
    {"gogarchCov",       (DL_FUNC) &gogarchCov,        2},
    {"gogarchcssigma",   (DL_FUNC) &gogarchcssigma,    2},
    {"gogarchSigma",     (DL_FUNC) &gogarchSigma,      2},
    {"radicalrot",       (DL_FUNC) &radicalrot,        2},
    {"tvbetacokurt",     (DL_FUNC) &tvbetacokurt,      3},
    {"tvbetacoskew",     (DL_FUNC) &tvbetacoskew,      3},
    {"tvbetacovar",      (DL_FUNC) &tvbetacovar,       3},
    {NULL, NULL, 0}
};

void R_init_rmgarch(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}