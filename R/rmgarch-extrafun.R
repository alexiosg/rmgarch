.information.test <- function(LLH, nObs, nPars) 
{
    AIC = (-2 * LLH)/nObs + 2 * nPars/nObs
    BIC = (-2 * LLH)/nObs + nPars * log(nObs)/nObs
    SIC = (-2 * LLH)/nObs + log((nObs + 2 * nPars)/nObs)
    HQIC = (-2 * LLH)/nObs + (2 * nPars * log(log(nObs)))/nObs
    informationTests = list(AIC = AIC, BIC = BIC, SIC = SIC, 
                            HQIC = HQIC)
    return(informationTests)
}