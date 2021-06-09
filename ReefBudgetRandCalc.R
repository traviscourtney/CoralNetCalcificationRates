## Modified from SUPPLEMENTARY MATERIAL for Possolo, Merkatas, & Bodnar (2019)
## FILE:              AsymmetricUncertainties2019Jun26-RCode.R
## AUTHORS:           Antonio Possolo, Christos Merkatas, Olha Bodnar
## MODIFICATION DATE: 2019 Jun 26
ReefBudgetRandCalc = function (x, lwr, upr, n,
                          pX=0.5, pL=0.0000001, pR=0.9999999,
                          wX=1, wL=1, wR=1)
{
  ## INPUTS
  ## x  = Measured value   : x is the 100*pX percentile
  ## lwr = Left uncertainty : x - uL is the 100*pL percentile
  ## upr = Right uncertainty: x + uR is the 100*pR percentile
  ## n = number of values to generate from distribution
  ## wX, wL, wR = Weights for the errors made when attempting to
  ## reproduce x, x-uL, and x+uR as percentiles of a skew-normal
  ## distribution 
  ## OUTPUT
  ## Vector with the values of xi, omega, and alpha for the best
  ## fitting skew-normal distribution
  uL=x-lwr
  uR=upr-x
  if (!require(sn)) {stop("R package 'sn' is not installed") }
  if (any(c(wX, wL, wR) < 0)) {
    stop(paste("ERROR in parSkewNormal: Weights wL, wX, and wR",
               "must all be positive")) }
  if (!((pL < pX) & (pX < pR))) {
    stop(paste("ERROR in parSkewNormal:",
               "Probabilities must be such that",
               "pL < pX < pR")) }
  fSkewNormal = function (theta, L, X, R, pL, pX, pR, wL, wX, wR)
  { xi = theta[1]; omega = theta[2]; alpha = theta[3]
  return(sum(c(wL,wX,wR) *
               (qsn(c(pL,pX,pR),
                    xi=xi, omega=omega, alpha=alpha, 
                    solver="RFB") - c(x-uL,x,x+uR))^2)) }
  o = try(optim(par=c(x, if (abs(pR-pL) < 0.75) {(uL+uR)/2
  } else {(uL+uR)/4}, 2),
  fn=fSkewNormal, method="Nelder-Mead",
  L=x-uL, X=x, R=x+uR, pL=pL, pX=pX, pR=pR,
  wL=wL, wX=wX, wR=wR))
  if (class(o) == "try-class") {
    stop("Optimization failed")
  } else {
    theta = o$par
    names(theta) = c("xi", "omega", "alpha")
    dist=rsn(n=n, xi=theta[1], omega=theta[2], alpha=theta[3])
    return(dist) }
}