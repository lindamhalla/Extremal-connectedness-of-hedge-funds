########################################################################################################
###                               Bivariate Husler-Reiss spectral density                            ###
########################################################################################################

########################################################################################################
### Log-likelihood and its 1st and 2nd derivatives (wrt its parameter \lambda) : 
### Needed for the NR step
########################################################################################################

#Below, w is the pseudo-observation (in the simplex) and alpha and beta the dependence parameters

#log-likelihood
log.h.HR <- function(x, lambda){
  (-(1/2))*log(2*pi) + log(lambda/(2*(-1 + x)^2*x)) - (2 + lambda^2*log(x/(1 - x)))^2/
    (8*lambda^2)
}

#1st derivative of the log-likelihood
#wrt \lambda 
HR.d1 <- function(x,lambda){
  (1/lambda^3 + 1/lambda - (1/4)*lambda*log(x/(1 - x))^2)
}

#2nd derivative of the log-likelihood
#wrt \lambda^2
HR.d2 <- function(x,lambda){
  (-((3 + lambda^2)/lambda^4) - (1/4)*log(x/(1 - x))^2)
}

#Maximize the log-likelihood (needed to set the initial parameter)
maxLikelihood.HR.spec <- function(data, init = NULL, maxit = 500, hess = F)
{
  
  lhood <- function(alpha)
  {
    out <- NULL
    for(i in 1:nrow(data)){
      out[i] <- log.h.HR(data[i,1],alpha)
    }
    return(-sum(out))
  }
  if(is.null(init))
  {
    init <- 0.5
  }
  
  min <- 1e-6
  max <- 1e6
  opt <- optim(init, lhood, method = "L-BFGS-B",lower = min, upper = max, control = list(maxit = maxit), hessian = hess)
  
  
  rtn <- list()
  rtn$par <- opt$par
  rtn$value <- opt$value
  rtn$convergence <- opt$convergence
  rtn$w <- data
  if(hess == F)
  {
    rtn$hessian <- 1
  }
  else
  {
    rtn$hessian <- opt$hessian
  }
  
  return(rtn)
}

########################################################################################################
### VGAM family
########################################################################################################

hr2d.vgam <- function (lshape1 = "loge", ishape1 = NULL, 
                       tol12 = 1e-04, parallel = FALSE, 
                       zero = NULL) 
{
  lshape1 <- as.list(substitute(lshape1))
  eshape1 <- link2list(lshape1)
  lshape1 <- attr(eshape1, "function.name")
  if (length(ishape1) && (!is.Numeric(ishape1))) 
    stop("bad input for argument 'ishape1'")
  if (!is.Numeric(tol12, length.arg = 1, positive = TRUE)) 
    stop("bad input for argument 'tol12'")
  
  new("vglmff", blurb = c("HR_2d distribution\n\n", "Links:    ", 
                          namesof("shape1", tag = FALSE), "\n"), 
      constraints = eval(substitute(expression({
        constraints <- cm.VGAM(matrix(1, M, 1), x = x, bool = .parallel, 
                               constraints = constraints)
        constraints <- cm.zero.VGAM(constraints, x = x, .zero, 
                                    M = M, predictors.names = predictors.names, M1 = 1)
      }), list(.parallel = parallel,.zero = zero))), infos = eval(substitute(function(...) {
        list(M1 = 1, Q1 = 1, expected = FALSE, multipleResponses = TRUE, 
             parameters.names = c("shape1"), lshape1 = .lshape1, 
             parallel = .parallel, zero= .zero)
      }, list(.parallel = parallel, .zero = zero, .lshape1 = lshape1))), 
      initialize = eval(substitute(expression({
        checklist <- w.y.check(w = w, y = y, Is.positive.y = TRUE, 
                               ncol.w.max = Inf, ncol.y.max = Inf, out.wy = TRUE, 
                               colsyperw = 1, maximize = TRUE)
        w <- checklist$w
        y <- checklist$y
        if (any((y <= 0) | (y >= 1))) stop("the response must be in (0, 1)")
        extra$multiple.responses <- TRUE
        extra$ncoly <- ncoly <- ncol(y)
        extra$M1 <- M1 <- 1
        M <- M1 * ncoly
        mynames1 <- param.names("shape1", ncoly)
        predictors.names <- c(namesof(mynames1, .lshape1, earg = .eshape1, tag = FALSE))[interleave.VGAM(M,M1 = M1)]
        if (!length(etastart)) {
          shape1.init <- .ishape1 
          shape1.init <- matrix(shape1.init, n, ncoly, 
                                byrow = TRUE)
          etastart <- cbind(theta2eta(shape1.init, .lshape1, 
                                      earg = .eshape1))[, interleave.VGAM(M, M1 = M1)]
        }
      }), list(.lshape1 = lshape1, .ishape1 = ishape1, .eshape1 = eshape1))), 
      linkinv = eval(substitute(function(eta,extra = NULL) {
        shape1 <- eta2theta(eta, .lshape1, 
                            earg = .eshape1)
        shape1
      }, list(.lshape1 = lshape1, .eshape1 = eshape1))), last = eval(substitute(expression({
        misc$link <- c(rep_len(.lshape1, ncoly))[interleave.VGAM(M, M1 = M1)]
        temp.names <- c(mynames1)[interleave.VGAM(M, M1 = M1)]
        names(misc$link) <- temp.names
        misc$earg        <- vector("list", M)
        names(misc$earg) <- temp.names
        misc$expected    <- FALSE
        misc$multiple.responses <- TRUE
        for (ii in 1:ncoly) {
          misc$earg[[M1 * ii]] <- .eshape1
        }
      }), list(.lshape1 = lshape1, .eshape1 = eshape1))),
      loglikelihood = eval(substitute(function(mu, y, w, residuals = FALSE, eta, extra = NULL, summation = TRUE) {
        shape1  <- eta2theta(eta, .lshape1, earg = .eshape1)
        ll.elts <- c(w) * log.h.HR(y, shape1)
        if (summation) sum(ll.elts) else ll.elts
      }, list(.lshape1 = lshape1, .eshape1 = eshape1))),
      vfamily = c("hr2d"), validparams = eval(substitute(function(eta, y, extra = NULL) {
        shape1 <- eta2theta(eta, .lshape1, earg = .eshape1)
        okay1 <- all(is.finite(shape1)) && all(0 < (shape1))
        okay1
      }, list(.lshape1 = lshape1, .eshape1 = eshape1))), 
      simslot = eval(substitute(function(object, nsim) {
        eta    <- predict(object)
        shape1 <- eta2theta(eta, .lshape1, earg = .eshape1)
        rmevspec(n=nsim*length(shape1), d=2, sigma=matrix(c(0,shape1^2,shape1^2,0),byrow=TRUE,ncol=2), 
                 model="hr")[,1]
      }, list(.lshape1 = lshape1, .eshape1 = eshape1))), deriv = eval(substitute(expression({
        shape1 <- eta2theta(eta, .lshape1, earg = .eshape1)
        dshape1.deta <- dtheta.deta(shape1, link = .lshape1, earg = .eshape1)
        d2shape1.deta2 <- d2theta.deta2(shape1, link = .lshape1, earg = .eshape1)
        dl.dshape1 <- HR.d1(y,shape1)
        dl.deta <- c(w) * cbind(dl.dshape1 * dshape1.deta)
        dl.deta
      }), list(.lshape1 = lshape1, .eshape1 = eshape1))), 
      # weight = eval(substitute(expression({
      #   ned2l.dshape11 <- -HR.d2(y,shape1)
      #   wz <- array(c(c(w) * cbind((ned2l.dshape11* dshape1.deta^2)-
      #                                (HR.d1(y,shape1)*d2shape1.deta2))), 
      #               dim = c(n, M/M1, 1))
      #   wz <- arwz2wz(wz, M = M, M1 = M1)
      #   wz
      # }), list(.lshape1 = lshape1, .eshape1 = eshape1, 
      #          .tol12 = tol12)))
      weight = eval(substitute(expression({
        
        ned2l.dshape11 <- 2*(1+(2/(shape1^2)))/(shape1^2)
        
        wz <- array(c(c(w) * ned2l.dshape11 * dshape1.deta^2),
                    dim = c(n, M / M1, 1))
        
        wz <- arwz2wz(wz, M = M, M1 = M1)
        wz
      }), list(.lshape1 = lshape1, .eshape1 = eshape1, .tol12 = tol12)))
  )
}

########################################################################################################
### Examples of use
########################################################################################################
if(!require(VGAM)) install.packages("VGAM") 
if(!require(mev)) install.packages("mev") 

library("VGAM")
library("mev")