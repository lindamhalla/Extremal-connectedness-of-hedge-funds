library(quantreg)
library(evd)
library(VGAM)
library(lubridate)

#load the dataframe containing the transformed negative returns (to the Frechet scale) and the set of predictors 
#for the extremal connectedness modeling
load("df.multivariate.Rdata")

#source the script with the needed functions to fit a Husler--Reiss spectral density
source("spectral_HR.R")


################################################
#number of pairs of strategies
unique.substrat <- unique(df.multivariate$strategy)
nb.pairs <- combn(length(unique.substrat),2)

#remove Fund of funds (3rd strategy)
index.fund.of.funds <- NULL
for(j in 1:ncol(nb.pairs)){
  if((nb.pairs[1,j] ==3)|(nb.pairs[2,j] ==3)){
    index.fund.of.funds <- c(index.fund.of.funds,j)
  }
}

id.pairs.finished <- (1:78)[-index.fund.of.funds]
################################################


################################################
###### model the extremal dependence of the i-th pair by considering observations 
#with a radial component exceeding its thd-th quantile
fct.boot <- function(i,thd){
  y12.std         <- NULL
  pair            <- nb.pairs[,i]
  unique.substrat <- unique(df.multivariate$strategy)

  dat1 <- df.multivariate[which((unique.substrat[pair[1]]==df.multivariate$strategy)),]
  dat2 <- df.multivariate[which((unique.substrat[pair[2]]==df.multivariate$strategy)),]
  
  dat.all <- sort(intersect(dat1$date,dat2$date))
  for(j in 1:length(dat.all)){
    id.time1 <- which(dat1$date==as.Date(dat.all[j],origin = "1970-01-01"))
    id.time2 <- which(dat2$date==as.Date(dat.all[j],origin = "1970-01-01"))
    dat1.j   <- dat1[id.time1,]
    dat1.j   <- dat1.j[order(dat1.j$frech,decreasing = TRUE),]
    dat2.j   <- dat2[id.time2,]
    dat2.j   <- dat2.j[order(dat2.j$frech,decreasing = TRUE),]
    id.time  <- min(nrow(dat1.j),nrow(dat2.j)) 
    if(id.time>0)
      y12.std   <- rbind(y12.std ,cbind(dat1.j$frech[1:id.time],
                                        dat2.j$frech[1:id.time],
                                        dat1.j$VIX[1:id.time],
                                        dat1.j$MSCI[1:id.time],
                                        dat1.j$EPU[1:id.time],
                                        dat1.j$stress[1:id.time],
                                        dat1.j$date[1:id.time]))
  }
  
  y1.std <- y12.std[,1]
  y2.std <- y12.std[,2]
  
  ################################
  db       <- data.frame("y1"=y1.std,"y2"=y2.std,"VIX"=y12.std[,3], "MSCI"=y12.std[,4],
                         "EPU"=y12.std[,5], "stress"=y12.std[,6], "time"=as.Date(y12.std[,7], origin="1970-01-01"))
  #keep rows with non-NA values
  db<-db[rowSums(is.na(db)) == 0, ]
  
  #fit only joint exceedances
  R          <- db$y1 + db$y2
  db$R       <- R
  quantreg.R <- rq(R~time, tau=thd, na.action = na.omit, data = db)
  uu         <- quantreg.R$fitted.values
  z.data     <- data.frame("y"=(db$y1/(db$y1+db$y2))[R>uu], 
                           "VIX"=db$VIX[R>uu],
                           "MSCI"=db$MSCI[R>uu],
                           "EPU"=db$EPU[R>uu],
                           "stress"=db$stress[R>uu],
                           "time"=db$time[R>uu])
  par.init   <- maxLikelihood.HR.spec(cbind(z.data$y,1-z.data$y))$par
  
  pair1.hr.vgam.1 <- VGAM::vgam(y~sm.ps(VIX, outer.ok=TRUE)+
                                  sm.ps(MSCI, outer.ok=TRUE)+
                                  sm.ps(EPU, outer.ok=TRUE)+
                                  sm.ps(stress, outer.ok=TRUE),family=hr2d.vgam(ishape1=par.init),
                                data=z.data,trace=TRUE,Maxit.outer=100,maxit=50)
  
  if(!pair1.hr.vgam.1@ospsslot$magicfit$gcv.info$fully.converged)
    pair1.hr.vgam.1 <- VGAM::vgam(y~VIX+MSCI+EPU+stress,family=hr2d.vgam(ishape1=par.init),
                                  data=z.data,trace=TRUE,Maxit.outer=100,maxit=50, outer.ok=TRUE)
  
  return(pair1.hr.vgam.1)
}

#run the fit over all pairs
vgam.pair.2keep <- list()
for(i in id.pairs.finished){
  vgam.pair.2keep[[i]] <- fct.boot(i,thd = .95)
}

###remove NULL elements in list of fitted vgam (smooth)
vgam.pair.2keep       <- vgam.pair.2keep[!sapply(vgam.pair.2keep, is.null)]

############### fix time period over which to predict all tail coef: 1994-01-01 to 2017-05-31
load("time_cov_2predict.RData")

#########################################################################
############### Compute the fitted and predicted extremal coefficients 
############### over the time period 1994-01-01 to 2017-05-31
#########################################################################
coef.ext.hr.1 <- list()
for(i in id.pairs.finished){
  #predict for all possible dates and hence values of covariates in time.cov.2.predict
  coef.ext.hr.1[[i]] <- vector(mode="numeric",length=nrow(time.cov.2.predict))
  
  blabla         <- VGAM::loglink(predictvglm(vgam.pair.2keep[[i]],newdata=time.cov.2.predict,
                                              se=FALSE), inverse=TRUE)
  for(l in 1:nrow(time.cov.2.predict)){
    coef.ext.hr.1[[i]][l] <- 2-2*abvevd(dep=blabla[l],
                                        model="hr",plot=FALSE)
  }
  
}

#remove pairs with fund of funds (null as we didn't fit the dependence)
coef.ext.hr.1 <- coef.ext.hr.1[!sapply(coef.ext.hr.1, is.null)]

#########################################################################
### saving data
#########################################################################

tail.coef.pairs <- list()

tail.coef.pairs$strategy  <- unique.substrat
tail.coef.pairs$pairs     <- nb.pairs[,-index.fund.of.funds]
tail.coef.pairs$tail.coef <- list()

for(k in 1:length(coef.ext.hr.1)){
  tail.coef.pairs$tail.coef[[k]] <- data.frame("time"=time.cov.2.predict$time,"ext.coef"=coef.ext.hr.1[[k]])
}

tail.coef.pairs$tail.coef <- tail.coef.pairs$tail.coef[!sapply(tail.coef.pairs$tail.coef, is.null)]

##############################################################################
############### Plots of the mean (unconditional) chi coef. ##################
##############################################################################

library(igraph) 
library(RColorBrewer)
library(fields)

substrat.finished <- 1:13
w <- as.vector(unlist(lapply(tail.coef.pairs$tail.coef, function(x) mean(x[,2]))))
w[(w<quantile(w,.5))] <- NA

vec.col <- colorRampPalette(brewer.pal(11,"Spectral"))(length(w[!is.na(w)]))[rank(w, na.last = NA)]

############ 
set.seed(413)
blue_dag <- make_full_graph(n=12,directed = FALSE)

bla.all <- rep(NA, length(w))
bla.all[!is.na(w)] <- vec.col[1:length(w[!is.na(w)])]

plot(blue_dag, edge.width = exp(5*w)/8, edge.color=bla.all,
     vertex.color="black", vertex.size=5,
     vertex.frame.color="black", vertex.label.color="black", vertex.label.dist=2,
     main=" ",
     vertex.label=unique.substrat[substrat.finished[-3]],
     vertex.label.cex=.9, layout=layout_in_circle(blue_dag,c(4,9,2,6,8,3,1,11,7,12,10,5)),
     mark.groups=list(c(4,9,2), c(6,8,3), c(1,11), c(7,12,10,5)), mark.col=rep("grey85",4), mark.border=NA)
image.plot(legend.only=T, zlim=range(w,na.rm=TRUE), 
           col=colorRampPalette(brewer.pal(11,"Spectral"))(length(w)),
           horizontal = TRUE, 
           legend.lab=expression(paste("E{", chi, "(", bold(x)[t], ") | ", t%in% tilde(T), "=[1994-01-10,2017-05-31]}")),
           legend.line = 2.3)

##############################################################################
############### compute the ECoVaR for all pairs
##############################################################################
estimate.C.EV <- function(u,v,lambda){ #A is the function returning the Pickands dependence function
  u1   <- log(u)
  u2   <- log(v)
  w    <- u1/(u1+u2)
  A.HR <- function(w, lambda){
    (1-w)*pnorm(lambda+(log((1-w)/w)/(2*lambda))) + 
      w*pnorm(lambda+(log(w/(1-w))/(2*lambda)))
  }
  return(exp((u1+u2)*A.HR(w, 1/lambda)))
}

fct.uniroot <- function(x,alpha, beta, lambda){
  1-beta-x+ estimate.C.EV(beta,x, lambda) - (1-alpha)*(1-beta)
}

####
ecovar_mat <- matrix(NA, ncol=length(tail.coef.pairs$tail.coef), nrow=nrow(tail.coef.pairs$tail.coef[[1]]))
for(k in 1:length(tail.coef.pairs$tail.coef)){
  time     <- tail.coef.pairs$tail.coef[[k]]$time
  chi.pair <- tail.coef.pairs$tail.coef[[k]]$ext.coef
  
  ### Need to compute the Husler--Reiss parameter from the extremal coefficient
  lambda.pair <- 1/qnorm(1-(chi.pair/2))
  
  ### ECoVaR level
  alpha=0.975
  ### VaR level
  beta=0.975
  
  ecovar <- NULL
  for(i in 1:length(lambda.pair)){
    u_alpha <- uniroot(fct.uniroot,interval = c(0.1,0.99999),alpha, beta, lambda.pair[i])$root
    ecovar  <- c(ecovar, u_alpha)
  }
  ecovar_mat[,k] <- ecovar
}

ecovar_list <- list("ecovar"=ecovar_mat,
                    "time"=tail.coef.pairs$tail.coef[[1]]$time,
                    "strategy"=tail.coef.pairs$strategy,
                    "pairs"=tail.coef.pairs$pairs)
col_names <- NULL
for(i in 1:ncol(ecovar_list$pairs)){
  col_names <- c(col_names, paste0(ecovar_list$strategy[ecovar_list$pairs[1,i]],"_vs_",
                                   ecovar_list$strategy[ecovar_list$pairs[2,i]]))
}
colnames(ecovar_mat) <- col_names

