#-------------------------------------------------------------------------------------------
#
# Code supplement for article
# 
# Choosy grazers: Dependence of forage selection by three cattle breeds on plant traits
# by Caren M. Pauler, Johannes Isselstein, Matthias Suter, Joël Bérard, Thomas Braunbeck, Manuel K. Schneider
#
# Correspondence to: Manuel K. Schneider, Reckenholzstr. 191, 8046 Zurich, Switzerland, manuel.schneider@agroscope.admin.ch
#-------------------------------------------------------------------------------------------

# Libraries
library(nlme)
library(rjags)
library(coda)

#-------------------------------------------------------------------------------------------
# Read datasets
#-------------------------------------------------------------------------------------------
setwd("  ") #Set appropriate working directory
y <- as.matrix(read.table("SuppData_Differences_biomass_proportion.txt",
                           sep="\t", h=T, stringsAsFactors = F))
surveyInfo <- read.table("SuppData_Info_surveys.txt", sep="\t", h=T)
taxa <- read.table("SuppData_Species_traits.txt",
                    sep="\t", h=T, stringsAsFactors = F)

surveyInfo$Area <- factor(surveyInfo$Area)
surveyInfo$Survey <- 1:nrow(surveyInfo)

#-------------------------------------------------------------------------------------------
# Local model fitted using LME
#-------------------------------------------------------------------------------------------

# Choose trait
var <- "TRY_P"

# Transform data into long table
y2 <- data.frame(do.call("rbind", lapply(1:nrow(y), function(i) na.exclude(cbind(i, y[i,], taxa[,var])))))
colnames(y2) <- c("Survey", "y_unsc", "tr_unsc")
y2 <- cbind(y2, surveyInfo[as.numeric(y2$Survey),])
y2$tr <- scale(y2$tr_unsc)
y2$y <- scale(y2$y_unsc)

# Set up model

# Example model works with diagonal positive-definite variance-covariance matrices pdDiag.
# Where possible and computation time allows, use pdLogChol (general positive-definite variance-covariance matrices)
M1 <- lme(y ~ tr + tr*Breed + tr*Rotation+ tr*Area, random=list(Paddock=pdDiag(~tr), Subplot=pdDiag(~tr), Survey=pdDiag(~tr)), 
             corr=corCompSymm(form=~1) , contrasts = list(Rotation="contr.sum", Area="contr.sum"), data=y2, 
             control=list(maxIter=1e100, msMaxIter=1000000, msVerbose=F, niterEM=50, tolerance=1e-100, msTol = 1e-100, opt="optim"))

# A second model is used to test the effect of HC vs. OB
contrastHC <- list(Rotation="contr.sum", Area="contr.sum", Breed=contr.treatment(levels(surveyInfo$Breed), base=2))
M2 <- update(M1, contrasts = contrastHC, data=y2)

# Plot results
cols <- c("#1b9e77","#d95f02", "#7570b3")

par(mfrow=c(1,1))
ylm <- c(-1, 1)
plot(1,1, t="n", xlim=range(y2$tr, na.rm=T), ylim=ylm, axes=F, ylab="", xlab="")
xticks <- pretty(range(y2$tr, na.rm=T)*attr(y2$tr, "scaled:scale") + attr(y2$tr, "scaled:center"))
axis(1, at=(xticks - attr(y2$tr, "scaled:center"))/attr(y2$tr, "scaled:scale"), labels=xticks)
yticks <- pretty(ylm*attr(y2$y, "scaled:scale") + attr(y2$y, "scaled:center"))
axis(2, at=(yticks - attr(y2$y, "scaled:center"))/attr(y2$y, "scaled:scale"), labels=yticks, las=2)
box()
for(i in surveyInfo$Survey){
  br <- surveyInfo$Breed[i]
  sel <- y2$y[y2$Survey==i,]
  vals <- y2$tr[y2$Survey==i,]
  xrange<- range(vals, na.rm=T)
  xvals <- seq(xrange[1], xrange[2], length.out=10)
  #GLMM pediction
  lines(xvals, predict(M1, newdata=data.frame(surveyInfo[i,c("Survey", "Area", "Paddock", "Subplot","Rotation","Breed")], tr=xvals, row.names = NULL)),
        col=cols[br], lty=2)
} 

abline(a=fixef(M1)["(Intercept)"], b=fixef(M1)["tr"], col=cols[1], lwd=3)
abline(a=sum(fixef(M1)[c("(Intercept)", "BreedHC")]), b=sum(fixef(M1)[c("tr", "tr:BreedHC")]), col=cols[2], lwd=3)
abline(a=sum(fixef(M1)[c("(Intercept)", "BreedOB")]), b=sum(fixef(M1)[c("tr", "tr:BreedOB")]), col=cols[3], lwd=3)
legend("bottomleft", levels(surveyInfo$Breed),col=cols, lwd=3, h=T, xpd=NA, cex=1, bty="n")

p.values <- anova(M1)[c("tr","tr:Breed","tr:Rotation","tr:Area"),  "p-value"]
symb <- as.character(cut(p.values, breaks=c(0,0.001,0.01,0.05,0.1,1), labels=c("***","**","*","°","ns")))
legend("topleft", col="white", adj=0, legend=c("trait", "trait x breed", "trait x rotation", "trait x area"),
       bty="n", x.intersp=.9, inset=c(0.05,0))
legend("topleft", col="white", adj=0, legend=symb, bty="n", x.intersp=.9, inset=c(.25,0))
p.values <- c(summary(M1)$tTable[c("tr:BreedOB", "tr:BreedHC"),"p-value"], summary(M2)$tTable[c("tr:BreedOB"),"p-value"])
legend("topleft", col="white", adj=0, legend=c("AxH vs. OB", "AxH vs. HC", "OB vs. HC"),
       bty="n", x.intersp=.9, inset=c(.5,0))
symb <- as.character(cut(p.values, breaks=c(0,0.001,0.01,0.05,0.1,1), labels=c("***","**","*","°","ns")))
legend("topleft", col="white", adj=0, legend=symb[1:3], bty="n", x.intersp=.9, inset=c(.7,0))

#-------------------------------------------------------------------------------------------------------
# Fitting the global model
#-------------------------------------------------------------------------------------------------------

# Choose trait
var <- "TRY_P"

# Standardize trait values
  trait <- taxa[,var]
  tr.sd  <- sd(trait, na.rm=T)
  tr.mean <- mean(trait, na.rm=T)
  tr.min <- min(trait, na.rm=T)
  tr.max  <- max(trait, na.rm=T)
  tr <- (trait-tr.mean)/tr.sd
  tr.min.stand  <- min(tr, na.rm=T) #Standardized minima used for sampling missing trait values
  tr.max.stand  <- max(tr, na.rm=T)


  jags.data <- list(N=nrow(y),
                    J =ncol(y),
                    breed = samp$Breed,
                    num.lv = 2, #number of latent variables
                    y=y,
                    tr=tr,
                    tr.min.stand=tr.min,
                    tr.max.stand=tr.max,
                    v=as.numeric(samp$Subplot),
                    V=length(levels(samp$Subplot)),
                    p=as.numeric(samp$Paddock),
                    P=length(levels(samp$Paddock)),
                    a=as.numeric(samp$Area),
                    A=3,
                    r=as.numeric(samp$Rotation),
                    R=3
  )
  
  # Model building
  J.jags <- "
model{

# Level 1: modeling observations i
for(j in 1:J) {
   for (i in 1:N){
      y[i, j] ~ dnorm(mu[i, j], tau) # Distribution for observations # oder y[i, j] ~ dnorm(mu[i, j], tau * w[i])
      mu[i, j] <- beta[j, breed[i]] +   # No intercept here because of categorical x which automatically generates intercept
                inprod(lam.v[j, 1:num.lv], z.v[v[i],]) +   #latent variable for subplot (66)
                inprod(lam.p[j, 1:num.lv], z.p[p[i],]) +   #latent variable for paddock (9)
                inprod(lam.a[j, 1:num.lv], z.a[a[i],])  +   #latent variable for area (3)
                inprod(lam.r[j, 1:num.lv], z.r[r[i],])     #latent variable for rotation (3)
   } #i (=observations)
  } #j (=species)

tau <- 1/ (sd * sd)
sd ~ dunif(0, 100)

# Levels 2: modeling beta paramater
 for (k in 1:3){
    for (j in 1:J){
      beta[j, k] ~ dnorm(mu.beta[j,k], tau.beta[k]) # Distribution for beta per species and breed
      mu.beta[j, k] <- alpha0[k] + alpha1[k] * tr[j]       # Linear predictor for beta
    }
  alpha0[k] ~ dnorm(0, .0001)
  alpha1[k] ~ dnorm(0, .0001)
  tau.beta[k] <- 1/(sd.beta[k]*sd.beta[k])
  sd.beta[k] ~ dunif(0, 100)
 }

#Latent variable approach. Definition of dimensions of latent variables and application of constraints
#='Prior for latent variables'
# Vegetation subplot V
for(i in 1:V){ for(lv in 1:num.lv){ z.v[i,lv] ~dnorm(0,1) } } ## Prior for the latent variable. Build matrix for z.v and fill with random values
for(i in 1:(num.lv-1)) { for(lv in (i+1):num.lv) { lam.v[i,lv] <- 0 } } ## Constraints to 0 on upper diagonal
for(i in 1:num.lv) { lam.v[i,i] ~ dunif(0,20) } ## Sign constraints on diagonal elements
for(i in 2:num.lv) { for(lv in 1:(i-1)) { lam.v[i,lv] ~ dnorm(0,0.05) } } ## Free lower diagonals
for(i in (num.lv+1):J) { for(lv in 1:(num.lv)) { lam.v[i,lv] ~ dnorm(0,0.05) } } ## All other elements
# Paddock P
for(i in 1:P){ for(lv in 1:num.lv){ z.p[i,lv] ~dnorm(0,1) } } 
for(i in 1:(num.lv-1)) { for(lv in (i+1):num.lv) { lam.p[i,lv] <- 0 } } 
for(i in 1:num.lv) { lam.p[i,i] ~ dunif(0,20) } 
for(i in 2:num.lv) { for(lv in 1:(i-1)) { lam.p[i,lv] ~ dnorm(0,0.05) } } 
for(i in (num.lv+1):J) { for(lv in 1:(num.lv)) { lam.p[i,lv] ~ dnorm(0,0.05) } } 
# Area A
for(i in 1:A){ for(lv in 1:num.lv){ z.a[i,lv] ~dnorm(0,1) } } 
for(i in 1:(num.lv-1)) { for(lv in (i+1):num.lv) { lam.a[i,lv] <- 0 } } 
for(i in 1:num.lv) { lam.a[i,i] ~ dunif(0,20) } 
for(i in 2:num.lv) { for(lv in 1:(i-1)) { lam.a[i,lv] ~ dnorm(0,0.05) } } 
for(i in (num.lv+1):J) { for(lv in 1:(num.lv)) { lam.a[i,lv] ~ dnorm(0,0.05) } } 
# Rotation R
for(i in 1:R){ for(lv in 1:num.lv){ z.r[i,lv] ~dnorm(0,1) } } 
for(i in 1:(num.lv-1)) { for(lv in (i+1):num.lv) { lam.r[i,lv] <- 0 } } 
for(i in 1:num.lv) { lam.r[i,i] ~ dunif(0,20) } 
for(i in 2:num.lv) { for(lv in 1:(i-1)) { lam.r[i,lv] ~ dnorm(0,0.05) } } 
for(i in (num.lv+1):J) { for(lv in 1:(num.lv)) { lam.r[i,lv] ~ dnorm(0,0.05) } } 

# Distribution for unobserved traits
for(j in 1:J){
    tr[j] ~ dunif(tr.min.stand,tr.max.stand) #we always take uniform values between min and max (on the standardized scale)
  }
}" #End model

# Initial values (not required for all)
inits <- function() list(sd= runif(1, 0, 10))

# MCMC settings
nb <- 2000	# no. of burn-in iterations
ns <- 10000	# no. of sampling iterations
nt <- 10  # thinning
nc <- 1	    # no. of chains (only 1, because sign could change)

# Parameters monitored
params <- c("alpha0", "alpha1", "beta", "sd", "sd.beta", "tr")

# Set-up sampler
Jag <- jags.model(textConnection(J.jags), data = jags.data, n.chains = nc,
                  n.adapt=0, inits=inits)

# Adaptative phase to maximize the efficiency of the MCMC samplers (this is NOT a burnin)
adapt(Jag, 1000, end.adaptation=T)	# adaptation period (not a Markov chain)

# Burnin
update(Jag, nb)				# nb = burn-in period

# Generate posterior samples
samples <- coda.samples(Jag, params,  n.iter=ns, thin=nt) # generate samples from the posterior

#Plotting
#-----------------------------------------------------------------
cols <- c("#1b9e77","#d95f02", "#7570b3")

  alphas <- samples[,which(grepl("alpha",colnames(samples[[1]])))] #samples alpha parameters
  alpha.summary <- summary(alphas)[[2]] #2.part of summary -> quantiles of alpha parameters
  sd.summary <- summary(samples[,which(grepl("sd",colnames(samples[[1]])))])[[2]]
  tr.summary <- summary(samples[,which(grepl("tr",colnames(samples[[1]])))])[[2]]
  beta.summary <- summary(samples[,which(grepl("beta",colnames(samples[[1]])))])[[2]] 
  specmat <- matrix(beta.summary[1:(ncol(y)*3),3], ncol(y),3) #beta of all species for 3 breeds
  a_params <- alpha.summary[1:6,3] #"3" -> median
  Nobs <- apply(y, 2, function(x) tapply(x, surveyInfo$Breed, function(z) sum(!is.na(z))))
  
  xseq <- seq(tr.min, tr.max, length.out=20)
  xvals <- tr.summary[,3]*tr.sd+tr.mean
  
  plot(1, 1, t="n", xlim=c(min(xvals, na.rm=T), max(xvals, na.rm=T)), ylim=c(-.25,.25), xlab=var, mgp=c(1.5,.3,0), ylab="", axes=F)
  for(i in 1:3){
    bound <- apply(sapply(((xseq-tr.mean)/tr.sd), function(x) as.matrix(alphas[[1]])[,i] + as.matrix(alphas[[1]])[,i+3] * x),2,quantile, c(.025, .975))
    polygon(c(xseq, rev(xseq)), c(bound[1,], rev(bound[2,])), col=paste0(cols[i],50), border=NA)
    points(xvals, specmat[,i], pch=21, col=cols[i], bg=c("white",paste0(cols[i],50))[as.numeric(!is.na(tr))+1], cex=Nobs[i,]/15)
    }
  for(i in 1:3) lines(xseq, a_params[i]+ a_params[i+3]*(xseq/tr.sd-tr.mean), col=cols[i], lwd=3)

  axis(1, mgp=c(1,.3,0), tcl=-.3)
  axis(2, las=2, mgp=c(1,.3,0), tcl=-.3)
  box()

  legend("bottomleft", c("Angus x Holstein", "Highland Cattle", "Original Braunvieh")[c(1,3,2)],col=cols[c(1,3,2)],
         yjust=.5, pt.bg=paste0(cols,50)[c(1,3,2)], pch=21, lwd=2, h=T, bty="n")

  # Calculate probabilities for alpha values
  # Note: probabilities stabilize with the number of samples and maybe variable at small sample numbers
  getStats <- function(a){
  am <- as.matrix(a[[1]])
  comp <- cbind(am[,4]- am[,5], am[,4]- am[,6], am[,6]- am[,5])
  out <- c(ptrait = mean(rowMeans(am[,4:6])>0),
           pbreed = mean(rowSums(comp)>0),
           AxH_vs_OB = mean(am[,4]- am[,6]>0),
           AxH_vs_HR = mean(am[,4]- am[,5]>0),
           OB_vs_HR = mean(am[,6]- am[,5]>0))
  out[which(out>0.5)] <- 1-out[which(out>0.5)]
  return(out)
}
  symb <- as.character(cut(getStats(alphas), breaks=c(-1,0.001,0.01,0.05,0.1,1),
                           labels=c("***","**","*","°","ns")))
  legend("topleft", col="white", adj=0, legend=c("trait", "trait x breed"),
         bty="n", x.intersp=.9, inset=c(0.05,0))
  legend("topleft", col="white", adj=0, legend=symb[1:2], bty="n", x.intersp=.9, inset=c(.3,0))
  legend("topleft", col="white", adj=0, legend=c("AxH vs. OB", "AxH vs. HC", "OB vs. HC"),
         bty="n", x.intersp=.9, inset=c(.5,0))
  legend("topleft", col="white", adj=0, legend=symb[3:5], bty="n", x.intersp=.9, inset=c(.7,0))
  


