#' Plots a regression nomogram showing covariate distribution.
#' @description 
#' \code{regplot} plots a regression nomogram. Covariate distributions are superimposed on nomogram scales and the plot
#'	is animated to allow on-the-fly changes to distribution representation and to 
#'	enable outcome calculation. 
#' @param reg A regression object of either class \code{glm}, \code{glm.nb}, \code{lm}, \code{survreg} or \code{coxph}
#'  
#' @param dummies \code{TRUE} to treat dummy indicators of factor variables as distinct binary variables, with their own nomogram panel. 
#' Otherwise different categories of the factor are represented on the same panel. 
#' @param points If \code{FALSE} the regression scores \eqn{\beta}\eqn{x}  are shown. 
#' Otherwise the scale is represented by a 0-100 "points"  scale.  
#' @param center  Produces plot in which regression score contributions of continuous data are  plotted with respect 
#' to mean values of non-factors.
#' @param observation  An observation, data frame,  whose values
#'  are superimposed on the plot. Must include all variables used
#'   in the regression formula. 
#' @param title A heading title written to the plot
#' 
#' @param failtime   Specifies the cut-off time(s) for
#' plotting the risk  nomogram of a \code{coxph} or \code{survreg} regression 
#' (if \code{failtime=NULL}, cut-off is the median of time variable)
#' @param prfail For survival models. \code{TRUE} if 
#' probability of failure before \code{failtime} is summarised, otherwise after \code{failtime}.
#'  
#' @param baseS  Only relevant for \code{coxph} regressions and, if non-null, specifies the baseline 
#' survival probability, or probabilities, corresponding to the specified \code{failtime} 
#' (otherwise this probability is computed within \code{regplot}using function \code{basehaz}).
#' @param droplines Draws faint vertical lines showing score contributions to an observation.
#' @param nsamp  The size of a random sample of data to plot distributions 
#' (if  huge data plotting may be slow). 
#' @param showi Whether interactions (if present) are to be shown as a panel of the plot.
#' @param showP Whether P-value regression coefficient asterisk codes are to be displayed.
#' @param odds For probability outcomes (logistic, cox and survreg models), 
#' the  output nomogram scale is odds rather than probability.
#' @param rank Allows  panels to be ranked by importance to the regression. Importance is
#' measured by standard deviation along nomogram scales 
#' (equivalent to standardized regression coefficients for non-factors) 
#' \code{NULL} for no ranking (ordered as the regression), otherwise \code{"decreasing"}
#'  or \code{"increasing"}.
#' @param subticks Whether to put intervening minor tick marks on axes. 
#' @param interval  draws an interval (confidence or prediction) for outcomes associated with the 
#' \code{observation} (see Details)
#' @param colors A \code{list} of colors that will override the default colors. May include: \code{dencol} color fill of density plots, 
#' \code{boxcol} color fill of frequency boxes, \code{obscol} color of superimposed observation,
#' \code{spkcol} color of spikes.
#'  
#' @author Roger Marshall <rj.marshall@@auckland.ac.nz> The University of Auckland, New Zealand
#' 
#' @return If \code{points=TRUE} a \code{points.tables} object is returned 
#' as a list of dataframes, each frame giving  points per covariate, and a
#' total points-to-outcome look-up table.
#' 
#' @details Creates a nomogram  representation of a fitted regression. 
#' The distribution of covariates in the model, and of the total regression score, can be superimposed on the 
#' nomogram scales. Also the values of a particular observation can be overlaid, with  outcome calculated. 
#' If \code{interval} is specified a (\code{95\%}) confidence interval on the outcome is displayed. 
#' For \code{lm} and \code{glm} OLS regressions a prediction interval can be requested.
#' The plot is active for mouse input 
#' allowing on-the-fly changes to distribution type (boxes, spikes, box plot, density, empirical cdf, violin and bean plots)
#' and also to observation values, making it a regression calculator. 
#' The regression object \code{reg} parameter must have been fitted by either \code{glm},  \code{lm},  
#' \code{survreg}, \code{coxph}  or \code{glm.nb}. For \code{glm},
#' the supported family/link pairings are: gaussian/identity, binomial/logit, and  poisson/log. 
#' For \code{survreg} the distribution may be lognormal, gaussian, weibull, exponential or loglogistic.
#' For \code{glm.nb} (negative binomial regression) only log-link is allowed. 
#' @examples
#' library(survival)
#'  data(pbc) 
#'  pbc$catbili <- cut(pbc$bili,breaks=c(-Inf, 2, 4, Inf),
#'                 labels=c("low","medium","high"))
#'  pbc$died <- pbc$status==2

#'  pbccox <-  coxph(formula = Surv(time,died) ~  age + catbili + sex + 
#'                   copper +stage + trt,data=pbc)
#'  ## Plot a Cox survival model, showing data for the first observation.
#'  ## Display risk for 730, and 1825 days
#'   regplot(pbccox,observation=pbc[1,],  failtime = c(730,1825), prfail = TRUE ) 
#'  ## Plot a Weibull model 
#'   pbcweib <-  survreg(formula = Surv(time,died) ~  age + catbili + sex + 
#'                   copper +stage + trt,dist="weibull", data=pbc)
#'  regplot(pbcweib,observation=pbc[1,], failtime = 1825, prfail = TRUE ) 
#'  ## Plot a logistic regression, showing odds scale and confidence interval
#'  pbcglm <- glm(died ~  age + catbili + sex + copper,family = "binomial", data=pbc )
#'  regplot(pbcglm, observation=pbc[1,], odds=TRUE,interval="confidence")

regplot <- function(reg,  dummies=FALSE,
            center=TRUE,observation=NULL,title=NULL,points=TRUE,
            failtime=NULL, prfail=TRUE,baseS=NULL,odds=FALSE,
            droplines=FALSE,nsamp=5000,showi=TRUE,showP=TRUE,
            rank=NULL, subticks=FALSE, interval=NULL,
            colors=NULL)
  {
  
  
  if(is.null(reg)) {return(message("Error: NULL regression object"))}
  
  FIRSTRUN <- TRUE
  ANIMATE  <- TRUE
  
  
  ## distribution plots: default to none 
  cplot="none"
  dplot="none"
  plot_select <- 2 

  cannotclick <- FALSE
  new_obs <- TRUE

  ## is a RED (or obscol) "person" to be added?
  person <- !is.null(observation)
  if(person){
    if(nrow(observation)>1){
      message(paste("\"observation\" has >1 row. The last row provides plotted values"))
      observation <- observation[nrow(observation),]
    }

 if(!all(!is.na(observation))) {
 return(message("NA values in \"observation\"  not allowed  "))  
 } 
    
  }
  
  
  
  adddata <- NULL
  weighted <- FALSE 
 
  nudist <- TRUE
  
  pointscale <- points
  points_values <- list()

  rtype   <- class(reg)[1]
  cox     <- (rtype=="coxph") 
  lm      <- (rtype=="lm" )
  glm     <- (rtype=="glm")
  survreg <- (rtype=="survreg")
  glm.nb  <- (rtype=="negbin")
    
  if(!( glm | lm | cox | survreg | glm.nb)){
    return(message(paste("Cannot do  regplot for class",
    paste("\"",rtype,"\"",sep=""), 
      "models. Only \"lm\", \"glm\", \"coxph\", \"survreg\" \"negbin\" supported")) )
  }
  
##if(survreg | cox){library(survival)}
  
  if(survreg){
    dist <- reg$dist
    ## check on survreg() distribution types
    if(!(dist=="lognormal"   | dist=="loglogistic" | dist=="weibull"|
         dist=="exponential" |  dist=="gaussian" ) ){
    return(message(paste("Cannot do for survreg distribution:",dist)))
    }
  }
  
 
  if(glm | glm.nb| lm ){
    
  ## check on weights model  for glm, lm,
  lw <- ncol(reg$model)
  if(colnames(reg$model)[lw]=="(weights)"){
    message("Replicate weights assumed")
    weighted <- TRUE 
    W <- floor(reg$model[,lw])
    }
## extract Y variable from formula 

  yvar <- names(reg$model)[1]
  
  }
  
   if( cox | survreg ){
     
     if(cox){
 ## check on  whether Surv() is non-interval
       ## (note survreg unnecessary as it doesn't support
       ##  stop-start Surv objects)
     
     Yatt <- attributes(reg$y)
     scurvy <- unlist(Yatt[[2]][2])

     ## this should be either ""time" "status" 
     ## or "start" "stop" "status"
     if(!is.element("time",scurvy)){
       if(is.element("start",scurvy)) {
         return(message("start-stop time type Surv() objects are not supported in regplot"))
         }
         else
          { return(message("mispecified Surv() in formula"))}
     }
      if(Yatt$type != "right"){
       return(message("not a right censored Surv() object"))}
     
     }
     
    
    
  ## check on weights model for cox and survreg models
      
  if(!is.null(reg$weights)){
  message("Replicate weights assumed")
  weighted <- TRUE
  W <- floor(reg$weights)
   }
 
Sobject <- names(attr(reg$terms,"dataClasses")[1])
## extract time variable name  from  "Surv(time,fail)" object structure
##  messy, should be better way??
 stripped  <- sub(pattern="Surv\\(", replacement="", x=Sobject )
 stripped  <- sub(pattern="\\)",     replacement="",x=stripped)
 stripped  <- sub(pattern=" ",       replacement="",x=stripped)
 
 if(Sobject == stripped) {
## Sobject must be a Surv() obect specified outside of formula
   message("Surv() object ",paste(Sobject)," should preferably be specified in formula")
 }
 
 outcome_names <- strsplit(stripped,",")
 yvar <- unlist(outcome_names)[1]
 deadvar <-  unlist(outcome_names)[2]

 
 ##browser()
 ##if(is.na(deadvar)){message("Note: no censoring in model")}
 
 
 ## deadvar may be a string with values eg  "status == 1"
 ## pick out first word using this little function
 deadvar <- string_fun(deadvar)

     } ##if(cox |survreg)
  
coefficients <- reg$coefficients
   
  NAcoef <- which(is.na(coefficients))
  if(length(NAcoef)>0){
   

    return(message("model has NA coefficients. Reformulate without NA terms"))
 ## what happens if allow carry on??  
    ## make zero 
   ##coefficients[NAcoef] <- 0
  ##OKcoef <-  which(!is.na(coefficients)) 
  ##coefficients <- coefficients[OKcoef] 
  ##browser()
  }
  
  if(cox){
    ## There is possibility of exp(coef)=Inf from coxph
    if(!all(!is.infinite(exp(coefficients)))){
    return(message("coxph model with at least one exp(coef)=Inf. Probably a bad model"))
    }}
  ## patch in possibility of zero beta.  Throws factor variables!
  coefficients <- ifelse(coefficients==0,0.0001,coefficients)
 
  ## default graph title, if not specified
  if(is.null(title)) {title <- paste(rtype,"regression")
  if(survreg){title <- paste(dist,rtype)}
  }

##  check aon colors  list of parameters 

  if(!is.null(colors)){
  if(!all(is.element(names(colors),c("obscol","dencol","spkcol","boxcol") )) ){
    message("note: mispecified or incomplete \"colors\" list")
  }
  }
  ## v0.0 old colours dencol=#DDEEDD, boxcol=#DDDDEE spkcol=#444444
  if( is.null(colors$dencol) ) {dencol <- "#FFFF99"}  else {dencol<- colors$dencol}
  if( is.null(colors$boxcol) )  {boxcol <- "#A2D9CE"}    else {boxcol<- colors$boxcol}
  if( is.null(colors$obscol) ) {obscol <- "#FF0000"}      else {obscol<- colors$obscol}
  if( is.null(colors$spkcol) ) {spkcol <- "#669999"}     else {spkcol<- colors$spkcol}
  if( is.null(prfail) ) {fail <- TRUE}         else {fail <- prfail}
  if( is.null(failtime) ) {tcut <- NA}         else {tcut <- failtime}
  ntimes <- 1
  if(!is.na(tcut[1])) {ntimes <- length(tcut) }

   s5 <- vector(length=ntimes)
  if( is.null(droplines)) {droplines <- FALSE} else {droplines <- droplines}
  if( is.null(nsamp)) {nsamp <- 5000}          else {nsamp <- nsamp}
  if( is.null(showi)) {showi <- TRUE }         else {showi <- showi}
  if( is.null(baseS)) {baseS <- NULL }         else {baseS <- baseS}
  
  if(!is.null(baseS)){ nbaseS <- length(baseS)
    if(nbaseS != ntimes){
      return(message("non-NULL baseS and failtime must have same length"))}
  }
  if( is.null(showP)) {showP <- TRUE }         else {showP <- showP}
  if( is.null(odds))  {odds <- FALSE }         else {odds <- odds}
  if( is.null(rank))  {rank <- FALSE }         else {howrank <- rank
                                                 rank <- TRUE
                                                  if(howrank=="decreasing"){decrease <- TRUE}
                                                  if(howrank=="increasing"){decrease <- FALSE}
                                                  if(howrank!="decreasing" & howrank!="increasing"){
                                                  return(message("mis-specifed rank option"))
                                                  rank <- FALSE}
  }
     
 if( is.null(subticks)) {showsubticks <- FALSE } else {showsubticks <- subticks}
if( is.null(interval)) {interval <- "none"} 
   if(!all(is.element(interval,c("confidence","prediction","none") )) ){
    return(message("mispecified interval parameter"))
  }


 
  if( !is.null(reg$offset) ){
 ##   regplot not operational for offsets---------------------              
  return(message("regplot unable to make plots with an offset in regression"))
  }
  nai <- names(reg$coefficient[1])
  nointercept <- FALSE
  
  ## model with no intercept?  If Cox then has no intercept anyway.
  
  if(nai != "Intercept" & nai != "(Intercept)" & !cox ) {nointercept <- TRUE}

  ## interpet only models not done.
  if(nointercept==FALSE & length(reg$coefficient)==1){
   return(message("regplot pointless for intercept-only model"))
    
  }
  ##-----------------------------------------------------------------
  ## is formula ok? Checks: 
 variable_names <- attr(reg$terms,"term.labels")
 num <- length(variable_names)

  # ## possibility of  :: package prefix?. Strip away
 ##  why? :: operator interfers with : interaction 
 ## need to add prohibition because, clicking obs
 ## doesnt work.  reg$terms contains the package name prefix
  dblcolon <- grep("::", variable_names) 
 if(length(dblcolon>0)){
     return(message("package prefix :: detected in formula. Reformulate regression without"))
   }
 
 ## use getX function to test for arithmetic I() operators in variable names
  gX <- getX(variable_names)
 ## functions in last half of  returned gX
 fX <- gX[((length(gX)/2)+1):length(gX)]
 fX <- unique(fX)
 vX <- gX[1:(length(gX)/2)]

 functions_in_formula <- !all(fX=="")
 
 if(functions_in_formula){
 message("Note: function(s) specified within model formula: ",paste0( fX[which(fX!="")],collapse=" " ))

 ##toavoid <- c("poly","bs","cut")
 # isin <- is.element( toavoid, fX )
 # if(any(isin)) {
 #  
 #   message("Try avoid using: ", paste( toavoid[which(isin==TRUE)] ) )
 # }

if(interval != "none" &  person){
  ## need to shut down interval calcultaions if functions
  ## specified. Why? 
  ## Because predict() is used to calculate intervals and 
  ## newdata must be specified with raw observation data.
  ## e.g. if log(age) is used log(age) when clicked is obtained
  ## age itself must be updated
  ## at present not programmed to invert function. To update age.
  message("Function(s) in formula prohibit interval calculation")
  interval <- "none"

}


   if(is.element("cut",fX) & person){
   return( message(" \"observation\" prohibited with cut() in formula"))
     }
}  ##if(functions_in...
##------------------------------------------------------------ 
  
 
  
  
  ##===============================================================
  ##   Start  idling loop (ANIMATE) for  mouse clicks and Esc here
  ##===============================================================
  
   while(ANIMATE){
 
 if(FIRSTRUN){
   ## FIRSTRUN creates base plot and other computations on
   ## first cylcle of the idling
   actualdata <- model.frame(reg) }
 
     if(person){


if( glm |  lm | glm.nb ) { 

      
    if(FIRSTRUN){
        vactual <- colnames( actualdata )
        
        ## check whther variables (stripped of function) are in observation
        ## vX may have interactions. Remove these from vX 
  #     inter <-grep(":",vX) 
  #     vX <- setdiff(vX,vX[inter])
  # browser()
  #       inc <- is.element(vX,names(observation))
  #       if(!all(inc) ){
  #       message(paste("Variable(s)in glm/lm model but not in \"observation\":",sep=""))
  #      return(message(paste(vactual[which(inc==FALSE)],sep=",",collapse="  ")))
  #       }
      
       
       
      adddata <- model.frame(reg, data=observation)  
      
      vadded <-  colnames( adddata )
      included <- is.element( vactual, vadded )
      vin <- vactual[which( !included )]

      # 
      # if(length(vin)>0){
      # message(paste("Variable(s)in glm/lm model but not in \"observation\":",sep=""))
      #  return(message(paste(unique(vin),sep=",",collapse="  ")))
      # }
      
      
      
   }  ##if(FIRSTRUN)
  
 
  both <- rbind(actualdata,adddata)
  MM <- model.matrix(reg$terms,data=both)
  newX <- MM[nrow(MM),]

## strip away the intercept term
if(!nointercept){newX <- newX[-1]}

}  ##if(glm | glm.nb|  lm) 


if( cox  | survreg ){
## this is within an if(person)...  
  if(FIRSTRUN){
#  ## check on time in observation 
 
  if(!all(is.element( yvar, names(observation) ) ) ){
   return( message(paste(yvar),"  missing in \"observation\" "))
  }
    

   actualdata <- model.frame(reg)
  #remove first Surv(time, censor) element from vactual
   vactual <- colnames( actualdata )[-1]
   # 
   #  hascut <- grep("cut\\(",vactual)
   # if(length(hascut)>0){
   # return( message("non-NULL \"observation\" not allowed with cut() in formula"))
   #   }

  adddata <- model.frame(reg, data=observation)
  vadded <-  colnames( adddata )[-1]
  
   }  ##if(FIRSTRUN cox


  both <- rbind(actualdata,adddata)
  MM <- model.matrix(reg,data=both)
  ##  use last row of combined data for plotting red points 
  newX <- MM[nrow(MM),]

  if(colnames(MM)[1]=="Intercept" | colnames(MM)[1]=="(Intercept)" ){newX <- newX[-1]}
 
  
}  #endif(cox  | survreg)
 

  }  #endif(person)


##==================================================================  
if(FIRSTRUN){

  
##-------------------------------------------------------- 
## need to get out actual variable corresponding to dummy 
## variables of factors. Do this irrespective of whether dummies=TRUE
## -------------------------------------------------------- 

 variable_names <- attr(reg$terms,"term.labels")
 num <- length(variable_names)
 ##functions was here
 
 
  ols <- FALSE
  logistic <- FALSE
  poisson <- FALSE
  negbin <- FALSE
  if(glm | glm.nb){
    ## has glm done an allowed family/link?
    ols <-      (reg$family[1]=="gaussian" & reg$family[2]=="identity")
    logistic <- (reg$family[1]=="binomial" & reg$family[2]=="logit")
    poisson <-  (reg$family[1]=="poisson"  & reg$family[2]=="log")
    negbin <-  (glm.nb & reg$family[2]=="log")
    if(!(ols | logistic | poisson | negbin)){
      return(message(paste("Cannot do for glm family/link combination:",
            reg$family[1],reg$family[2])))
      }
      }
  
    
  if(!(cox | survreg | lm | ols | logistic | poisson |negbin )){
    return(message("regplot cannot do this regression model"))}
 
  if(!(ols|lm) & interval=="prediction"){
  message("Cannot do prediction interval for ",
    paste(rtype),".  Confidence done instead")
  interval <- "confidence"
}

  
  if(is.null( reg$formula )){ reg$formula <- formula(reg) }


nregcoef <- length(coefficients)
npars <- nregcoef

 ns <- ""

 if(glm | glm.nb| lm | cox | survreg){
   summry <- summary(reg)
   
 if(glm | glm.nb| lm) { Pval <- summry$coefficients[,4] }
#  P-values in different places, by model type
 # P values in 5th col of  coxph class
 if(cox)      { Pval <- summry$coefficients[,5] }
 if(survreg)  { Pval <- summry$table[,4] }  
 
 ns <- ifelse( Pval<.05   ,"*","" )
 ns <- ifelse( Pval<.01  ,"**",ns ) 
 ns <- ifelse( Pval<.001,"***",ns ) 
 #frig to patch if suppressing P values other option  showP=FALSE
 if(!showP){ns <- ifelse(Pval>0,"","")}
 }

  
  if(cox | survreg){
 ## need actualdata - may not have been done if !person
 if(!person) {actualdata <- model.frame(reg)}
   SU <- actualdata[,1]
 if(is.null(tcut[1]) | is.na(tcut[1]) ){
   ##get Surv object in SU extract first element which is time 
   ## and make time cutoff the median
   tcut <- median(SU[,1])
 
  }
  
  else
  {    if(max(tcut) > max(SU[,1]) | min(tcut) < min(SU[,1])){
     return( message("failtime ",paste(failtime)," out of range of observed: ",
       paste(signif(min(SU[,1]),4),signif(max(SU[,1]),4) ) ) )
    }
    
   }
  } 

  
  if(nointercept){message("Note: Regression model without an intercept")}
  
  if(cox| survreg ){ 
  
    if(cox){x <- model.matrix( reg,data=reg$model )}
    else {x <- model.matrix( reg )}
  

  }
  else  ##if(cox|survreg)
  { 
    x <- model.matrix( reg$formula, data=reg$model )
  }
  
 ## if replicate weighted model, expand out x matrix to unit records
  if(weighted){
   
  x <- x[rep(row.names(x),times =W), 1:ncol(x)]
  }
 ##---------------------------------------------------------------- 

  ##  large data?  Base distributions on subsample size nsamp (default 5000)
  ## if data huge, take subsample. Do so on first pass through


  if(nrow(x) > nsamp){
    nobs <- nrow(x)
   ## repeatable subsample  
   set.seed(1234)
    
    subs <- sample(1:nrow(x), size=nsamp)
    

    x  <- x[ subs, , drop = FALSE]
 
 
   message(paste("Distributions estimated using  nsamp=",nsamp,
  " random subsample of", nobs))}


 ## need to mess about with parameters. coxph class has no intercept.

 if(cox){ 

   betas <- coefficients[1:npars]
## fix up if npars=1, when x is a vector, coerce to matrix
   x <- as.matrix(x)
   X <- x[,1:npars, drop=FALSE]
   intercept <- 0
  }
 else   ##if(cox)
  {
   nai <- names(coefficients[1])
  nointercept <- FALSE
  if(nai !="Intercept" & nai !="(Intercept)" ){
 ## model with no intercept
    nointercept <- TRUE
    intercept <- 0
    betas <- coefficients[1:npars] 
    
    ##browser()
    X <- x[,1:npars,drop = FALSE] 
        }
        else
        {
      ## model with intercept
   betas <- coefficients[2:npars]
   ns <- ns[2:npars]
   intercept <- coefficients[1]
  X <- x[,2:npars,drop = FALSE]
 npars <- npars -1

        }

 }  ##if(cox




isinteraction <- vector(length=npars)
for (i in 1:(npars)){
   
  ## obtain last factor that is NOT an interaction (i.e. assume interactions all contain : and tagged last  
  if(length(grep(":", names( betas[i] ) ) ) == 0 ) {
  isinteraction[i] <- FALSE}
  else
  {
   
    isinteraction[i] <- TRUE}
  
}



## check for logical variables in the data and  coerce them into appearing
## as factors, by augmenting the  reg$xlevels list, which  gives levels of
## factor variables. 

## keep reg$xlevels and copy as "xlevels"
xlevels <- reg$xlevels
nvars <- length(variable_names)


  for(k in 1:nvars){
#  
  whichv <- which(names(actualdata)==variable_names[k])
  type_v <- class(actualdata[, whichv ] )
 
  
  if(type_v[1]=="logical"){
    
# message(paste("note: logical variable ",variable_names[k]," coerced to factor"))
# ### try coercing  to look like a factor   
     xxx <- c(xlevels,list(c("FALSE","TRUE")))
     names(xxx)[length(xxx)] <- variable_names[k]
     xlevels <- xxx
     
     

   }
#   
   }  ## end for(k in 1:nvars)


##   factorvarnames  does not include interactions at this point
factorvarnames <- names(xlevels)
n_interaction <- sum(isinteraction)

firstfactor <- NA
isfact <- FALSE

 ## which interactions are interactions of factors?
  ## try to fudge in xlevels for factor interactions.
  DF <- vector(length=nvars)
  isafactor <- vector(length=nvars)
  names(isafactor) <- variable_names
  names(DF) <- variable_names
  irow <- 1
## determine type of variables in the  raw data 
## get names of variables used from actualdata
  datavars <- names(actualdata)
  dataclass  <- lapply(actualdata,class)
  ## stop if special class variables found e.e. "poly"

  factorvar <- ifelse(dataclass=="factor" ,TRUE,FALSE)

  lvls <- vector(length=length(datavars))
  names(lvls) <- datavars
 ## determine variable levels (## catogories or polynomial)
 ## and  fix up logical varaibles to behave as factors
  for(v in 1:length(datavars)){
    if(factorvar[v]){lvls[v] <- nlevels(actualdata[,v])}
    if(dataclass[v]=="logical"){
      factorvar[v] <- TRUE
      lvls[v] <- 2}
  }

## fix up possibility of a character class. 
## make it like a factor. 
  
    for(v in 1:length(datavars)){
    if(dataclass[v]=="character"){
      factorvar[v] <- TRUE
      lvls[v] <- length(unique(actualdata[,v]))}
  }

  ## relies on ordering of effects - interaction last
  mixed_int <- FALSE
  is.mixedint <- logical(length=nvars)
  n_interaction_vars <- 0
   parameter_names <- names(betas) 
nmain_effects <- 0 
  for (k in 1:nvars){
    isfact <- FALSE
    nway <- length(grep(":", variable_names[k] ))
##  is an interaction term:
   
    if(nway>0){
     
      
     n_interaction_vars <- n_interaction_vars +1
     DF[k] <- pars_per_var(variable_names[k],parameter_names)
     ivars <- unlist(strsplit(variable_names[k],":"))
     
   ## probititions of poly() and bs() interactions. 
   ##  why?  Because  falls apart in some circumstances. 
     for(v in 1:length(ivars)){
     iv <- which(names(dataclass)==ivars[v])
     class <- unlist(dataclass[iv])[1]
     if(class=="poly"|class=="bs"){
     return(message("interactions that include class ",paste(class), "() prohibited"))}
     }
     
    
     
      isafactor[k] <- FALSE
      if(all(factorvar[ivars])){
   ##  all are factors in the interaction     
  ##establish how many parameters for this combo (degrees of freedom)
   ## it is a factor-interaction of order "nv"
        ## unclear whether this bit of code necessary
        ## because labels of factor interactions 
        ## obtained  from the paramater name. Anyway, retain as is. 
    EE <- expand.grid(xlevels[ivars], stringsAsFactors=FALSE)
   
    cross <- do.call(paste0, EE)
   
    xxx <- list(xxx=cross)
    names(xxx) <- variable_names[k]
    xlevels <- c(xlevels,xxx)
    isfact <- TRUE
    isafactor[k] <- TRUE
   ## browser()
      }
      else  ##if(all(isafactor
      {
        ##  numeric*factor interactions 
       ## prohibition of bs() and poly() interactions
        classes <- unlist(dataclass[ivars])
        if(is.element("poly",classes) & is.element("bs",classes) )
         {return(message("poly() with bs() interaction not allowed"))
          }
        if(DF[k] > 1){mixed_int <- TRUE
      is.mixedint[k] <- TRUE}
      }
    }
    else  ##if(nway >0
    { 
## not an interaction   
  nmain_effects <- nmain_effects + 1 
      if(factorvar[variable_names[k]]){

       DF[k] <- lvls[variable_names[k]] -1 
       isfact <- TRUE
      isafactor[k] <- TRUE
      }
      else
      {DF[k] <- 1 
      isafactor[k] <- FALSE
      }
     }   
    ## fudge in the additional parameter inserted with nointercept
    ## model when there is a factor
     if(nointercept & isfact & is.na(firstfactor) ) {

       firstfactor <- k
       DF[k] <- DF[k] +1}


} ## end for(k in 1:nvars)
##-----------------------------------------------------
##  check on whether main effects for interactions
## has bearing on Cox model: basehaz() which won't
## work unless all main effects included
if(n_interaction_vars>0){
nomain <- NULL
main_effects <- variable_names[1:nmain_effects]
for (k in (nmain_effects +1):nvars){
   
  ivars <- unlist(strsplit(variable_names[k],":"))
  for(v in 1:length(ivars)){
  if(!is.element(ivars[v],main_effects)) {nomain <- c(nomain," ",ivars[v])}
  }
}
if(!is.null(nomain)){
  nomain <- unique(nomain)
  message("Note: interaction but no main effects for:",paste(nomain,collapse=" "))
 if(cox & is.null(baseS)){
   message("=> may fail. Known coxph problem (see message)! Respecify model or supply baseS parameter")}
}

}
##-----------------------------------------------
## for mixed interactions (factor*numeric) in which the factor
##  has more than one level, need to 
## fudge variable_names, and update other arrays, after slotting
## in the names of the associated parameters

   if(mixed_int){
     new_variable_names <- NULL
     new_DF <- NULL
     new_isafactor <- NULL

    i <- 0
    for(k in 1:nvars){

      if(is.mixedint[k]){
       new_variable_names <- 
      c(new_variable_names, parameter_names[(i+1):(i+DF[k])]) 
      n_interaction_vars <- n_interaction_vars +DF[k] -1
      new_DF <- c(new_DF,rep(1,times=DF[k]) )
      new_isafactor <- c(new_isafactor,rep(FALSE,times=DF[k]))
      }
      else{
        new_variable_names <- 
         c(new_variable_names, variable_names[k])
        new_DF <- c(new_DF,DF[k])
        new_isafactor <- c(new_isafactor,isafactor[k])
      }
      
      i <- i +DF[k]
    }
    variable_names <- new_variable_names
    nvars <- length(variable_names)
    isafactor <- new_isafactor
    DF <- new_DF
  }    
  ## end of mixed terms   

#----------------------------------------------
 ## polynomial terms in model  
  
  if(is.element("poly",unlist(dataclass))){
   
    if(!dummies) return(message("\"poly\" class covariates detected. Needs dummies=TRUE"))
    if(person){
      message("\"poly\" class covariates detected. observation not actived")

    cannotclick <- TRUE
    }
  
  is.poly <- logical(length=nvars)
  is.poly <- rep(FALSE,times=nvars)
   for(v in 1:length(datavars)){

    if(unlist(dataclass[v])[1]=="poly"){
      kv <- which(variable_names==datavars[v])
      is.poly[kv] <- TRUE  }
  }

  expanded <- expand_vars(is.poly, variable_names,parameter_names,DF,isafactor,nvars)

  
  variable_names <- expanded$variable_names

  nvars <- expanded$nvars
  DF <- expanded$DF
  isafactor <- expanded$isafactor

  }
 ## end poly() 

 ##------------------------------------------------------ 
 ## spline bs() terms in model? Adopt same as for poly fudge  
  
  if(is.element("bs",unlist(dataclass))){
    

   
    if(!dummies) return(message("\"bs\" class covariates detected. Needs dummies=TRUE"))
    if(person){
      message("\"bs\" class covariates detected. observation not actived")

    cannotclick <- TRUE
    }
  
  is.bs <- logical(length=nvars)
  is.bs <- rep(FALSE,times=nvars)
   for(v in 1:length(datavars)){

    if(unlist(dataclass[v])[1]=="bs"){
      kv <- which(variable_names==datavars[v])
      is.bs[kv] <- TRUE  }
  }

  expanded <- expand_vars(is.bs, variable_names,parameter_names,DF,isafactor,nvars)

  
  variable_names <- expanded$variable_names

  nvars <- expanded$nvars
  DF <- expanded$DF
  isafactor <- expanded$isafactor

  }

## end bs() section-----------------------------
##   put in sum check
   if(sum(DF) != npars){
     message("Warning: checksum fail sum(DF) ",paste(sum(DF))," !=npars ",paste(npars),
       "Contact developer")
     }

if(!dummies){
  n_interaction_rows <- n_interaction_vars
}
else
{n_interaction_rows <- n_interaction }
  

##factorvarnames <- names(xlevels)
## number of main effect parameters
n_main_pars  <- npars - n_interaction

n_factor_rows <- sum(isafactor)
n_non_factor_rows <- sum(DF[which(!isafactor)])

## set up the rows

  
 
n_non_factors <- npars - sum(DF,na.rm=TRUE)

## make sure factors main effcts (they might not be, could be just interaction)
index_factornames <- which(is.element(factorvarnames,variable_names))

levels <-        xlevels[ index_factornames]
## number of factors
n_main_factors <- length(levels)
## if dummies requiring sum(xlevels) - n_factors  (ie. less on nlev per  dummy) rows 
# 
# n_main_factor_pars <- length(unlist(levels)) - n_main_factors
# ## ie less on category for each
# ## number of on-factor main effects
# n_nonfactor_main_pars <- n_main_pars - n_main_factor_pars 

if(n_main_pars ==0){
  message("This model has no main effects")
}
n_factors <- length(xlevels)
if(dummies){
 nrows <- npars
 row_names <- names(betas)
 
}
else
{
nrows <- nvars
row_names <- variable_names
}


    if(FIRSTRUN & n_interaction >0 & !showi){
    message("Model has interactions. Option showi=FALSE suppressed")
    
  }

 
vact <- vector(length = npars  )
isadummy <- vector(length = npars )
vact_names <- vector(length = npars  )


vact <- seq(1:npars)
isadummy <- rep(TRUE,times=npars)

i <- 0


k <- 0


for(k in 1:nvars ){
  
 
    nlev <- DF[k]
 ## need to create vact() giving variable number of each parameter (i)
  for(j in 1:nlev){
    i <- i + 1
    vact[i] <- k
  
    #if(isafactor[k] |  DF[k]>1 ) {
     if(isafactor[k]) {
      
      isadummy[i] <- TRUE
      vact_names[i] <- variable_names[k]
    }
    else
    {isadummy[i] <- FALSE
     vact_names[i] <- variable_names[k]
     
    }
  }
}  

  S <- X
  
 
## -------------------------------------------------------- 
## extract scores in a loop 
## -------------------------------------------------------- 


  for (i in 1:npars){
  
 
    S[,i] <- X[,i]*betas[i]

  }


##------------------------------------------------------------------------------------
## center the  scores and data if center

  
  if(center){m <- colSums(X)/nrow(X)}
  else
 ## default align  to minimum values?
    ## do I want this?? Later. No  align to zero  so that actual
    ## beta*X scores or shown 
  {  #m <- apply(X,2,min)
    m <- rep(0,times=ncol(X))
    }
  ## exclude factor variables from centering, make 0 always the reference
     m[which(isadummy)] <- 0
  
  X <- scale(X,center=m,scale=FALSE)
  bm <- m*betas
  S <- scale(S,center=bm,scale=FALSE)
  intercept <- intercept +sum(bm, na.rm=TRUE)
  
## if dummies=FALSE, need to make rows of plot !dummies rather than dummy variables
## of factors. Do this by integrating over the dummy variables of factors.
if(!dummies) {


SS <- matrix(nrow=nrow(X),ncol=nrows)
XX <- matrix(nrow=nrow(X),ncol=nrows)
# SS <- matrix(nrow=nrow(X),ncol=npars)
# XX <- matrix(nrow=nrow(X),ncol=npars)
  
i <- 1
k <- 1

while (i <= npars) {
  
## if(!isinteraction[i]){
  
  dummys <- which(vact==k)
##browser()
  
  if(length(dummys) >1) {
   SS[,k] <- rowSums( S[,dummys], na.rm=TRUE)
   XX[,k] <- SS[,k]
  
  
  i <- i + length(dummys)
  k <- k+1 
  }
  else
  {
  SS[,k] <-  S[,i] 
  XX[,k] <-  X[,i]
  
  i <- i+1
  k <- k+1  
  }
  
 
  
  
}  ##while (i <= npars) 


colnames(SS) <- row_names

## replace  S and X by collapsed values

S <- SS
X <- XX


}  ##if(!dummies)
  else
  {nrows <- nregcoef 
  ##take off the intercept if necessary
    if( !nointercept & !cox  ) {nrows <- nregcoef - 1}
  }
  ## number panels show. Normall =nrows, except when interactions 
  ## are supressed with showi=FALSE
  npanels <- nrows
  if(n_interaction>0 & !showi) {npanels <- nrows - n_interaction_rows}

##------------------------------------------------------------
## set graphic boundaries. Max and min of beta-scores or fraction thereof??
L <- min(S[,],na.rm=TRUE)
M <- max(S[,],na.rm=TRUE)
diff <-  (M-L)
tickl <- (M-L)/50
tot_score <-  rowSums(S, na.rm=TRUE) 

}  ##@if(FIRSTRUN)
##==================================================================
 
 
make_space <- 2 +(ntimes-1)*.6
## graphics margin + 4 
par(mar=c(0,0,0,0) + 4,mai=c(0,0,0,0),bg = "#FFFFFF") 

if(FIRSTRUN) { plot.new()
  plot(c(L-0.3*(M-L),M+0.1*(M-L)), c(0.5,npanels+2  + make_space),col="white") 

##nrows <- npars
##  delta controls positioning of scale annotation.
##  one delta unit below  
delta <- 0.05 + (0.2-0.05)*(npanels)/10
## long and short tick length
ltick <- 0.35*delta
stick <- 0.2*delta

##factr <- vector(length=nrows)
isinteraction_row <-  vector(length=nrows)

}

else
{if(nudist){
  plot(c(L-0.3*(M-L),M+0.1*(M-L)), c(0.5,npanels+2 + make_space),col="white" ) }
  else
 {replayPlot(bareplot)}
} 


sumcheck <- 0
aspect <- 1.4*(M-L)/(npanels+2 + make_space -0.5) 

# #try scaling point axis 0,-100  (L,M)

lscale <- max(S[,],na.rm=TRUE)- min(S[,],na.rm=TRUE)

if(FIRSTRUN){
## set up ranking permutation for rank=TRUE
if(rank){
 Rng <- vector(length=nrows) 
  for(i in 1:nrows){
    
    # MinS <- min(S[,i])
    # MaxS <- max(S[,i])
    # Rng[i] <- MaxS- MinS
    ## ranking criterion:  sd of S. Equivalnent to standardised beta?
    ## or range? 
    Rng[i] <- sd(S[,i])
  }
  
  irank <- order(Rng,decreasing=decrease)
  

}
  else
  {irank <- seq(1:nrows)}
  
  row_of_panel <- vector(length=nrows)
  
  

}  ##if(FIRSTRUN)


#=================================================================
## main graphics loop over potential nrows of the graphic. 
##  Those shown are "panels"  
#=================================================================
if(nudist){
  
none <- (cplot=="none" & dplot=="none")
ipos <- 1
npanel <- 0
 for (row_num in 1:nrows){
   
   
   i <- irank[row_num]
   
   
   ## deal with stat significance, choose largest if !dummies
   ## note: max("","*","**","***") returns "***"!
   
   if(!dummies){
     iv <- which(vact==i)
     ss  <- max(ns[iv], na.rm=TRUE)
    isint <- isinteraction[iv[1]]
        }
   else
   {ss <- ns[i]
   isint <- isinteraction[i] }
   
   isinteraction_row[i] <- isint
   

   ## patching switch to omit showing interaction panels
    showpanel <- ( !isinteraction_row[i] | showi )

    
   size <- 1
   ##v <- colnames(S)[i]
   v <- row_names[i]
   
   if(showpanel) {npanel <- npanel + 1
   row_of_panel[npanel] <- row_num}
   ipos <- npanel+0.5 +make_space
  
if(showpanel){
  ##  add the variable/parameter to l.h.s  of plot
  if(isint){
   ##reduce font size for complexity of interaction (number of colons)??
   ncolons <- length(unlist (gregexpr(":",v) ))
   size <- 1- ncolons*0.1
   text(x=L-0.3*(M-L),y= ipos+0.4,paste(v,ss,sep=""),cex=size, font=3,adj=c(0,1),col="#2471A3")

   } 
   else  #if(isint) 
   {
   text(x=L-0.3*(M-L),y= ipos+0.4,paste(v,ss,sep=""),cex=size, adj=c(0,1),col="#2471A3")
   }
 
}  
  
 
  #    
     if(dummies) { 
     numeric <- !isafactor[vact[i]]}
     else
     {numeric <- !isafactor[i]}
#====================================================================
 ### numeric (non-factor variable) or interaction
 ##  make a  distribution type
 ##====================================================================

  
  if(numeric){
   
    
   uniq <- sort(unique(S[,i]))
   Xuniq <- sort(unique(X[,i]))
   luniq <- length(Xuniq)
   
  
  if(!dummies){
    B <- betas[which(vact==i)]
  }
  else
  {B <- betas[i]}
  
  
  ##browser()
    
   ## need to reverse direction for -ve beta
  if(B < 0) {Xuniq <- sort(unique(X[,i]),decreasing=TRUE)
    }
 

    
        
     if(!none & showpanel){
      bandwidth <- 0.01*(max(S[,i])-min(S[,i])) 
     if(cplot=="violin"){

   box <-  vioplot(S[,i], add=TRUE,at=ipos+0.3,horizontal=TRUE,wex=0.5, border="blue", 
             h=bandwidth, colMed="black", pchMed=21, col=dencol)
            }

   if(cplot=="bean" ){
 
  box <-  beanplot(S[,i], add=TRUE,at=ipos+0.3,horizontal=TRUE, maxwidth=0.5,
       what = c(0,1,0,1),  bw=bandwidth, log="", overallline="median",beanlinewd=0, 
      kernel="epanechnikov",  axes=FALSE, border="blue", col=dencol,method="overplot",
      ll=0.1 , frame.plot=FALSE,  pars = list(boxwex = 0.8, staplewex = 0.5, 
                                              outwex = 0.5,cex=0.5))  
  }
  
  if(cplot=="density") {
    
    dens <- density(S[,i], bw=bandwidth)
    maxden <- max(dens$y)
    dens$y <- 0.5* dens$y/maxden +ipos
    
    zero <- c(rep(ipos,times=length(dens$x)) ) 
    
    ## do fill-in with handy bit of code grabbed from somewhere! 
    polygon(c(dens$x, rev(dens$x)), c(dens$y, zero),
      col = dencol, border = NA)
    lines(dens$x,dens$y,col="blue")
    segments(min(dens$x),ipos,max(dens$x),ipos,col="black")
  }
       
       
        if(cplot=="cumul") {
 
  orderS <- sort(S[,i])
  norderS <- length(orderS)
  ecdy <- (1:norderS)/norderS
  yspace <- 0.65
  ecdy <- yspace* ecdy +ipos
  ##browser()
  ## position vertical scale at about 20th percentile 
  n20 <- floor(norderS*0.2)
  ypos <- orderS[n20]
 lines(orderS,ecdy,type="s", col="blue")

 segments(ypos,ipos,ypos,ipos+yspace)
 text(ypos - tickl,ipos+yspace, paste("1"),cex=0.6,adj=0)
 text(ypos - tickl,ipos+0.5*yspace, paste("0.5"),cex=0.6,adj=0)


  segments(ypos, ipos+yspace, ypos +tickl,ipos+yspace)
  segments(ypos, ipos+0.5*yspace, ypos +tickl,ipos+0.5*yspace)
  segments(ypos, ipos+yspace, orderS[norderS],ipos+yspace,col="gray")
  segments(ypos, ipos+yspace*0.5, orderS[norderS],ipos+yspace*0.5,col="gray")
  }

 
    if(cplot=="boxplot"){
 
    box <-  boxplot(S[,i], add=TRUE,at=ipos+0.3,horizontal=TRUE,  border="blue",col=dencol,
             pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.5,cex=0.5))
         }
   
  
  ## spikes, frequencies of unique values
  
  if(cplot=="spikes"){
    ## continuous frequency table - will be ordered by  increasing S 
    
    frq <- as.numeric(table(S[,i]))
    
    maxfreq <- max(frq)
    
      
    dy <- 0.6*frq/maxfreq 
    dx <- dy * aspect
    xleft <-  uniq - dx
    xright <- uniq + dx
    ybottom <- c(rep(ipos,times=luniq)) 
    ytop <- c(rep(ipos,times=luniq)) + dy
    
  ##spikes of discrete 
    for(nspike in 1:length(uniq)){
     
      segments( uniq[nspike], ybottom[nspike], uniq[nspike], ytop[nspike], lwd=3,col=spkcol)
    }
    ## also add an axis 
    segments( min(uniq),ipos,max(uniq),ipos) 
  }
         
         
  
     } ##if(!none 



## now add the variable scale ruler
##  extract beta value for this numeric variable. 
if(!dummies){

  Beta <- betas[which(vact==i)]
}
else
{
  Beta <- betas[i]
}



  if(!dummies){  
    
   ## browser()
    mean <- m[which(vact==i)] }
 
  else
  {
    mean <- m[i]
    }
  
     ## nmax parameter used in pretty()     
  nmax = round(min(10, 1+10*(max(S[,i]) -min(S[,i]))/lscale ),0)
 
  ticks     <- pretty(X[,i]+mean,n=nmax) 
  ticks_pos <- ticks*Beta -mean*Beta 
  

## add the scale for this variable panel
if(showpanel){
segments(min(ticks_pos),ipos,max(ticks_pos),ipos,col="black",)
segments(ticks_pos,ipos,ticks_pos,ipos - ltick,col="black",)
if(showsubticks) subticks(ticks_pos,ticks_pos,ipos,stick)
text(x=ticks_pos,y= ipos-delta,paste(ticks),cex=0.65,col="black",)
}
     
 
 if(pointscale){    
## points corresponding to  main ticks  (for output, no other use):
tick_points <- round(100*(ticks_pos -L)/(M-L))
## assume non-negative X values not allowed
positv <- which(ticks >=0)

  pointsframe <-  as.data.frame( cbind(ticks[positv], tick_points[positv]) )
 
 colnames(pointsframe) <- c(v,"Points") 
  points_values[[nrows+1-i]] <-  pointsframe 

 }  ##if(pointscale)
    
  
   }  ##if(numeric)

   else  ##if(numeric)
   { 
 ##================================================
 ## draw  frequencies as boxes or spikes    
##================================================     
 
 
## frequency table - will be ordered by  increasing S 
     
     
   uniq <- sort(unique(S[,i]))
   luniq <- length(uniq)
   
  
   
     frq <- as.numeric(table(S[,i]))
   
      tot <- sum(frq)  
   

    if(!dummies){
      
       par_nos <- which(vact_names==variable_names[i])      
       B  <- betas[par_nos]
    
      ## make the labels for the  categories 
       if(!isint){
         ##use levels of factor variables unless an interaction
       cats <- unlist(xlevels[variable_names[i]])
       }
       else
       {
         ## interaction term: re-make Xuniq as the (character)  labels 
         ## using  betacoefficient names as labels
         ##  and set a baseline "ref" label
        lpar_nos <- length(par_nos)
        cats <- vector(length=lpar_nos)
        cats[2:(lpar_nos+1)] <- names(B)
        cats[1] <- "ref"
        }
      
      ## augment  B with baseline value zero
      if(nointercept & firstfactor==i){
       
          beta_values <- B}
      
      else
       
       {
      beta_values <- c(0,B)
      ## artificial beta=0 for ref category 
      ## pick up baseline category name
      names(beta_values[1]) <- names(cats[1])
      }
       ##  Need to order with 0 added to bring in line with order of frq table
      obetas <- order(beta_values)
     
      cats <- cats[obetas]
      
      luniq <- length(cats)
      
      ##val <- Xuniq[o]
      
      val <- cats
       
  ##  but uniq also needs to be ordered by  magnitude    
      
      ofact <- order(uniq)
      uniq <- uniq[ofact]
      
      
    }  ##if(!dummies)
    
    else  ##of if(!dummies) 
      
    { 
      
      ## must be 0 - 1 categories, need to soet by sign of beta
      
        ## need to reverse direction for -ve beta
      if(betas[i] < 0) {cats <- sort(unique(X[,i]),decreasing=TRUE)}
      else
      {cats <- sort(unique(X[,i]),decreasing=FALSE)}

     
     val <- signif(cats+m[i], 2) }
              
  
o <- order(frq,decreasing=TRUE)
frq <- frq[o]
uniq  <- uniq[o]
val <- val[o]


 if(showpanel){
   
    if(dplot !="none" ){
      ## draw squares so that  largest is first, smaller superimposed
      if(dplot=="spikes"){
              dy <- 0.6*frq/tot 
              dx <- 0.001* diff
              xleft <-  uniq - dx
              xright <- uniq + dx
              ybottom <- c(rep(ipos,times=luniq)) 
              ytop <- c(rep(ipos,times=luniq)) + dy
             
      }
      else  ##if(spikes)
      {
        
   dy <- 0.35*sqrt(frq/tot) 
   dx <- dy * aspect
   xleft <-  uniq - dx
   xright <- uniq + dx
   ybottom <- c(rep(ipos,times=luniq)) - dy
   ytop <- c(rep(ipos,times=luniq)) + dy
   
   
      }
    }
   else   ##if(dplot !="none")
   {dx <- 0.01
    dy <- 0.01
    xleft <-  uniq - dx
    xright <- uniq + dx
    ybottom <- c(rep(ipos,times=luniq)) - dy
    ytop <- c(rep(ipos,times=luniq)) + dy
   }
   
    
   ## is possible that there is no such combination in interactions
   if(length(uniq) >0) {
   

## also add a faint axis ??
       segments( min(uniq),ipos,max(uniq),ipos,col="gray") 
   
   if(dplot!="none"){
    colfill <-  boxcol
    
    ## do I want to distinguish  actual factors from pseudo factors
    ## ie with fewer than ndisc   values?  How? 

     if(dplot=="spikes"){
       for(nspike in 1:length(uniq)){
        
       segments( uniq[nspike], ybottom[nspike], uniq[nspike], ytop[nspike], 
         lwd=3,col=spkcol)
       }
        ## also add an axis 
       segments( min(uniq),ipos,max(uniq),ipos) 
       
      ## put black labels at top of spike
 
       text(uniq,ipos+dy+0.15,paste(val),cex=0.7,col="black")
     }
     else
     {
       ##boxes
   
     rect(xleft, ybottom, xleft+2*dx, ytop,border="blue",col=colfill)
     points(0.5*(xleft+xright),(ybottom+ytop)*0.5,cex=0.3, col="black")
     ## put labels at bottom of box ( -dy)
     
      ##   no distribution "none"
     ## try avoid overlay, by alternating above below point?
     ## alternating 1,-1,1,-1 sequence
    plmi <- rep(c(1,-1), times=length(val))
    plmi <- plmi[1:length(val)]
    plmi <- plmi[o]
    posn <- ipos  - (dy+0.8*delta)*plmi
    
     
     
     text(uniq,posn,paste(val),cex=0.7,col="black")
     }
     
   }  ##if(!none) 
   
    
   else  ##if(!none) 
     
   {
  
   
 
     ## try avoid overlay, by alternating above below point?
     ## alternating 1,-1,1,-1 sequence
    plmi <- rep(c(1,-1), times=length(val))
    plmi <- plmi[1:length(val)]
    plmi <- plmi[o]
    posn <- ipos  - delta*plmi
 ## browser()
  if(is.null(val)) val <- " "
   text(uniq,posn,paste(val),cex=0.7,col="black")
   size <- 0.3
   
   ## draw vertical up/down ticks 
   ## also add an  axis 
    segments( min(uniq),ipos,max(uniq),ipos,col="black") 
    segments(uniq, ipos, uniq, ipos + 0.4*delta*plmi)
 
        }
   ##  put in point at center of box
 
 
 
}  ##if(length(uniq) >0)

}  ##showpanel
   
     if(pointscale){

        tick_points <- round(100*(uniq-L)/(M-L))
        pointsframe <-  as.data.frame( cbind(val,tick_points) )
        colnames(pointsframe) <- c(v,"Points") 
  ## make sure output ordered list same as graphic 
        points_values[[nrows + 1-i]] <-  pointsframe
               
     }
     
}  ##else  of if(numeric
}   ## end of main graphics loop for (row_num in 1:(nrows))

#-----------------------------------------------------
# add beta X  or points scale  contributions at top of graphic
#-----------------------------------------------------

tickp <- pretty(S[,],n=8)
# values of the ticks, keep in tickval

tickval <- tickp
# change output scale to points, 
if(pointscale) {
  
 
  tickval <- 100*(tickp -L)/(M-L)
  tickval <- pretty(tickval,n=8)
  tickval <- tickval[which(tickval >= 0 & tickval <=100 )]
  
  tickp <- ((M-L)*tickval)/100 +L 

}
npos <- npanels+1 + make_space
ypos <- npos+0.4
segments(min(tickp),ypos,max(tickp),ypos)
segments(tickp,ypos,tickp,ypos-ltick)
cexv <- 0.9
if(npars >15) {cexv <- 0.7}
text(x=tickp,y= ypos- delta,paste(signif(tickval,3)),cex=cexv) 

# subticks?

if(showsubticks) subticks(tickp, tickp,ypos,stick)



if(pointscale) {ex <- "Points  "}
  else
      {if(center){ex <- expression(italic(paste(beta,"(X-m) terms  ")))}
    else
      {ex <- expression(italic(paste(beta,"X terms  ")))}
  }

text(min(tickp), ypos , adj=1, ex ,cex=0.8, font=3) 

#write titleheading above  upper scale 
titlepos <- min(ypos + 2.5*delta,npars+2+ make_space)

text((M+L)/2,titlepos, title ,cex=1.05 ,font=4) 


# =======================================================
# add  Total points scale and its distribution
##use yref_pos as reference point of the caption "Total-points-to outcome"
#========================================================
 yref_pos <- make_space + 0.5 
allfact <- all(isafactor)
Total_distribution(L,M,pointscale, tot_score,none, 
  dplot,cplot,nrows,ltick, stick, center,delta,dencol,boxcol,
  spkcol,tickl,aspect,yref_pos,showsubticks,allfact,plot_select)
#==========================================================

npos <- 1
ypos <- npos+0.5
##total score axis is a position 1.5 (in Total_distribution)
totscore_ypos <- 1.5

Max_tot_score <- max(tot_score)
Min_tot_score <- min(tot_score)
Range <- Max_tot_score-Min_tot_score

  #if pointscale adjust scale to output pretty points


tickvalues <- pretty(tot_score,n=6)

if(!pointscale){
tickv <- tickvalues
tickp <- (M-L)* (tickv -Min_tot_score)/Range + L


## fix positions of the scale relative to values of pretty values
## draw line to fit L-M boundary, ticks having  screen coordinates
}
else
{

## require "pts" for pointscale
  tickv <- tickvalues
  
  ## convert to points, but need to account for there being sum of npars parameters
  tickv <- 100*(tickv -L*nrows)/(M-L)
  tickv <- pretty(tickv,n=8)
  ## convert back into actual total scores
  pts  <- ((M-L)*tickv)/100 +L*nrows
  tickp <- (M-L)* (pts -Min_tot_score)/Range + L 
  
  
} 

##==============================================================
## add score-to-probability nomograms at the bottom             
##==============================================================
#

pos_of_nomo <- yref_pos-1.8

# ------logistic  nomogram-----------------------------------------  
if(logistic){
    Logistic_scale(tot_score,intercept,L,M,Range,Min_tot_score,delta,
    ltick,stick,pos_of_nomo,yvar,odds,showsubticks)
 }

# ------lm  nomogram-----------------------------------------  
if(lm | ols ){
   Lm_scale(tot_score,intercept,L,M,Range,Min_tot_score,delta,
   ltick,stick,pos_of_nomo,yvar,showsubticks)
      }

# ------poisson  nomogram-----------------------------------------  
if(poisson | negbin){
   Poisson_scale(tot_score,intercept,L,M,Range,Min_tot_score,delta,
  ltick,stick,pos_of_nomo,yvar,negbin,showsubticks)
 }

# ------Cox nomogram-----------------------------------------  


for(time in 1:ntimes){
pos_of_nomo <- yref_pos - 1.8  -(time -1)*0.6 

  
  ftime <- tcut[time]
  base <- baseS[time]
if(cox){
  ## use basehaz() to get hazards at X=0. Only do FIRSTRUN as it
  ## is time consuming

  if(FIRSTRUN & is.null(baseS) & time==1) {
  h <- basehaz(reg, centered=FALSE)
  
  }
  
  cxp <- Cox_scale(tot_score,intercept,L,M,Range,Min_tot_score,delta,
  ltick,stick,pos_of_nomo,yvar,odds,base,h,s5,betas,center,m,bm,ftime,
  fail,showsubticks)
  s5[time]  <- cxp$s5
 
}

## survreg nomogram -------------------------------------------
if(survreg){
  p <- 1/reg$scale
  Survreg_scale(tot_score,intercept,L,M,Range,Min_tot_score,delta,
  ltick,stick,pos_of_nomo,yvar,betas,m,ftime,fail,p,dist,odds,showsubticks) 
}
}
##================================================
## finalise last total points table 
 if(FIRSTRUN & pointscale){
# make table for first [1]  element of s5 and tcut
    points_values[[nrows + 1]] <- 
    points_tables(survreg,dist,p,cox,s5[1],fail,tcut[1],poisson,
      negbin,lm,ols,logistic,pts,intercept,npars,yvar,tickv)  
        }
##======================================================
##   ANIMATE <- !is.null(observation)
}  #if(nudist)


   if(ANIMATE){
    nudist <- FALSE 
    
    ##save plot without any red observation data
     bareplot <- recordPlot()
     
     
  ##-========================================================
  ## start overlay  red (obscol) observation dots, text  and lines   
  #==========================================================
     if(person){
     npanel <- 0
     for (row_num in 1:(nrows)){
       i <- irank[row_num]
     

#  add in  RED points  of observation and value 
    
    if(dummies){
      
        addX <- (newX[i]-m[i])*betas[i]
        
      
        
    }  ##if(dummies)
    
   else
      
    {
      
 ##??     dummys <-  which(vact_names==row_names[i])
      dummys <-  which(vact==i)
   
      if(isafactor[i]){
        
          
          
          value_one <- which(newX[dummys] ==1)
          
 
          levels <- unlist(xlevels[variable_names[i]])
          
 
          if(length(value_one)==0){
            
            level <- levels[1]
            
            ## clicked on reference category
            addX <- 0
            
          }
          else
          {   
            if(nointercept & firstfactor==i)
            {level <- levels[value_one]}
            else
            {
            level <- levels[1+value_one]}
            addX <-  betas[names(value_one)]           
          }  
         #j <- which(Xuniq[o]==level)
 
          ##write label text in bold (font=2) 
  
          
    plmi <- rep(c(1,-1), times=length(levels))
    ilev  <- which(levels==level)
    plmi <- plmi[ilev]

    #browser()
    ## rather messy trying to overwrite black with red. 
    ## do I need red text?  Maybe fixed position
     posn <- ipos  + 1.5*delta
     ##or suprress altogether??
   ##text(addX,posn,paste(level),cex=0.7,font=2,col=obscol)
   size <- 0.3
     
     ##     text(addX,ipos+0.15,paste(level),cex=0.7,col=obscol,font=2)}
          
          
         #browser()
      }  ##if(isafactor
    else
    {
      
        #!dummies & !isafactor
      addX <-  betas[dummys]*(newX[dummys]-m[dummys])

    }
      
     ##browser()
    
    } ##if(!dummies  
       
## add RED dot point at observation value
 npos <- npanels+1 + make_space
   if(showi | !isinteraction_row[i]){ 
       npanel <- npanel +1
       ipos <- npanel+0.5 +make_space

    cex <- 1
    col <- obscol
    if(isinteraction_row[i]){
      cex <- 1
     ##browser()
      ##pch=4 is a multiplication sign
      points(addX,ipos, col=obscol,cex=cex, pch=4, lwd=2, bg=obscol)
  ## make dot a little smaller if many rows 
        }
   else
     {
    ## if(nrows>10){cex <- 0.8}
    ##browser()
      points(addX,ipos, col=obscol,cex=cex, pch=21, bg=obscol)
      }
   
## also add red point to top of beta X points scale
 #join by vertical line??
    if(droplines){
      ##lighten doesnt seem to work!!
      segments(addX,ipos+0.02,addX,npos+0.38, col="pink") 
    }



}#if(showi | !interaction_row[i]

    if(pointscale){
      
      sumcheck <- sumcheck + 100*( addX-L )/( M-L )
    }
    else
    { sumcheck <- sumcheck +addX}

       

 cex <- 1
 
   
  if(isinteraction_row[i]){
    ##browser()
    ##add a X in red
     points(addX,npos+0.4, col=obscol,cex=cex, pch=4,lwd=2,  bg=obscol)}
 else
 {
     ##browser()
   points(addX,npos+0.4, col=obscol,cex=cex, pch=21, bg=obscol)}
 
     
 }  ## end of loop over row_num
  
       
      
 
 npos <- 1

  if(length(newX) != length(betas)){
## this may be a problem with NAs in beta. 
## Pick out only newX elements referring to beta    
  newX <- newX[is.element(names(newX),names(betas))]
##message("Data inconsistency . Possible beta NAs, attempt to fix")
   }
 
  totS <- sum((newX-m) *betas, na.rm=TRUE) 

  
    
    totSpos <- (M-L)* (totS -Min_tot_score)/Range + L
 ## drawing pointer arrow scale in red
 ## filled diamond pch=23
    sz <- 1.3
    if(npars > 10) {sz <- 0.8}
    
    
  ##position of total score scale
  Total_pos <- yref_pos - 1 
  pos_of_nomo <- yref_pos-1.8
  points(totSpos,Total_pos, col=obscol,cex=sz, pch=23, bg=obscol)
 
  ## draw downward arrow from total scale, arrow cursor 
  ypos <- pos_of_nomo - (ntimes-1)*0.6
  segments(totSpos,Total_pos, totSpos, ypos, col=obscol,cex=1.)
  # w1 <- 0.01*Range
  # w2 <- w1
  ## my own arrowhead function, with intention of varaying width
  ##arrowhead(totSpos,ypos,w1,w2,delta, obscol)
 points(totSpos, ypos + 2*ltick  , col=obscol,cex=sz, pch=25, bg=obscol)
 
   ## total in red obscol
  
  exval <- paste(signif(totS,3),sep="")
  
  if(pointscale) {
    ## need to frig around to get the outputted points value
     exval <- paste(signif(sumcheck,3), sep="")
    }
  ## output total score - exval - in red of red diamond  on the total score axis   
 
 text(totSpos,Total_pos + delta, exval,cex=0.9 ,font=3,col=obscol ) 

 
 ##  add in  red output probabilities/means
 ## loop over ntime - the number of scales. ntimes=1 except
 ## for survival models where different failtimes may have
 ## been requested. 
 for(time in 1:ntimes){ 
  ftime <- tcut[time] 
   
   pos_of_nomo <- yref_pos - 1.8  -(time -1)*0.6 
     ## value position just above the pos-of-nomo-scale
   redscore_ypos <- pos_of_nomo + 0.75*delta
 if(cox){
   ftime <- tcut[time]
  if(odds){
    pred_mean <- (s5[time]^exp(totS))/(1-s5[time]^exp(totS))
    if(fail){pred_mean <- 1/pred_mean}
    }
   else
     {
  ##S(t)^exp(BX)=1-P  
  pred_mean <- s5[time]^exp(totS)
  if(fail) { pred_mean <- 1-pred_mean }
     }
## pos=2  to the left of
  text(totSpos,redscore_ypos, paste(signif(pred_mean,3)),pos=2,cex=0.9,font=3,col=obscol  ) 


  ## must make observation have specified failtime
 
  observation[yvar] <- ftime

    if( interval != "none") pcheck(new_obs,reg,interval, observation,intercept, lm,ols,cox,logistic,poisson,negbin,survreg,pred_mean,
      odds,L,M,Range, Min_tot_score,delta,pos_of_nomo,obscol,dist,ftime,fail,p,s5[time])

  
  
  }  ##if(cox)
##-------------------------------------------------------- 
 
if(survreg){
  if(dist=="loglogistic"){
  pred_mean <- 1/ (1+(exp(-totS-intercept)*ftime)^p) }
 
  if(dist=="lognormal"){
  pred_mean <- 1- pnorm(p*log( exp(-totS-intercept)*ftime) )}
  
  if(dist=="gaussian"){
  pred_mean <- 1- pnorm(p*(-totS-intercept + ftime) )}
  
  if(dist=="weibull"| dist=="exponential"){
  lam <- exp(-totS-intercept)
  pred_mean <- exp(  -(lam*ftime)^p )
  }
  
  if(fail) { pred_mean <- 1-pred_mean }
  if(odds) {pred_mean <- pred_mean/(1-pred_mean)}

  text(totSpos,redscore_ypos, paste(signif(pred_mean,3)),pos=2,cex=0.9,font=3,col=obscol  ) 


 if( interval != "none") pcheck(new_obs,reg,interval, observation,intercept, 
  lm,ols,cox,logistic,poisson,negbin,survreg,pred_mean,odds,
   L,M,Range, Min_tot_score,delta,pos_of_nomo,obscol,
  dist,ftime,fail,p)

  
  }  ##if(survreg)
 
 
   
 if(logistic){
   if(odds){
  pred_mean <- exp(totS+intercept)
    }
   else
   {
  pred_mean <- exp(totS+intercept)/(1 + exp(totS+intercept))
   }
  text(totSpos,redscore_ypos, paste(signif(pred_mean,3)),pos=2,cex=0.9,font=3,col=obscol  ) 
  
     

  if( interval != "none") pcheck(new_obs,reg,interval, observation,intercept, 
  lm,ols,cox,logistic,poisson,negbin,survreg,pred_mean,odds,
   L,M,Range, Min_tot_score,delta,pos_of_nomo,obscol,
  dist,ftime,fail,p)

 
   }  ##if (logistic)



   if(poisson | negbin){
      pred_mean <- exp(totS+intercept)
     text(totSpos,redscore_ypos, paste(signif(pred_mean,3)),pos=2,cex=0.9,font=3,col=obscol  ) 
  
     

if( interval != "none") pcheck(new_obs,reg,interval, observation,intercept, 
  lm,ols,cox,logistic,poisson,negbin,survreg,pred_mean,odds,
   L,M,Range, Min_tot_score,delta,pos_of_nomo,obscol,
  dist,ftime,fail,p)

      }
  
  if(lm | ols){
    pred_mean <- totS+intercept
    text(totSpos,redscore_ypos, paste(signif(pred_mean,4)),pos=2,cex=0.9,font=3,col=obscol  ) 
  if( interval != "none") pcheck(new_obs,reg,interval, observation,intercept, 
  lm,ols,cox,logistic,poisson,negbin,survreg,pred_mean,odds,
   L,M,Range, Min_tot_score,delta,pos_of_nomo,obscol,
  dist,ftime,fail,p)
   
      }

 }  ##for(time in 1:ntimes) 

     
 }  
          
  ##-========================================================
  ## end of  overlay  red (obscol) observation section  
  #==========================================================

##-----------------------------------------------------------------
       

     
   if(FIRSTRUN) {
     if(person){message('Click on graphic expected.  To quit click Esc or press Esc')}
     else
     {message('Click on top menu bar expected.  To quit click Esc or press Esc')
     }
     }
     
  ##=========================================================
  ##function clicked_action  to return action of a graphic mouse click
  ##==========================================================
     
  act <- clicked_action(npanels,row_of_panel,nrows,L,M,S,cplot,dplot,dummies,isadummy,
  isinteraction_row, isafactor,xlevels,row_names,vact,
  vact_names, nointercept,firstfactor,betas,m,person,adddata, nudist,irank,cannotclick,
    make_space,variable_names)

     
 
   esc <- act[[1]]
   cplot <- act[[2]]
   dplot <- act[[3]]
   adddata <- act[[4]]
   nudist <- act[[5]]
   plot_select <- act[[6]]
   Rvar <- act[[7]]
   Rval <- act[[8]]
   new_obs <- FALSE
   if(person & !is.na(Rval) & !is.na(Rvar) & interval!="none"){
     ## update
   if( is.element(Rvar,names(observation)) ){
  ## possibility that if  functions in variable observation cannot be updated   
   ##this possibility not actually possible. Will arise if Rvar is 
  ## a  function 
     
 ## need to  ensure logical kept as logical, i.e. not seeing "TRUE" as char strinbg

     
     if(Rval=="TRUE" | Rval=="FALSE"){
       observation[Rvar] <- Rval=="TRUE"
     } else {
       observation[Rvar] <- Rval
       }
     
     new_obs <- TRUE
      }

   }
##--------------------------------------------
## Escape to close the program. 

if(esc){

  text((L+M)/2,titlepos, title ,cex=1.05 ,font=4) 
  if(points) {  output <- list(points.tables=points_values)  }
  else
  {    output <- paste("points tables not constructed unless points=TRUE ")  }
 return(output)

}  #if(esc)

  
 
  FIRSTRUN <- FALSE
  
}  ##if(ANIMATE)

}  ##WHILE(ANIMATE)
  
}  ##=================end of regplot function==========================

##------------------------------------------------------------------------

subticks <- function(maintickpos, maintickvalue,npos,stick=0.05, func=NA,par, decide=FALSE){
l <- length(maintickvalue)
## add 4 subtick marks between major ones, unless decide = TRUE.
##  Accounts for function if  func is not NA
## stick is short tick length
n4 <-  c(1,2,3,4)
mtv <- maintickvalue
if( !is.na(func) ) {
  if(func == "exp")      {mtv <- log( maintickvalue )}
  if(func == "expit")    {mtv <- log( maintickvalue / (1-maintickvalue) )}
  if(func == "cox")      {mtv <- log( log(maintickvalue) / log(par[1]) )}
  if(func=="loglogistic"){
    mtv <- -log(( ( (1-maintickvalue)/maintickvalue)^(1/par[1]) ) / par[2])}
  
 
 
 if(func=="lognormal"){
 z <- qnorm(1-maintickvalue)
 mtv <- log(par[2]) -z/par[1]
 } 
 if(func=="gaussian"){
   z <- qnorm(1-maintickvalue)
  mtv <-  par[2] -z/par[1]  }
  
 
if(func=="weibull"| func=="exponential"){
 
 lam <- ((-log(maintickvalue))^(1/par[1]))/par[2]
mtv <- -log(lam)
 
 } 
 }

##  do it piecemeal, filling between major  ticks
for(i in 1:(l-1) ) {
  
  if(decide){
    stv <- pretty(maintickvalue[i:(i+1)])
    
    }
  else
  {
    gap <- ( maintickvalue[i+1] - maintickvalue[i] )/5
    stv <- maintickvalue[i] + n4*gap
  }
  

  
  if(!is.na(func)) {
  if(func == "exp")   {stv <- log(stv)}
  if(func == "expit") {
    ## disallow possible 1 or 0 values

    stv <- stv[which(stv>0 & stv <1)]
    stv <- log(stv / (1-stv))}
    
  if(func == "cox")   {
    
    ## disallow possible 1 or 0 values
    stv <- stv[which(stv>0 & stv <1)]
    stv <- log( log(stv) / log(par[1]) )}        
  
  
  
  
   if(func == "loglogistic"){
  stv <- stv[which(stv>0 & stv <1)]
stv <- -log(( ( (1-stv)/stv)^(1/par[1]) ) / par[2])
 } 
 
 
 if(func=="lognormal"){
  stv <- stv[which(stv>0 & stv <1)]
 z <- qnorm(1-stv)
 stv <- log(par[2]) -z/par[1]
 } 
  
 if(func=="gaussian"){
   stv <- stv[which(stv>0 & stv <1)]
   z <- qnorm(1-stv)
   stv <-  par[2] -z/par[1]  
}
  
 
if(func=="weibull"| func=="exponential"){
  stv <- stv[which(stv>0)]
 lam <- ((-log(stv))^(1/par[1]))/par[2]
 stv <- -log(lam)
} 

  }  
  
  
  
    R <- (maintickpos[i+1] - maintickpos[i])/(mtv[i+1]-mtv[i])
    subtickpos <-  R*(stv - mtv[i]) + maintickpos[i] 
  segments( subtickpos,npos,subtickpos,npos-stick)
}
  
} 
##===========================================================
## pretty values of a probability scale
myprettyP <- function(P){ 
  
  xP <- max(P)
  mP <- min(P)
  
  prettyprobs <- pretty(P,n=8)
  ##browser()

  prettyprobs <- prettyprobs[which(prettyprobs>0 & prettyprobs<1)]
  nticks <- length(prettyprobs)
   while(prettyprobs[1] - mP > 0.0007) {
   
   xpand <- pretty( c(mP, prettyprobs[1:2] ) )
   prettyprobs <- c(xpand,prettyprobs[3:nticks]) 
   prettyprobs <- prettyprobs[which(prettyprobs>0 & prettyprobs<1)]
   nticks <- length(prettyprobs)
   }

 
while(xP -prettyprobs[nticks] > 0.0007) {
   xpand <- pretty( c(prettyprobs[(nticks-1):nticks],xP) )
   prettyprobs <- c(prettyprobs[1:(nticks-2)],xpand) 
  prettyprobs <- unique(prettyprobs)
    prettyprobs <- prettyprobs[which(prettyprobs>0 & prettyprobs<1)]
    nticks <- length(prettyprobs)
   
    }

 return(prettyprobs)
}

##=======================================================
## pretty  vsalues of an Odds output scale

myprettyO <- function(O){ 
  
  xP <- max(O)
  mP <- min(O)
  
  prettyodds <- pretty(O,n=8)
  ##browser()

  prettyodds <- prettyodds[which(prettyodds>0 )]
  nticks <- length(prettyodds)
   while(prettyodds[1] - mP > 0.0007) {
   
   xpand <- pretty( c(mP, prettyodds[1:2] ) )
   prettyodds <- c(xpand,prettyodds[3:nticks]) 
   prettyodds <- prettyodds[which(prettyodds>0)]
   nticks <- length(prettyodds)
   }


 return(prettyodds)
}



##================================================================== 
 
 

rfactor <- function(p,fnames,n){
  
  ## function to simulate n obsevations of a categorical 
  ##  variable, by using the multinomial function
  ## 
ncat <- length(p)

x <-  rmultinom(n=n, size=1, prob=p)
## produces matix the wrong way roun. transpose (t)
x <- as.matrix(t(x))

colnames(x) <-  fnames
Fx <- vector(length=n)

for(i in 1:ncat){
  j <- which(x[,i]==1)
  Fx[j] <- colnames(x)[i]
}


Fx <- factor(Fx, levels = fnames)
return(Fx)
}
##===================================================
## function to extract first word [1:1} of a string 
string_fun <- function(x) {
  ul = unlist(strsplit(x, split = "\\s+"))[1:1]
  paste(ul,collapse=" ")
}

##======================================================================

Total_distribution <- function(L,M,pointscale, tot_score,none, 
  dplot,cplot,nrows,ltick, stick, center,delta,dencol,boxcol,spkcol,tickl,aspect,
  yref_pos,showsubticks,allfact,plot_select){
 
# -------------------------------------------------------- 
# function for total score  distribution
# -------------------------------------------------------- 
 ##  add  Total score axis 
  ## font=4 is bold italic =3 is italic
  
  

  if(pointscale){
  text(L,yref_pos,paste("Total-points-to-outcome nomogram:"),font=4)}
  else
  {text(L, yref_pos,paste("Total-score-to-outcome nomogram:"),font=4)}
  
  
## now force the position of the box plot to coincide with L, M at the end points

Max_tot_score <- max(tot_score)
Min_tot_score <- min(tot_score)
Range <- Max_tot_score-Min_tot_score



scale_score <- (M-L)* (tot_score -Min_tot_score)/Range + L

##position of distribution at npos+1=2
## position of Total scale at npos+0.5=1.5
##  make npos the centre line of each  box, violin plot 
npos <- yref_pos - 0.5

# 
#   
#   if(plot_select == 2) {cplot <- "none" 
#                         dplot <- "none"}
#   if(plot_select == 6)  cplot <- "boxplot"
#   if(plot_select == 9)  cplot <- "bean"
#   if(plot_select == 8)  cplot <- "violin"
#   if(plot_select == 7)  cplot <- "cumul"
#   if(plot_select == 3)  cplot <- "density"
#   if(plot_select == 4)  dplot <- "boxes"
#   if(plot_select == 5) {dplot <- "spikes" 
#                         cplot <- "spikes"}          
## use  just click plot selection

if(plot_select !=2){
   
 # 
 ## use boxes of total score if all factors

  if(plot_select==4){
    frq <- as.numeric(table(scale_score))
    
    tot <- sum(frq)
    score_uniq <- as.numeric(names(table(scale_score)))
    luniq <- length(score_uniq)
 
    
    dy <- 0.35*sqrt(frq/tot) 
    #dx <- dy*scale
    dx <- dy * aspect
    xleft <-  score_uniq - dx
    xright <- score_uniq + dx
    ybottom <- c(rep(npos,times=luniq)) - dy
    ytop <- c(rep(npos,times=luniq)) + dy

    # also add a faint axis ??
    ipos <- (ytop+ybottom)/2
       segments( min(xleft),ipos,max(xleft),ipos,col="gray") 

    
    rect(xleft, ybottom, xleft+2*dx, ytop,border="blue",col=darken(boxcol))
    
    
    points(0.5*(xleft+xright),(ybottom+ytop)*0.5,cex=0.3, col="black")
     
    
  } ##if(plot_select==4 )  
  
  
  else 
    
  {  
    #browser()
  score_uniq <- unique(scale_score)
  luniq <- length(score_uniq)
 
  bandwidth <- 0.01*(max(scale_score)-min(scale_score))
    if(cplot=="violin"){
    
    
    box <- vioplot(scale_score, add=TRUE,at=npos,horizontal=TRUE,  wex=0.7, h=bandwidth, 
                   colMed="black", pchMed=21, border="blue",col=darken(dencol))
  }
   if(cplot=="bean") {
     

     box <-  beanplot(scale_score, add=TRUE,log="", at=npos,horizontal=TRUE,wex=0.8, bw=bandwidth, 
                    what = c(0,1,0,1),maxwidth=0.5,ll=0.1,overallline="median",beanlinewd=0, method="overplot",
                     kernel="epanechnikov",  axes=FALSE, border="blue", col=darken(dencol))
             pars = list(boxwex = 0.4, staplewex = 0.5, outwex = 0.5,cex=0.5) }                  
  if(cplot=="boxplot"){

  box <-  boxplot(scale_score, add=TRUE,at=npos,horizontal=TRUE,  border="blue",col=darken(dencol),
                  pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.5,cex=0.5))
      
  }
  
  if(cplot=="density") {
    ## kernel density
    dens <- density(scale_score,kernel="epanechnikov", bw=bandwidth)
    dens$y <- 0.6* dens$y/max(dens$y) +npos-0.5
 ##   lines(dens$x,dens$y,col="blue")
    zero <- c(rep(npos-0.5,times=length(dens$x)) ) 
    
  ## fill in with handy bit of code 
    
    tcol <- darken(dencol)
    polygon(c(dens$x, rev(dens$x)), c(dens$y, zero),
            col = tcol, border = NA)
    lines(dens$x,dens$y,col="blue")
    
   
  }
    
if(cplot=="cumul") {
  ipos <- npos-0.5
  orderS <- sort(scale_score)
  norderS <- length(orderS)
  ecdy <- (1:norderS)/norderS
  yspace <- 0.7
  ecdy <- yspace* ecdy +ipos

   n20 <- floor(norderS*0.2)
  ypos <- orderS[n20]
 lines(orderS,ecdy,type="s", col="blue")

 segments(ypos,ipos,ypos,ipos+yspace)
 text(ypos - tickl,ipos+yspace, paste("1"),cex=0.6,adj=0)
 text(ypos - tickl,ipos+0.5*yspace, paste("0.5"),cex=0.6,adj=0)

 segments(ypos, ipos+yspace, ypos +tickl,ipos+yspace)
  segments(ypos, ipos+0.5*yspace, ypos +tickl,ipos+0.5*yspace)
  segments(ypos, ipos+yspace, orderS[norderS],ipos+yspace,col="gray")
segments(ypos, ipos+yspace*0.5, orderS[norderS],ipos+yspace*0.5,col="gray")

    }

  
  
  if(cplot=="spikes"){
    ## frequency table - will be ordered by  increasing S 
    ##browser()
    frq <- as.numeric(table(scale_score))
    
    score_uniq <- as.numeric(names(table(scale_score)))
   
 ##  before I used this to get values, but found instances where
 ## table() and unique() return different length.
 ## score_uniq <- sort(  unique(scale_score)  )
    ##(total scores)
    maxfreq <- max(frq)
   ## if(maxfreq==1){message("Note: total scores all unique, value frequency=1")}
    
    dy <- 0.7*frq/maxfreq 
    dx <- dy * aspect

    xleft <-  score_uniq - dx
    xright <- score_uniq + dx
    luniq <- length(score_uniq)
    ybottom <- c(rep(npos-0.5,times=luniq)) 
    ytop <- c(rep(npos-0.5,times=luniq)) + dy
    
    ## spikes of continuous  
    for(nspike in 1:luniq){
     
      segments( score_uniq[nspike], ybottom[nspike], score_uniq[nspike], 
        ytop[nspike], lwd=3,col=spkcol)
    }
  
  }  ##if(cplot=="spikes")
  
 
  }  ## else of ##if(dplot=="boxes" )
    
  
  
}

#if pointscale adust scale to output pretty points
tickvalues <- pretty(tot_score,n=6)

if(!pointscale){
tickv <- tickvalues
tickp <- (M-L)* (tickv -Min_tot_score)/Range + L
## fix positions of the scale relative to values of pretty values
## draw line to fit L-M boundary, ticks having  screen coordinates
}
else
{


  tickv <- tickvalues
  
  ## convert to points, but need to account for there being sum of nrows parameters
  tickv <- 100*(tickv -L*nrows)/(M-L)
  tickv <- pretty(tickv,n=8)
  ## convert back into actual total scores
  pts  <- ((M-L)*tickv)/100 +L*nrows
  tickp <- (M-L)* (pts -Min_tot_score)/Range + L 
} 

Total_pos <- yref_pos - 1 
  limits <- which(tickp > L-0.25*(M-L) & tickp < M+0.15*(M-L) )
  tickp_lim <- tickp[limits]
  tickv_lim <- tickv[limits]
segments(min(tickp_lim),Total_pos,max(tickp_lim),Total_pos)
segments(tickp_lim,Total_pos,tickp_lim,Total_pos - ltick)
text(x=tickp_lim,y= Total_pos -delta ,paste(tickv_lim),cex=0.7)
if(showsubticks) subticks(tickp_lim, tickv_lim,Total_pos,stick)

tickvalue <- tickv

if(pointscale) {ex <- "Total points "}
else
{if(center) 
     {ex <- expression(italic(paste("Total ", beta,"(X-m) ")))}

   else
      {ex <- expression(italic(paste("Total ", beta,"X "))) }
  }

text(min(tickp_lim),adj=1,Total_pos, ex,cex=0.8 ,font=3 ) 

}## end of function Total_distribution
##====================================================

clicked_action <- function(npanels,row_of_panel,nrows,L,M,S,cplot,dplot,dummies,isadummy,
  isinteraction_row, isafactor,xlevels,row_names,vact,
  vact_names,nointercept,firstfactor,betas,m,person,adddata, nudist,irank,cannotclick,
  make_space,variable_names){
##--------------------------------------------------------
## function clicked_action
## function to return action required for a mouse click
## if clicked on new  values adddata is returned
##------------------------------------------------------- 
     Rvar <- NA
     Rval <- NA
     
     gap <- (npanels+2 + make_space)/100
     ## position just down from top
  ## create "menu bar" of graphic choices 
     ##browser()
     dialogpos <- npanels+2+ make_space + 1.5*gap
     
     xgap <- (M-L)/10
     rect( L, dialogpos-gap, L+xgap*9.5, dialogpos+gap,border="blue",col="white")
     
     text( L+1*xgap,dialogpos,  paste("Esc"),cex=0.7,col="blue")                          
     text( L+2*xgap,dialogpos,  paste("none"),cex=0.7,col="blue")
     text( L+3*xgap,dialogpos,  paste("density"),cex=0.7,col="blue")
     text( L+4*xgap,dialogpos,  paste("boxes"),cex=0.7,col="blue")
     text( L+5*xgap,dialogpos,  paste("spikes"),cex=0.7,col="blue")
     text( L+6*xgap,dialogpos,  paste("boxplot"),cex=0.7,col="blue")
     text( L+7*xgap,dialogpos,  paste("ecdf"),cex=0.7,col="blue")                          
     text( L+8*xgap,dialogpos,  paste("violin"),cex=0.7,col="blue")
     text( L+9*xgap,dialogpos,  paste("bean"),cex=0.7,col="blue")
     
    ## click on graphic to locate one (n=1) x,y coordinate
    
    XY <- locator(n=1)
    
    panel <- floor(XY$y -  make_space ) 
   
  esc <- FALSE
  plot_select <- 2
   ##browser()
  
  #Press Esc to escape
   if(length(panel)==0){
     esc <- TRUE
                    }
  
  else

  {
   
 ## plot_select <- floor( (dialogpos-XY$y)/gap + 1.5)
    plot_select <- floor((XY$x-L)/xgap +0.5)

  
 ## click on dialog area
    
    currentplot <- c(cplot,dplot)

  if(XY$y > dialogpos - gap & plot_select <=9 & plot_select >= 1 ){ 
    
  if(plot_select == 2) {cplot <- "none" 
                        dplot <- "none"}
  if(plot_select == 6)  cplot <- "boxplot"
  if(plot_select == 9)  cplot <- "bean"
  if(plot_select == 8)  cplot <- "violin"
  if(plot_select == 7)  cplot <- "cumul"
  if(plot_select == 3)  cplot <- "density"
  if(plot_select == 4)  dplot <- "boxes"
  if(plot_select == 5) {dplot <- "spikes" 
                        cplot <- "spikes"}                     
 
  
  # if(plot_select == 3) {dplot <- "boxes"
  #                       if(!distribution ) {cplot <- "boxplot"}
  # }
  if(plot_select == 1)   {  esc <- TRUE }
   
##   if(cplot =="violin"){library(vioplot)} 
##   if(cplot =="bean")  {library(beanplot)} 
  nudist <- !identical(currentplot,c(cplot,dplot) )
  nudist <- TRUE
  
 
  }  ##if(XY$y > dialogpos - gap & plot_select <=9 & plot_select >= 1 
  #XXX
  else
  ## click within rows of plot if person 
    
   {
     
    
   if(!person  | cannotclick  ){
  
  if( panel <=0 ){ message("Cannot click in Total nomogram area") }  
  if(!person){message("observation is NULL. Clicking variable has no effect")}
  if(cannotclick){message("Cannot click for new observation")}
     

   }  
     else
       
   {  
     person <- TRUE 
    if( panel <=0 | panel>npanels){
   message("Cannot click here")   
      }
     
     else
    
     {
       row <- row_of_panel[panel]
       row <- irank[row]
 

       if(row <= nrows ){
    
     #browser()
    if(!isinteraction_row[row]){
  
   if(!dummies){
    
    
    
    if(isafactor[row]) {
     
    Xuniq <- unlist(xlevels[row_names[row]])
   ## browser()
    if(nointercept & firstfactor == row){
      beta_values <- betas[which(vact_names==row_names[row])]
      
###      eps <- abs(XY$x-beta_values )
  ##   eps <- XY$x-beta_values
      eps <- abs(XY$x-beta_values )
   ##browser()
      
    }
    else  ##if(nointercept & firstfactor == row
    {
    
    beta_values <- betas[which(vact_names==row_names[row])]
    ##  Need to order with 0 added
     
    ## find closest point clicked i.e.  snap to point
    eps <- abs(XY$x-c(0,beta_values) )
    
    
    } 
    
    
    oeps <- order(eps)
    
  ## update adddata with clicked value   
    
    adddata[row_names[row] ]  <- Xuniq[oeps[1]]
 ##  output clicked on value    
  message("clicked: ",paste(row_names[row],"=",Xuniq[oeps[1]]))
    ##browser()
 ##return with varaible and value   
    Rvar <- row_names[row]
    Rval <- Xuniq[oeps[1]]
   
    }
    
    else  ##if(isafactor[row])
    
      { nfact <- which(vact==row)
      value <- XY$x/betas[nfact]
   ##update adddata with clicked value   
       var <- row_names[row]
      value <- value +m[nfact]

      ###THIS FAILS if(factr[nfact]) value <- round(value)

      adddata[var] <- value
       message("clicked: ",paste(var,"=",signif(value,3) ))
       
       Rvar <- var
       Rval <- value
       ##browser()

   } 
    
   
}

else
 ##if(!dummies)................... 
{
  ##do for dummies=TRUE
  
    value <- XY$x/betas[row]

    value <- value +m[row]
    


    if(isafactor[vact[row]]){
      
## must be either 0 or 1, whichever is "nearest"
      if(value>0.5){
        value <- 1}
      else {value <-0} 
      
      
      }
      ##var <- colnames(S)[row]
    var <- row_names[row]
  
if(isadummy[row]){
  var_name <- variable_names[vact[row]]
  level <- substring(var,nchar(var_name)+1,nchar(var))
  
  
  ## if clicked on 0, the reference level has been chosen, which is first item of unlisted xlevels
  
 if(value == 0) level <- unlist(xlevels[var_name])[1]
 ## observation[var_name] <- level
 adddata[var_name] <- level
message("clicked: ",paste(var_name,"=",level))

  Rval <- level
  Rvar <- var_name
}
 else
 {
  
   
 ##if(factr[row]) value <- round(value)  
 adddata[var] <- value
  
message("clicked: ",paste(var,"=",signif(value,3)))
 Rvar <- var
 Rval <- value
 
 }

}  ##else of if(!dummies) 

}  ##if(!interaction

else
  
{
  
  message("Cannot click on an interaction") 
  

}

}  ##if(row>=1 & <=npars
else
  
{
  
 ## esc <- TRUE
  ## do nothing  
  
}


} 
       
       
} # else of if(!person)
}  ##else of XY$x < L & plot_select <=5 & plot_select >= 1
}  ##if(length(row)==0 i.e. pressed Esc
  
  ## erase dialog menu bar with a white-over box slightly bigger

  rect( L-0.1*xgap, dialogpos-gap, L+10*xgap, dialogpos+gap,
  border="white",col="white")


 act <- list(esc,cplot,dplot,adddata,nudist,plot_select,Rvar,Rval)
return(act)
   

} ##end of clicked_action function

##======================================================================
Poisson_scale  <- function(tot_score,intercept,L,M,Range,Min_tot_score,delta,
  ltick,stick,ypos,yvar,negbin,showsubticks){
 ## tickvalues are pretty values of total beta score
  ## corresponding   exp(betXa+b0) 
  mean <- exp(tot_score + intercept)
   ## make nice pretty correponding values  of the exp() scale
  prettymean <- pretty(mean,n=20)
  
  ## try filling in lower end of scale
   if(prettymean[1]==0){ 
     prettymean[1] <- signif(min(mean),1)
   }

  prettyscores <- log(prettymean) -intercept
  
    
  ## corresponding plotting positions on L,M
  tickpos <- (M-L)* (prettyscores  - Min_tot_score)/Range + L

  ## also limit by upper value of counts 
   mxcount <- exp(intercept+max(tot_score))
   tickpos <-  tickpos[which(prettymean< mxcount)]  
   prettymean <- prettymean[which(prettymean< mxcount)]

    
 ## tickpos <-  tickpos[limits]  
  ##prettymean <- prettymean[limits]
  
sievePr <- sieve(tickpos,prettymean, 0.05*(M-L))

Xpos <- unlist(sievePr[1])
Pr <-   unlist(sievePr[2])

if(!showsubticks) tickpos <- Xpos

  
  limits <- which(tickpos > L-0.3*(M-L) & tickpos < M+0.1*(M-L) )
  
  segments(min(tickpos),ypos,max(tickpos),ypos)
  segments(tickpos,ypos,tickpos,ypos-ltick)
  if(showsubticks) subticks( tickpos,prettymean,ypos,stick,func="exp")


  
##sievePr <- sieve(tickpos,prettymean, 0.05*(M-L))


text(x=Xpos,y= ypos - delta,paste(signif(Pr,3)),cex=0.7)

 ## text(x=tickpos,y= npos+0.6,paste(signif(prettymean,3)),cex=0.7)  
  
 if(negbin)
 {text(min(tickpos),ypos , adj=1, paste("Negbin ",yvar," "),cex=0.8,font=3  )}
else
{text(min(tickpos),ypos , adj=1, paste("Poisson ",yvar," "),cex=0.8,font=3  )}
  ## text((max(tickpos)+min(tickpos))/2,ypos -2.5*delta , paste("Poisson mean",yvar),cex=0.8,font=3  ) 
      
  
}  ## Poisson_scale

##=====================================================================

Lm_scale <- function(tot_score,intercept,L,M,Range,Min_tot_score,delta,
  ltick,stick,ypos,yvar,showsubticks){
  ##npos <- 0.5
  tickmean <- tot_score+intercept
  prettyscores <- pretty(tickmean,n=8)
  
  
  ## tick positions  correspond to scale of the boxplot and its position
  ## after subtracting out the  intercept
  
  tickpos <- (M-L)* (prettyscores -intercept - Min_tot_score)/Range + L
  
      
  
limits <- which(tickpos > L-0.3*(M-L) & tickpos < M+0.1*(M-L) )
tickpos <-  tickpos[limits]
prettyscores  <-  prettyscores[limits]


  

  
  
  
  ## draw line to fit L-M boundary
   
  segments(min(tickpos),ypos,max(tickpos),ypos)

  segments(tickpos,ypos,tickpos,ypos-ltick)
  
 
  
  text(x=tickpos,y= ypos -delta,paste(signif(prettyscores,3)),cex=0.7)                                     
  text(min(tickpos), adj=1,ypos, paste("mean",yvar," "),cex=0.8 ,font=3 ) 
  
   if(showsubticks) subticks(tickpos, prettyscores,ypos,stick)
  
}  ## end Lm_scale
##=======================================================================
 Logistic_scale <- function(tot_score,intercept,L,M,Range,Min_tot_score,delta,
  ltick,stick,ypos,yvar,odds,showsubticks){

E <- exp(tot_score+intercept)
if(odds){
 prettyprobs <- myprettyO( E ) 
 
 prettyscores <- log(prettyprobs) - intercept
 
 
}
else
{
prettyprobs <- myprettyP( E/(1+E))

prettyscores <- log(prettyprobs/(1-prettyprobs)) -intercept
}

tickpos <- (M-L)* (prettyscores  - Min_tot_score)/Range + L
  
limits <- which(tickpos > L-0.3*(M-L) & tickpos < M+0.1*(M-L) )
tickpos <-  tickpos[limits]
prettyprobs  <-  prettyprobs[limits]


 
sievePr <- sieve(tickpos,prettyprobs, 0.05*(M-L))

Xpos <- unlist(sievePr[1])
Pr <-   unlist(sievePr[2])

if(!showsubticks) tickpos <- Xpos

segments(min(tickpos),ypos,max(tickpos),ypos)
segments(tickpos,ypos,tickpos,ypos-ltick)

## use sieve function to avoid text overlay on axis

##sievePr <- sieve(tickpos,prettyprobs, 0.05*(M-L))
#logistic
text(x=Xpos,y= ypos - delta,paste(signif(Pr,3)),cex=0.7)

if(odds){
text(min(tickpos), adj=1,ypos , paste("Odds(",yvar,") ") ,cex=0.8 ,font=3)  
if(showsubticks) subticks( tickpos,prettyprobs,ypos,stick, func="exp", decide=TRUE)

}
else
{  
text(min(tickpos), adj=1,ypos , paste("Pr(",yvar,") ") ,cex=0.8 ,font=3)
if(showsubticks) subticks( tickpos,prettyprobs,ypos,stick, func="expit", decide=TRUE)

}

}  ## end Logistic_scale 
 ##===============================================================

 Cox_scale <- function(tot_score,intercept,L,M,Range,Min_tot_score,delta,
  ltick,stick,ypos,yvar,odds,baseS,h,s5,betas,center,m,bm,tcut,fail,showsubticks){

  ## possibility that hazard not computable - NAs in betas? 
  if(is.null(baseS ) ){

   sumexp <- exp( sum(betas*m,na.rm=TRUE) )
    hz <- h$hazard*sumexp 
 
## user specified failtime. Select hazard "close" to tcut
    st <- sort(h$time)
    ibelow <- max(which(st<=tcut ))         
   ## h5 <- hz[c(ibelow, ibelow+1)]
    h5 <- hz[ibelow]
    s5 <- exp(-h5)
  
 ##  take default mean-time unit baseline ??
 ##s5 <- mean(exp(-h5))

  }
  
  else
  {## input baseline survival
    ## need to scale if centred graphic
  if(center){
    s5 <- baseS^exp(sum(bm))}
    else
    {s5 <- baseS}
}
 ## get pretty survival scores
 
 scut <- s5 ^ exp(tot_score)
 
 if(odds){
   ##  surviavl odds odds
   sodds <- (scut)/(1-scut)
   prettyX <- myprettyO(sodds)
   
   ##retain "prettyprobs" even though odds
   probs <- prettyX/(1+prettyX) 
   
   exp_score <- log(probs)/log(s5)
  ## corresponding  betaX scores are 
    prettyscore <- log(exp_score)

 }
 else
 {
 
 prettyX <- myprettyP( scut)

  ## use S(t)^exp(BX)=1-P to get scores corresponding to probabilities
 exp_score <- log(prettyX)/log(s5)
  ## corresponding  betaX scores are 
 prettyscore <- log(exp_score)
 }
 
 tickpos <- (M-L)* (prettyscore  - Min_tot_score)/Range + L
 
 
 limits <- which(tickpos > L-0.3*(M-L) & tickpos < M+0.1*(M-L) )
 tickpos <-  tickpos[limits]
 
 
sievePr <- sieve(tickpos,prettyX, 0.05*(M-L))

Xpos <- unlist(sievePr[1])
Pr <-   unlist(sievePr[2])

if(!showsubticks) tickpos <- Xpos
 

 segments(min(tickpos),ypos,max(tickpos),ypos)
 segments(tickpos,ypos,tickpos,ypos-ltick)
 
 
 prettyX  <-  prettyX[limits] 
 
 ## subticks on odds scale not yet implemented. 
 if(!odds){
   if(showsubticks) subticks(tickpos,prettyX,ypos,stick,func="cox",par=s5,decide=TRUE)
 }
 
## use sieve function to avoid text overlay on axis
# 
# sievePr <- sieve(tickpos,prettyX, 0.05*(M-L))
# 
# Xpos <- unlist(sievePr[1])
# Pr <-   unlist(sievePr[2])

 ## put values on the axis  Cox 
 if(fail){ chr <- "<" } else { chr <- ">" }
  if(!odds){
     if(fail) {
     text(x=Xpos,y= ypos - delta,paste(signif(1-Pr,3)),cex=0.7)
     }
     else
     {
      text(x=Xpos,y= ypos -  delta,paste(signif(Pr,3)),cex=0.7)
     }
        text(x=min(tickpos), ypos , adj=1, paste("Pr(",yvar,chr, signif(tcut,4),") " ), cex=0.8,font=3) 

  }  
    else ##if(!odds
      
    {
      if(fail){
      text(x=Xpos,y= ypos -  delta,paste(signif(1/Pr,3)),cex=0.7)}
      else
      {
      text(x=Xpos,y= ypos -  delta,paste(signif(Pr,3)),cex=0.7)}
  text(x=min(tickpos), ypos , adj=1, paste("odds(",yvar,chr, signif(tcut,4),") " ), cex=0.8,font=3) 
#    
    }

return(list(s5=s5))
}  ##end Cox_scale
##=====================================================================

Survreg_scale <- function(tot_score,intercept,L,M,Range,Min_tot_score,delta,
  ltick,stick,ypos,yvar,betas,m,tcut,fail,p,dist,odds,showsubticks) 
  
  {
  
 sumexp <- exp( sum(betas*m,na.rm=TRUE) )
   
## lp_lnorm <- predict(reg, type="linear") 
## grabbed formula  Klabfleisch & Prentice p25 (assuming p and lambda as for Weibull)

lambda <- exp(-tot_score-intercept)

if(dist=="loglogistic"){
scut <- 1/(1+(lambda*tcut)^p)}


if(dist=="lognormal"){
##scut <- 1- pnorm(p*log(lambda*tcut))
 scut <- 1- pnorm(p*(-tot_score - intercept +log(tcut)) )}

if(dist=="gaussian"){
scut <- 1- pnorm(p*(-tot_score-intercept+tcut))}


if(dist=="weibull" | dist=="exponential"){
 scut <-  exp(-(tcut*lambda)^p)
}

 ## get pretty survival scores and surviaval odds
 if(odds){
   if(fail) {
   prettyX <- myprettyO((1-scut)/scut )
   probs <- 1/(1+prettyX)
   }
   else
   {
   prettyX <- myprettyO(scut/(1-scut) )
   probs <- prettyX/(1+prettyX)
   }
}
else
{prettyX <- myprettyP( scut)
  probs <- prettyX}
 

 ## fix position of lower nomogram scale 
 ##npos <- 0.5
 
  if(dist=="loglogistic"){
 
 prettyscores <- -log(( ( (1-probs)/probs)^(1/p) ) / tcut)

 } 
 
 
 if(dist=="lognormal"){
 z <- qnorm(1-probs)
 prettyscores <- log(tcut) -z/p
 } 
 if(dist=="gaussian"){
   z <- qnorm(1-probs)
   prettyscores <-  tcut -z/p  }
  
 
if(dist=="weibull"| dist=="exponential"){
 
 lam <- ((-log(probs))^(1/p))/tcut
 prettyscores <- -log(lam)
 
 
 } 
 
 

 
 
 tickpos <- (M-L)* (prettyscores -intercept  - Min_tot_score)/Range + L
 
 
  limits <- which(tickpos > L-0.3*(M-L) & tickpos < M+0.1*(M-L) )
 tickpos <-  tickpos[limits]
 
 
 prettyX  <-  prettyX[limits]
 
  
sievePr <- sieve(tickpos,prettyX, 0.05*(M-L))

Xpos <- unlist(sievePr[1])
Pr <-   unlist(sievePr[2])

if(!showsubticks) tickpos <- Xpos



 ## main tick scale
 segments(min(tickpos),ypos,max(tickpos),ypos)
 segments(tickpos,ypos,tickpos,ypos-ltick)
 
 ##
 ## subticks on odds scale not yet implemented. 

 if(!odds) {
   if(showsubticks) subticks(tickpos,prettyX,ypos,stick,func=dist,par=c(p,tcut),decide=TRUE)
   }


## use sieve function to avoid text overlay on axis
# 
# sievePr <- sieve(tickpos,prettyX, 0.05*(M-L))
# 
# Xpos <- unlist(sievePr[1])
# Pr <-   unlist(sievePr[2])
# 
 
 # 
 if(fail){ chr <- "<" } else { chr <- ">" }


  if(!odds){
     if(fail) {
     text(x=Xpos,y= ypos - delta,paste(signif(1-Pr,3)),cex=0.7)
     }
     else
     {
      text(x=Xpos,y= ypos -  delta,paste(signif(Pr,3)),cex=0.7)
     }
      text(x=min(tickpos), ypos , adj=1, paste("Pr(",yvar,chr, signif(tcut,4),") " ), cex=0.8,font=3) 

  }  
    else ##if(!odds
      
    {
      
      text(x=Xpos,y= ypos -  delta,paste(signif(Pr,3)),cex=0.7)
      
      
  text(x=min(tickpos), ypos , adj=1, paste("odds(",yvar,chr, signif(tcut,4),") " ), cex=0.8,font=3) 
#    
    }

}  ## end Survreg_scale 
##=========================================================

 points_tables <- function(survreg,dist,p,cox,s5,fail,tcut,poisson,
   negbin,lm,ols,logistic,pts,intercept,
   npars,yvar,tickv) 

{
 ## point_table to create last total points to outcome scoring table.   
   if(logistic){
        lprob <- exp(pts+intercept)/(1 + exp(pts+intercept))
        pointsframe <-  as.data.frame( cbind(tickv,signif(lprob,4) ) )
        colnames(pointsframe) <- c("Total Points",paste("Pr(",yvar,")")) 
   }
   
   if(lm|ols){
     
        pointsframe <-  as.data.frame( cbind(tickv,signif(pts+intercept,4) ) )
        colnames(pointsframe) <- c("Total Points", paste("mean",yvar)) 
   }
   
   if(poisson | negbin){
     
        pred_mean <- exp(pts+intercept)
        pointsframe <-  as.data.frame( cbind(tickv,signif(pred_mean,4)) )
        if(poisson){
        colnames(pointsframe) <- c("Total Points",paste("Poisson mean",yvar)) 
        }
        else
        {colnames(pointsframe) <- c("Total Points",paste("Negbin mean",yvar))
        }
   }
   
   if(cox){
        if(fail){chr <- "<"} else {chr <- ">"}
        pred_mean <- signif(s5^exp(pts),4)
        if(fail) { pred_mean <- 1-pred_mean }
        pointsframe <-  as.data.frame( cbind(tickv,signif(pred_mean,4) ) )
        colnames(pointsframe) <- c("Total Points",paste("Pr(",yvar,chr, signif(tcut,4),")" )) 
   }
   
   if(survreg){
  if(fail){ chr <- "<" } else { chr <- ">" }      
  if(dist=="loglogistic"){
  pred_mean <- 1/ (1+(exp(-pts-intercept)*tcut)^p) }
   if(dist=="lognormal"){
  pred_mean <- 1- pnorm(p*(-pts - intercept +log(tcut)) )} 
  if(dist=="gaussian"){
  pred_mean <- 1- pnorm(p*(-pts-intercept + tcut) )}
  if(dist=="weibull"| dist=="exponential"){
  lam <- exp(-pts-intercept)
  pred_mean <- exp(  -(lam*tcut)^p )}

        if(fail) { pred_mean <- 1-pred_mean }
        pointsframe <-  as.data.frame( cbind(tickv,signif(pred_mean,4)) )
        colnames(pointsframe) <- c("Total Points",paste("Pr(",yvar,chr, signif(tcut,4),")" )) 
        
   }


   
   return( pointsframe)
 }  #en points_table
##=========================================================
## augment variable names  for  function regression poly() and bs() 
## Defunct
 augv <- function(func, vnames, cnames){
  
##NB:  grep doesnt seem to like "poly(" with (  because not regular expression 
##ANS: need \\C trick. 
  num <- length(vnames)
cgrep <- grep(pattern=func,cnames)
vgrep <- grep(pattern=func,vnames)
## treat all polynomial terms as a "variable". Current variable_names
## has  just "poly(x,n)". Need to slot in  terms for each polynomial
## this will probably fail with more than one poly() exression in the model!!
if(vgrep[1] >1){
vnames <- c(vnames[1:(vgrep[1]-1)],
                    cnames[cgrep],
                     vnames[(vgrep[1]+1):num])
}
else
{
  ##nb:   R works out x[1:0] to be x[1] so need separate expression
  vnames <- c(cnames[cgrep],
                   vnames[(vgrep[1]+1):num])
}

return(vnames)
}
##===================================================
##  strip away package name if attached to variables
## defunct function
strippkg <- function(v){
  
   for(i in 1:length(v)){
     if(length(grep("::",v[i],fixed=TRUE))>0){
      p <- instr( v[i], "::" )
      
      v[i] <- substr(v[i],p+2,nchar(v[i]) )
     }
   }
  return(v)
}
  
##-==================================================================
## Utility functions to  extract variable from  function e.g  "age" from "log(age)"
##  instr grabbed from a web page and adapted  find positions of "(" and ")"
instr <- function(str1, str2, startpos=1, n=1){
  
  aa <- unlist(strsplit(substring(str1,startpos),str2))
 
  if(length(aa) < n ) return(0);
  return(sum(nchar(aa[1:n])) + startpos+(n-1)*nchar(str2) )
}
##============================================================
getX <- function(in_str){
  
  ## extract  variable name from function e.eg "age" from "log(age)"
  ## need double backslash for special "(" and ")" characters
  ## str is a (possibly vector) of  characters
  ## special case "(xxx)" excluded (with left >1 condition)
  str <- in_str
  
   nelem <- length(str)
   func <- rep("",times=nelem)
  
  for(i in 1:nelem){
  
  left  <- instr( str[i], "\\(" )
  right <- instr( str[i], "\\)" )
  

if( right > left ) {str[i] <- substr(str[i], left+1, right-1)

## eg pick ip  log( in log(X)
## if  logical  (X>0) this is ( 
## paste in righthand )
  func[i] <-  paste0(substring(in_str[i],1,left),")") }
  
  }
  
  outgetX <- c(str,func)
  
  return(outgetX)
}
##======================================================================
## function to  delete (sieve out) too many tick axis marks 
## usage sieve(tickpos,prettymean, 0.05*(M-L))
## separation must be > c
sieve <- function(x,v,c){
  n <- length(x)
  i <- 2
  xp <- x[1]
  keep <- logical(length=n)
  
  while (i <=n) {
    
    d <- abs(x[i] -xp)
   
    if(d >= c) {  xp <- x[i]
                 keep[i] <- TRUE
                 }
    
    else
      {keep[i] <- FALSE}
    
    i <- i+1
    
    }
  keep[1] <- TRUE
  kept <- which(keep==TRUE)
  
  sievedv <- v[kept]
  sievedx <- x[kept]  
  
  return (list(sievedx,sievedv))
}
##============================================================
##grabbed function to darken/lighten  colors of "Total" distribution
darken <- function(color, factor=1.2){
    col <- col2rgb(color)
    col <- col/factor
    col <- rgb(t(col), maxColorValue=255)
    col
}

lighten <- function(color, factor=1.4){
    col <- col2rgb(color)
    col <- col*factor
    col <- rgb(t(as.matrix(apply(col, 1, function(x) if (x > 255) 255 else x))), maxColorValue=255)
    col
}

##======================================================
## function to slot in parameter names as variable names
## to create a new variable_names.
expand_vars <- function(X, variable_names, parameter_names,DF,isafactor,nvars){
     new_variable_names <- NULL
     new_DF <- NULL
     new_isafactor <- NULL

    i <- 0
    for(k in 1:nvars){
    
      if(X[k]){
## note: fixed parameter necessary 
      DF[k]  <- length(grep(variable_names[k],parameter_names,fixed=TRUE) )
      new_variable_names <- 
      c(new_variable_names, parameter_names[(i+1):(i+DF[k])]) 

      new_DF <- c(new_DF,rep(1,times=DF[k]) )
      new_isafactor <- c(new_isafactor,rep(FALSE,times=DF[k]))
      }
      else{
        new_variable_names <- 
         c(new_variable_names, variable_names[k])
        new_DF <- c(new_DF,DF[k])
        new_isafactor <- c(new_isafactor,isafactor[k])
      }
      
      i <- i +DF[k]
    }
    variable_names <- new_variable_names
    nvars <- length(variable_names)
    isafactor <- new_isafactor
    DF <- new_DF

    return(expanded=list(nvars=nvars,isafactor=isafactor,
      DF=DF,variable_names=variable_names))
} 
##=================================================
## extract DF, the number of parameters per variable
##  of an interaction variable
##  v is variable name eg  "sex:age:eth"
## p is list of parameter names which may include 
##  for example  "sexF:age60:ethM" , "sexM:age40:ethM"

pars_per_var <- function(v,p){
  ivars <- unlist(strsplit(v,":"))
  ##order of variable interaction
  ov <- length(ivars)
##establish order o of each parameter interaction term
  lp <- length(p)
  o <- numeric(length=lp)
  for(i in 1:lp){
  o[i] <- length( unlist( strsplit(p[i],":"))) 
  }
## 
   pars <- grep(ivars[1],p[which(o==ov)],fixed=TRUE)
     
  for(iv in 2:ov){
    parsi <- grep(ivars[iv],p[which(o==ov)],fixed=TRUE)
    pars <- intersect(pars,parsi)
  }
  DF <- length(pars)
 
  return(DF)
}
##=======================================================
## use predict() as a check on  regplot calculatred outcomes
## also to add in an interval estimate on the  outcome scale as a 
## red (obscol) overline.
pcheck <- function(new_obs,reg,interval, observation,intercept, 
  lm,ols,cox,logistic,poisson,negbin,survreg,pred_mean,odds,
   L,M,Range, Min_tot_score,delta,pos_of_nomo,obscol,
  dist,ftime,fail,p,s5){
  ##------------------------------------------------------
## for lm() regression,  interval is a predict() option 
    tick <- 0.7*delta
  ##NB: there is code repetition here, but done in interest of clarity. 

    if(lm){

    pr <- predict(reg,newdata=observation,interval=interval)
    if(abs((pred_mean-pr[1])/pred_mean) > 0.001){
      message("warning: regplot outcome differs from predict(). Contact developer")
      message(paste(signif(pred_mean,4))," v ",paste(signif(pred,4)))  }
    
    if(interval != "none"){
       Int <- "CI: "
      if(interval=="prediction") {Int <- "PI: "}
    
         ##  write out CI ? Used for checking
if(new_obs) message( paste(Int),signif(pred_mean,4),"(", paste0(signif(pr[2],4),"-",signif(pr[3],4)),")" )

    ##  get confidence or prediction interval and  draw on the nomogram axis
    lwr <-  (M-L)* (pr[2]-intercept -Min_tot_score)/Range + L
    upr <-  (M-L)* (pr[3]-intercept -Min_tot_score)/Range + L

    segments(lwr,pos_of_nomo, upr, pos_of_nomo, col=obscol,cex=2.,lwd=2)
    segments(lwr,pos_of_nomo+tick,lwr,pos_of_nomo-tick,col=obscol,lwd=2)
    segments(upr,pos_of_nomo+tick,upr,pos_of_nomo-tick,col=obscol,lwd=2)
    }
    }  #if(lm)
  ##--------------------------------------------------------
  ## for ols (glm object)  predict returns a list.
  if(ols){
    
    pr <- unlist(predict(reg,newdata=observation,se.fit=TRUE))
## 
    if(abs((pred_mean-pr[1])/pred_mean) > 0.001){
      message("warning: regplot outcome differs from predict(). Contact developer")
      message(paste(signif(pred_mean,4))," v ",paste(signif(pred,4)))  }
    
    ##  get confidence interval and  draw on the nomogram axis
    if(interval=="confidence" | interval=="prediction"){
      if(interval=="confidence"){
        Int <- "CI: "
      SE <- pr[2]}
      else
      {Int<- "PI: "
      SE <- pr[3]}

    lower <- pr[1] -1.96*SE
    upper <- pr[1] +1.96*SE
    
    
    ##  write out CI ? Used for checking
if(new_obs) message(paste(Int),signif(pred_mean,4),"(", paste0(signif(lower,4),"-",signif(upper,4)),")" )

    lwr <-  (M-L)* (lower-intercept -Min_tot_score)/Range + L
    upr <-  (M-L)* (upper-intercept -Min_tot_score)/Range + L

    segments(lwr,pos_of_nomo, upr, pos_of_nomo, col=obscol,cex=2.,lwd=2)
    segments(lwr,pos_of_nomo+tick,lwr,pos_of_nomo-tick,col=obscol,lwd=2)
    segments(upr,pos_of_nomo+tick,upr,pos_of_nomo-tick,col=obscol,lwd=2)
    }
}  ##if(ols)
##----------------------------------------------------------------
  
    ## 
  if(poisson | negbin ){
    
    pr <- unlist(predict(reg,newdata=observation,se.fit=TRUE))
## 

    if(abs((pred_mean-exp(pr[1]))/pred_mean) > 0.001){
      message("warning: regplot outcome differs from predict(). Contact developer")
      message(paste(signif(pred_mean,4))," v ",paste(signif(pred,4)))  }
    
    ##  get confidence interval and  draw on the nomogram axis
    
   if(interval=="confidence"){
      SE <- pr[2]
 
    low <- pr[1] -1.96*SE 
    upp <- pr[1] +1.96*SE 
    Su <- exp(upp)
    Sl <- exp(low)
    ##  write out CI ? Used for checking
if(new_obs) message("CI: ",signif(pred_mean,3),"(", paste0(signif(Sl,3),"-",signif(Su,3)),")" )

    lower <- pr[1] -1.96*SE -intercept
    upper <- pr[1] +1.96*SE -intercept
    
    
    lwr <-  (M-L)* (lower  -Min_tot_score)/Range + L
    upr <-  (M-L)* (upper  -Min_tot_score)/Range + L

    segments(lwr,pos_of_nomo, upr, pos_of_nomo, col=obscol,cex=2.,lwd=2)
    segments(lwr,pos_of_nomo+tick,lwr,pos_of_nomo-tick,col=obscol,lwd=2)
    segments(upr,pos_of_nomo+tick,upr,pos_of_nomo-tick,col=obscol,lwd=2)
    
  

}
}  ##if(poisson | negbin ){
##---------------------------------------------------------
## for ols (glm object)  predict returns a list.
  if(logistic ){
    
    pr <- unlist(predict(reg,newdata=observation,se.fit=TRUE))
## 
    Q <- pr[1]
    pred <- exp(Q)/(1+exp(Q))
    if(odds){pred <- exp(Q)}
    if(abs((pred_mean-pred)/pred_mean) > 0.001){
      message("warning: regplot outcome differs from predict(). Contact developer")
     message(paste(signif(pred_mean,4))," v ",paste(signif(pred,4)))   }
    
    ##  get confidence interval and  draw on the nomogram axis
   if(interval=="confidence"){
      SE <- pr[2]
## nwed - intercept for graphic position
    lower <- pr[1] -1.96*SE -intercept
    upper <- pr[1] +1.96*SE -intercept
    
    low <- pr[1] -1.96*SE 
    upp <- pr[1] +1.96*SE 
   
    if(odds){
    Su <- exp(upp)
    Sl <- exp(low)
    }
    else
    {
    Su <- exp(upp)/(1+exp(upp))
    Sl <- exp(low)/(1+exp(low))
    }
    ##  write out CI ? Used for checking
if(new_obs) message("CI: ",signif(pred,3),"(", paste0(signif(Sl,3),"-",signif(Su,3)),")" )

    
    lwr <-  (M-L)* (lower  -Min_tot_score)/Range + L
    upr <-  (M-L)* (upper  -Min_tot_score)/Range + L

    segments(lwr,pos_of_nomo, upr, pos_of_nomo, col=obscol,cex=2.,lwd=2)
    segments(lwr,pos_of_nomo+tick,lwr,pos_of_nomo-tick,col=obscol,lwd=2)
    segments(upr,pos_of_nomo+tick,upr,pos_of_nomo-tick,col=obscol,lwd=2)


}
}  ##if(logistic){
##------------------------------------------------------    
 if(survreg){

    
     
    pr <- unlist(predict(reg,newdata=observation,se.fit=TRUE,type="linear"))
      SE <- pr[2]   
    lower <- pr[1] -1.96*SE 
    upper <- pr[1] +1.96*SE 
  
    if(dist=="loglogistic"){
  pred <- 1/ (1+(exp(-pr[1])*ftime)^p) 
  Slower <- 1/ (1+(exp(-lower)*ftime)^p) 
  Supper <- 1/ (1+(exp(-upper)*ftime)^p) 
  
  }
 
  if(dist=="lognormal"){
  pred <-   1- pnorm(p*log( exp(-pr[1])*ftime) )
  Slower <- 1- pnorm(p*log( exp(-lower)*ftime) )
  Supper <- 1- pnorm(p*log( exp(-upper)*ftime) )
  }
  
  if(dist=="gaussian"){
  pred <- 1- pnorm(p*(-pr[1] + ftime) )
Slower <- 1- pnorm(p*(-lower + ftime) )
Supper <- 1- pnorm(p*(-upper + ftime) )
  }
  
  if(dist=="weibull"| dist=="exponential"){
  lam <- exp(-pr[1])
  lamL <- exp(-lower)
  lamU <- exp( -upper)
  pred <- exp(  -(lam*ftime)^p )
 Slower <- exp(  -(lamL*ftime)^p ) 
 Supper <- exp(  -(lamU*ftime)^p )
 }
  
  if(fail) { pred <- 1-pred 
  Slower <- 1- Slower
  Supper <- 1- Supper}

      if(odds){pred <- pred/(1-pred)
      Slower <- Slower/(1-Slower)
      Supper <- Supper/(1-Supper)
      }
    if(abs((pred_mean-pred)/pred_mean) > 0.001){
      message("warning: regplot outcome differs from predict(). Contact developer")
      message(paste(signif(pred_mean,4))," v ",paste(signif(pred,4)))  }
    
    ##  get confidence interval and  draw on the nomogram axis
   if(interval=="confidence"){
     

    lower <- pr[1] -1.96*SE -intercept
    upper <- pr[1] +1.96*SE -intercept
    lwr <-  (M-L)* (lower  -Min_tot_score)/Range + L
    upr <-  (M-L)* (upper  -Min_tot_score)/Range + L

    segments(lwr,pos_of_nomo, upr, pos_of_nomo, col=obscol,cex=2.,lwd=2)
    segments(lwr,pos_of_nomo+tick,lwr,pos_of_nomo-tick,col=obscol,lwd=2)
    segments(upr,pos_of_nomo+tick,upr,pos_of_nomo-tick,col=obscol,lwd=2)
    
    
    if(new_obs) message("CI: ",signif(pred,3),"(", paste0(signif(Slower,3),"-",signif(Supper,3)),")" )
 
 }
    
  
}  ## if(survreg
##---------------------------------------------------
    ## cox model
    
    if(cox ){
  
    
     
    pr <- unlist(predict(reg,newdata=observation,type="expected",se.fit=TRUE))
   ##survival
    pred <- exp(-pr[1])
  
  if(fail) { pred <- 1-pred }
   

      if(odds){pred <- pred/(1-pred)}
    if(abs((pred_mean-pred)/pred_mean) > 0.001){
      message("warning: regplot outcome differs from predict(). Contact developer")
    message(paste(signif(pred_mean,4))," v ",paste(signif(pred,4)))
      
      }
    
    ##  get confidence interval and  draw on the nomogram axis
   if(interval=="confidence"){
      SE <- pr[2]
      
    lower <- pr[1] -1.96*SE 
    upper <- pr[1] +1.96*SE 
    Slower <- exp(-upper)
    Supper <- exp(-lower)
##write out checking    
    if(fail){
       Sl <- 1- Supper
       Su<- 1- Slower}
    else
    {Sl <- Slower
    Su <- Supper}
    if(odds){
      Sl <- Sl/(1-Sl)
      Su <- Su /(1-Su)
    }
    
 ##  write out CI ? Used for checking
if(new_obs) message("CI: ",signif(pred,3),"(", paste0(signif(Sl,3),"-",signif(Su,3)),")" )
    
    if( log(Slower)/log(s5) >0 & log(Supper)/log(s5)>0 ){
    PL <- log(log(Slower)/log(s5))
    PS <- log(log(Supper)/log(s5)) 
   ##   prettyX <- myprettyP( scut)

  ## use S(t)^exp(BX)=1-P to get scores corresponding to probabilities
 ##exp_score <- log(prettyX)/log(s5)
  ## corresponding  betaX scores are 
 ##prettyscore <- log(exp_score)

      
      

    lwr <-  (M-L)* (PL  -Min_tot_score)/Range + L
    upr <-  (M-L)* (PS  -Min_tot_score)/Range + L

    segments(lwr,pos_of_nomo, upr, pos_of_nomo, col=obscol,cex=2.,lwd=2)
    segments(lwr,pos_of_nomo+tick,lwr,pos_of_nomo-tick,col=obscol,lwd=2)
    segments(upr,pos_of_nomo+tick,upr,pos_of_nomo-tick,col=obscol,lwd=2)
    }
    else
    {message("Warning: NaNs in interval computation")}
 }
    
  
} 
    
    
}  ##end of pcheck
