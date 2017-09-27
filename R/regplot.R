#' Plots a regression nomogram showing covariate distribution.
#' @description 
#' \code{regplot} plots a regression nomogram of \code{coxph}, \code{lm} and \code{glm}
#'  regressions. Covariate distributions are superimposed on nomogram scales and the plot
#'	is animated to allow on the fly changes to distribution representation and to 
#'	enable outcome calculation. 
#' @param reg A regression object of either class \code{glm}, \code{lm}, or \code{coxph}  from a fitted model
#'  
#' @param dummies \code{TRUE} to treat dummy indicators of factor variables as distinct binary variables, with their own nomogram row. 
#' Otherwise different categories of the factor are represented on the same row. 
#' @param points If \code{FALSE} the regression scores \eqn{\beta}\eqn{x}  are shown. 
#' Otherwise the scale is represented by "points".  
#' @param center  Produces plot in which regression score contributions of continuous data are  plotted with respect 
#' to mean values.
#' @param observation  An observation whose values
#'  are superimposed on the plot. A data frame whose column names must include variables used
#'   in the regression formula. 
#' @param title A heading title written to the plot
#' @param other A \code{list} of other optional parameters for features of the plot. \code{bvcol} is color fill of box and violin plots, 
#' \code{sqcol} is color fill of frequency boxes, \code{obscol} is color of observation dots,
#' \code{failtime}   specifies the cut-off time for
#' plotting the risk  nomogram of a \code{coxph} regression (if \code{failtime=NA}, cut-off is the median of time variable,
#'  with no account of censoring), \code{prfail=TRUE} if 
#' probability of failure before \code{failtime} is summarised, otherwise after \code{failtime}. 
#' If the number of unique values of covariate is less than
#' \code{ndisc},  discrete frequencies are plotted by scaled boxes. 
#'  \code{droplines} draws faint vertical lines showing score contributions to an observation.
#' @author Roger Marshall <rj.marshall@@auckland.ac.nz> The University of Auckland, New Zealand
#' 
#' @details Creates a nomogram  representation of a fitted regression. 
#' The distribution of covariates in the model, and of the total regression score, are superimposed on the 
#' nomogram scales. Also the values of a particular observation can be superimposed. The plot is active for mouse input 
#' allowing on-the-fly changes to distribution type (boxes, spikes, box plot, density, violin and bean plots)
#' and also to observation values. 
#' Assumes that a regression has been fitted by either \code{glm},  \code{lm} or \code{coxph}. For \code{glm},
#' the supported family/link pairings are: gaussian/identity, binomial/logit, and  poisson/log. 
#' 
#' @examples
#' ## Analysis of pbc data
#' library(survival)
#' data(pbc) 
#' pbc$catbili <- cut(pbc$bili,breaks=c(-Inf, 2, 4, Inf),
#'                labels=c("low","medium","high"))
#' pbc$died <- pbc$status==2
#' ## Fit a Cox survival model
#' pbccox <-  coxph(formula = Surv(time,died) ~  age + catbili + sex + 
#'                  log(copper+1) +stage + trt,data=pbc)
#' obs <- pbc[1,]
#' regplot(pbccox,observation=obs, other=list( failtime = 1825, prfail = TRUE) ) 

#' ## Fit a glm  logistic regression
#' pbcglm <- glm(formula = died ~  age + catbili + sex + copper +stage ,family = "binomial", data=pbc )
#' regplot(pbcglm)

regplot <- function(reg,  dummies=FALSE,
                    center=FALSE,observation=NULL,title=NULL,points=TRUE,
            other=list(bvcol="#DDEEDD", sqcol="#EEEEEE",obscol="red",failtime=NA,
            prfail=TRUE,ndisc=5,droplines=FALSE) )
  {
  
  
  options(locatorBell = TRUE)
  FIRSTRUN <- TRUE
  ANIMATE <- TRUE
  
  cplot="density"
  dplot="boxes"
  
  person <- !is.null(observation)
 
  weighted <- FALSE 
  
  while(ANIMATE){
  
  distribution <- TRUE 
  
  
  if(cplot=="none"){distribution <- FALSE}
   
   
  rtype <- class(reg)[1]
  cox   <- (rtype=="coxph") 
  lm    <- (rtype=="lm" )
  glm   <- (rtype=="glm")
  
  
  
 if(FIRSTRUN){
 
  if(glm | lm ){
    
  ## check on weights model,
    
 
  lw <- ncol(reg$model)
  if(colnames(reg$model)[lw]=="(weights)"){
    message("Replicate weights assumed")
    weighted <- TRUE 
    W <- floor(reg$model[,lw])
   }}
  
   if( cox ){
    
  ## check on weights model,
      
  if(!is.null(reg$weights)){
  message("Replicate weights assumed")
  weighted <- TRUE
  W <- floor(reg$weights)
  
    }}
  
 }  ##if(FIRSTRUN)
  pointscale <- points
  
  if(!( glm | lm | cox )){
    return(message(paste("Cannot do  regplot for class",
    paste("\"",rtype,"\"",sep=""), "Only \"lm\", \"glm\", \"coxph\" regression supported")) )
  }
  
  
  NAcoef <- which(is.na(reg$coefficients))
  
  if(length(NAcoef)>0){
    message("model has NA beta coefficients. No plot done, model assumed questionable")
    return()
  }
  
  
  if(is.null(title)) {title <- paste(rtype,"regression")}

  
  if( is.null(other$bvcol) ) {bvcol <- "#DDEEDD"} else {bvcol<- other$bvcol}
  if( is.null(other$sqcol) ) {sqcol <- "#DDDDEE"} else {sqcol<- other$sqcol}
  if( is.null(other$obscol) ) {obscol <- "red"} else {obscol<- other$obscol}
  if( is.null(other$prfail) ) {fail <- FALSE} else {fail <- other$prfail}
  if( is.null(other$failtime) ) {tcut <- NA} else {tcut <- other$failtime}
  if( is.null(other$ndisc)  ) {ndisc <- 5} else {ndisc <- other$ndisc}
  if( is.null(other$droplines)  ) {droplines <- FALSE} else {droplines <- other$droplines}
  if(!is.null(other)){
   allowed <- c("prfail","failtime","bvcol","sqcol","obscol","ndisc","droplines")

    test <- is.element(names(other),allowed)
    if(length( names(other)[which(!test)])>0 ){
      message(paste("These \"other\" parameters are inadmissible and ignored"))
      message(paste(names(other)[which(!test)],collapse="  "))
    } 
  }
    nomo <- TRUE

  offset <- FALSE
  
  if( !is.null(reg$offset) ){
    
    offset <- TRUE
  ## only way to extract offset var seems to be via reg$call
  ##  item 4 for glm and 3 for lm.  This needs check for always being correct!
    
    if(glm) {offset_varname <- as.character(reg$call)[4]
      }
    else
         {if(lm) {offset_varname <- as.character(reg$call)[3]
          }
          else
            { if(cox){
              message("no offset for cox")
            }
}
         }      

    if(FIRSTRUN) {message("There is an offset in model: " , paste(offset_varname) )
                  
    return(message("Unable to make a plot with an offset in regplot0.1"))
    
                  
    
    if(is.na(offset_varname) ){
    return(message('... but value is NA. Please use "offset=" argument in regression build') )
    }
    
    
    GX <- getX(offset_varname)
    
    ##getX return two halves: first half are variable, second is function
    ## get actaul_vars stripped of possible function (e.g log())
    ##Stripped_variable_names <-  GX[1:(length(GX)/2)]
    
    func_offset_varname  <- GX[(length(GX)/2 +1):length(GX)]
    ##offset_varname  <- GX[1:(length(GX)/2)]
    if(func_offset_varname!="" & person){
      return(message("offset specified as a ", paste(func_offset_varname),
          " function. Please create actual offset variable for regression build") )   
    
    }
    }
    
  }
  
  nai <- names(reg$coefficient[1])
  nointercept <- FALSE
  
  ## model with no intercept?  If Cox then has no intercept anyway.
  
  if(nai != "Intercept" & nai != "(Intercept)" & !cox ) {nointercept <- TRUE}

  
 ##----------------------------------------------------------------- 
 ## RED observation (person) to be added, indicated by non-NULL observation
 ##  ---------------------------------------------------------------
 if(person){
    
    if(nrow(observation)>1){
      message(paste("\"observation\" has >1 row. The last row provides plotted values"))
      observation <- observation[nrow(observation),]
      
          
    }

if( glm |  lm ) { 

      
      if(FIRSTRUN){
      actualdata <- model.frame(reg)

       vactual <- colnames( actualdata )
       adddata <- model.frame(reg, data=observation)
       vadded <-  colnames( adddata )
    
        if(offset){ vactual[which(vactual=="(offset)") ] <- offset_varname } 
    
          included <- is.element( vactual, vadded )
    
           vin <- vactual[which( !included )]
    
    
      if(length(vin)>0){

        message(paste("Variable(s)in glm model but not in \"observation\":",sep=""))
        return(message(paste(unique(vin),sep=",",collapse="  ")))
      
      }
      
      }  ##if(FIRSTRUN glm)
  
  both <- rbind(actualdata,adddata)
  
  ## note  reg$formula works here for glm, but no lm. reg$terms works for both!!
  MM <- model.matrix(reg$terms,data=both)
  ##  cannot do if NAs in observation, when MM has no rows 
  
  if(nrow(adddata) ==0) {
     message("NA values in  \"observation\".  \"observation\" xx")
  }

  newX <- MM[nrow(MM),]

## strip away the intercept term
if(!nointercept){newX <- newX[-1]}

}  ##if(glm |  lm) 


if( cox){
  
  if(FIRSTRUN){
  fvariable_names <- names(reg$assign)
  
   actualdata <- model.frame(reg)
  #remove Suv() element. Dunno how to deal with it!
  vactual <- colnames( actualdata )[-1]
  adddata <- model.frame(reg, data=observation)
  vadded <-  colnames( adddata )[-1]
  included <- is.element(vactual,vadded)
  vin <- vactual[which(!included)]
  
  if(length(vin)>0){
    message(paste("Variable(s)in Cox model but not in \"observation\":",sep="")) 
    return(message(paste(unique(vin),sep=",",collapse="  ")))
   
  }  
  
  
  adddata <- model.frame(reg, data=observation)

  ## doesnt like which() of survival time!! i.e. ignore 1st te
  FNA <- which(adddata[2:length(adddata)]== -999)

  if(length(FNA) > 0 )  { adddata[FNA + 1] <- NA} 
  
  }  ##if(FIRSTRUN cox


  both <- rbind(actualdata,adddata)
  MM <- model.matrix(reg,data=both)
  ##  use last row of combined data for plotting red points 
  newX <- MM[nrow(MM),]

  if(colnames(MM)[1]=="Intercept" | colnames(MM)[1]=="(Intercept)" ){newX <- newX[-1]}
  
  if(nrow(adddata) ==0) {
    message("NA values in \"observation\".  \"observation\" is yyy")

  }
 
}  #endif(cox)
 

  }  #endif(!is.null(observation) )

##-------------------------------------------------------- 
## need to get out actual variable corresponding to dummy 
## variables of factors. Do this irrespective of whether dummies=TRUE
## -------------------------------------------------------- 


if(!cox) {
   variable_names <- names(reg$model[2:length(reg$model)])
   
 
   ## first term of reg$model is the outcome Y variable
   ## other terms are X's whether no intercept model or not
   ## note: includes variable "(offset)" if offset in model
   ## reg$model doesn't include interactions if pressent
         }
else
        { variable_names <- names(reg$assign)
          
## for Cox  reg$assign also includes interactions. This is messy need to remove, leaving
## just variables used in the reegression.  Use presence of ":" and assume always appear last
          num <- length(variable_names)
          
          
          for (i in 1:(num)){
            
            ## obtain last factor that is NOT an interaction (i.e. assume interactions all contain : and tagged last  
            if(length(grep(":", variable_names[i] ) ) == 0 ) {no_int <- i}
          }
          
          variable_names <- variable_names[1:no_int]
          
}  ##else if(!cox
  
#---------------------------------------------------
  ols <- FALSE
  logistic <- FALSE
  poisson <- FALSE
  if(glm){
    ## has glm done an allowed family/link?
    ols <-      (reg$family[1]=="gaussian" & reg$family[2]=="identity")
    logistic <- (reg$family[1]=="binomial" & reg$family[2]=="logit")
    poisson <-  (reg$family[1]=="poisson"  & reg$family[2]=="log")
 
    if(!(ols | logistic | poisson)){
      return(message("regplot cannot do total score-to-outcome nomogram for this glm family/link"))
      }
      }
  
    
  if(!(cox | lm | ols | logistic | poisson )){
    return(message("regplot cannot do  total  nomogram for this regression"))}
  
  
  if(is.null( reg$formula )){ reg$formula <- formula(reg) }




  
npars <- length(reg$coefficients)
delta <- 0.05 + (0.2-0.05)*(npars-1)/(12-1)
   #browser()

  
 ns <- ""
 if(glm | lm | cox ){
   summry <- summary(reg)
   
 if(glm | lm) { Pval <- summry$coefficients[,4] }

 #P values in 5th col of  coxph class
 if(cox){ Pval <- summry$coefficients[,5] }
 
 ns <- ifelse( Pval<.05   ,"*","" )
 ns <- ifelse( Pval<.01  ,"**",ns ) 
 ns <- ifelse( Pval<.001,"***",ns ) 
 }
 
## extract Y variable from formula (nb also extracted above as yname from model)
  yvar <- as.character(reg$formula)[2]

if(FIRSTRUN ){
  
  
  
  if(nointercept){
    message("Regression model (not Cox) without an intercept")}
  
  if(cox ){ 
  
    x <- model.matrix( reg,data=reg$model )
  

  }
  else
  { x <- model.matrix( reg$formula, data=reg$model )}
  ## if replicate weighted model, expand out x matrix to unit records
  if(weighted){
  x <- x[rep(row.names(x),times =W), 1:ncol(x)]
  }
 ##---------------------------------------------------------------- 

  ##  large data?  Base distributions on subsample 5000
  ## if data huge, take subsample. Do so on first pass through


  if(nrow(x) > 5000){
    nobs <- nrow(x)
    
    set.seed(1234)
    
    subs <- sample(1:nrow(x), size=5000)
    

    x  <- x[ subs, , drop = FALSE]
 
 if(offset){ reg$offset <- reg$offset[subs]}
 if(distribution) { 
   message(paste("Distributions plotted based on random n=5000 subsample of",
    nobs))}
  }
 
} ##if(FIRSTRUN)
  

 ## need to mess about with parameters. coxph class has no intercept.

 if(cox){ 

   betas <- reg$coefficients[1:npars]
## fix up if npars=1, when x is a vector, coerce to matrix
   x <- as.matrix(x)
   X <- x[,1:npars, drop=FALSE]
   intercept <- 0

   
  }
 else   ##if(cox)
 {
 
  
 
  nai <- names(reg$coefficients[1])
  nointercept <- FALSE
  if(nai !="Intercept"& nai !="(Intercept)" ){
 ## model with no intercept
    nointercept <- TRUE
    intercept <- 0
    betas <- reg$coefficients[1:npars] 
    X <- x[,1:npars,drop = FALSE] 
  }
  else
  {
## model with intercept
   betas <- reg$coefficients[2:npars]
   ns <- ns[2:npars]
   intercept <- reg$coefficients[1]
  X <- x[,2:npars,drop = FALSE]
 npars <- npars -1

  }
 

##  regplot0.1 does not allow offset. Temporary code left in place
if(offset){
  ## offset included. Force in beta=1 value and column of offset values 
  
  
  npars <- npars +1
  betas <- c(betas,1)
  
   
  ## cant see how to pick up name of the offset variable. it is designated "(offset)"
  ##  in reg$model, so make this the coefficient name.
  names(betas)[npars] <- "(offset)"
  
 
  
  X <-  cbind(X,reg$offset)
  colnames(X)[npars] <- "(offset)"
  ns[npars] <- ""
  
  
  
}   ##if(offset

 }  ##if(cox


#------------------------------------------------------  




if(person  &  offset){
  
  ## for some reason appending using c() creates class(list). Need to unlist    
  
  ##   newX <- unlist(c(newX,observation[which(names(observation)=="(offset)")] ))
  newX <- unlist(c(newX,observation[offset_varname] ))
  
 
}


no_interaction <- npars

isinteraction <- vector(length=npars)
for (i in 1:(npars)){
   
  ## obtain last factor that is NOT an interaction (i.e. assume interactions all contain : and tagged last  
  if(length(grep(":", names( betas[i] ) ) ) == 0 ) {no_interaction <- i
  isinteraction[i] <- FALSE}
  else
  {isinteraction[i] <- TRUE
  }
}

num_vars <- length(variable_names) 

if(offset) { 
  
  ##remove "(offset)" variable name, to allow interaction terms to be
  ## added before the offset
 variable_names <-  variable_names[ which(variable_names !="(offset)")]
}

##  force in the  interaction terms (using beta names) as variable_names
  variable_names <- c(variable_names, names(betas[which(isinteraction==TRUE)]) )
  ## add the offset variable back in at the  top
  if(offset) {variable_names <- c(variable_names ,"(offset)") }

  num_vars <- length(variable_names) 
  if(weighted  & (glm | lm)){
    variable_names <- variable_names[1:(num_vars-1)]
    num_vars <- num_vars - 1 }
  #browser()
  
    if(FIRSTRUN & sum(isinteraction)>0 ){
    message("Model has interactions. ")
  
  }


## keep number of actual variable (ie in formula list in vector vact[]
vact <- vector(length=length(unlist(reg$xlevels)) +length(variable_names)  )
isadummy <- vector(length=length(unlist(reg$xlevels)) +length(variable_names)  )
isafactor <- vector(length=length(variable_names)  )
vact_names <- vector(length=length(unlist(reg$xlevels)) +length(variable_names)  )


i <- 0
minus1 <- 0
firstfactor <- NA
for(k in 1:num_vars){
  

  ## use reg$xlevels to determine whether is a factor variable 
  if(!is.null(unlist(reg$xlevels[variable_names[k]]) )) {
    
 ## need to frig around here for no intercept model. In this case
 ## the first factor var has a coefficient for ALL levels. 
    if(nointercept) {
    nlev <- length(unlist(reg$xlevels[variable_names[k]])) - minus1
    if(minus1==0) {firstfactor <- k}
    minus1 <- 1
  }
  else
  {nlev <- length(unlist(reg$xlevels[variable_names[k]])) - 1
  }
  
  
    isafactor[k] <- TRUE 
  }
  else
  {
    nlev <- 1
    isafactor[k] <- FALSE
  }
  
  for(j in 1:nlev){
    i <- i + 1
    vact[i] <- k
    
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
vact <- vact[1:i]
isadummy <- isadummy[1:i]


  S <- X
  
 
## -------------------------------------------------------- 
## extract scores in a loop 
## -------------------------------------------------------- 
  
 
 
  for (i in 1:(npars)){
  
 
    S[,i] <- X[,i]*betas[i]

  }


##------------------------------------------------------------------------------------
## center the  scores and data if center

  
  if(center){m <- colSums(X)/nrow(X)}
  else
 ## default align  to minimum values
  {m <- apply(X,2,min)}
  
  
  ## exclude factor variables from centering
     m[which(isadummy)] <- 0
  
  X <- scale(X,center=m,scale=FALSE)
  bm <- m*betas
  S <- scale(S,center=bm,scale=FALSE)
  intercept <- intercept +sum(bm, na.rm=TRUE)
  

 

## if dummies=FALSE, need to make rows of plot !dummies rather than dummy variables
## of factors. Do this by integrating over the dummy variables of factors.
if(!dummies) {
  

SS <- matrix(nrow=nrow(X),ncol=num_vars)
XX <- matrix(nrow=nrow(X),ncol=num_vars)
  
i <- 1
k <- 1
while (i <= npars) {
  
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

  
colnames(SS) <- variable_names

## replace  S and npars by collapsed values
npars <- num_vars
S <- SS
X <- XX



}  ##if(!dummies)
##------------------------------------------------------------



## set boundaries. Max and min of beta-scores or fraction thereof??
L <- min(S[,],na.rm=TRUE)
M <- max(S[,],na.rm=TRUE)

midpt <- (L+M)/2
diff <-  (M-L)
halfdiff <- diff/2

tot_score <-  rowSums(S, na.rm=TRUE) 


make_space <- 2
## graphics margin + 4 
par(mar=c(0,0,0,0) + 4,mai=c(0,0,0,0),bg = "#FFFFFF") 

if(FIRSTRUN) {plot.new()}

#------------------------------------------------------------
# main loop over predictors in the model, drawing graphic for each
# -------------------------------------------------------- 


sumcheck <- 0

nrows <- npars


factr <- vector(length=nrows)
isinteraction_row <-  vector(length=nrows)
## avoid ambigiuty by referring to each "row" of plot  

## set plot window 
plot(c(L-0.3*(M-L),M+0.1*(M-L)), c(0.5,npars+2 + make_space),col="white" ) 


aspect <- 1.4*(M-L)/(npars+2 + make_space -0.5) 

# #try scaling point axis 0,-100  (L,M)

lscale <- max(S[,],na.rm=TRUE)- min(S[,],na.rm=TRUE)





## vertical dashed grey zero line, if zero within the scale
## necessary ?? Probably not, Close off with FALSE
 

if(L<=0 & FALSE ){segments(0,make_space+0.5,0,nrows+1.5+make_space,col="#555555",lty=2)}

 for (i in 1:(nrows)){
   
   ipos <- i+0.5 +make_space
  ## unique doesnt sort, ordered in first appearance 
   
   uniq <- sort(unique(S[,i]))
   Xuniq <- sort(unique(X[,i]))
   luniq <- length(Xuniq)

  ## need to reverse direction for -ve beta
  B <- NULL
  if(!dummies){
    iv <- which(vact==i)
    if(length(iv)==1) {B <- betas[iv]}
  }
  else
  {B <- betas[i]}
  
  
  if( !is.null(B) ){ 
    
    ##browser()
  if(B < 0) {Xuniq <- sort(unique(X[,i]),decreasing=TRUE)
    }
  }
   
   
   ## deal with stat significance, choose largest if !dummies
   ## note: max("","*","**","***") returns "***"!
   
   if(!dummies){
     iv <- which(vact==i)
     ss  <- max(ns[iv])
     isint <- FALSE 
     if(length(iv)==1) {isint <- isinteraction[which(vact==i)]}
   }
   else
   {ss <- ns[i]
   isint <- isinteraction[i] }
   
   isinteraction_row[i] <- isint
   
   size <- 1
   ## smaller size font for interactions and split name in two. 
   
   v <- colnames(S)[i]
   
  
  if(isint){
    size <- 0.8 
  lengthv <- nchar(v)
   colon <- unlist (gregexpr(":",v) )
   v1 <- substring(v,1, colon-1)
  #replace : with an X
  v1 <- paste(v1," :",sep="") 
  v2 <- substring(v,colon+1, lengthv)
  text(x=L-0.3*(M-L),y= ipos+0.2,paste(v1),cex=size, adj=c(0,1),col="black")
  text(x=L-0.3*(M-L),y= ipos-0.04,paste(v2,ss,sep=""),cex=size, adj=c(0,1),col="black")
  
  
  ## if previous row NOT an interaction draw a dividing line
  
  if(!isinteraction_row[i-1]){  
    segments(L-0.3*(M-L),ipos-0.55,M,ipos-0.55,col="gray")}
    
  
   } 
   else
   {
   
   if(v=="(offset)") {v <- paste("offset(",offset_varname,")", sep="")
                      ##italicise font=3
   text(x=L-0.3*(M-L),y= ipos+0.25,paste(v,ss,sep=""),cex=size, font=3, adj=c(0,1))}
   
   else
   {
   text(x=L-0.3*(M-L),y= ipos+0.25,paste(v,ss,sep=""),cex=size, adj=c(0,1))
   }
   }
   
   
  ## plot box/bean/violin plot. If plot="none" do same to output just a scale
  ## n_do_a_box number of unique values to do a boxplot
   n_do_a_box <- ndisc
  ## force so that  will not do box plot if(!dummies & isafactor[i]) 
   if(!dummies & isafactor[i]){ n_do_a_box  <- 1000}  
   
   
   if(luniq >  n_do_a_box ) {
     factr[i] <- FALSE
     if(distribution){
       
     if(cplot=="violin"){
       
       bandwidth <- 0.01*(max(S[,i])-min(S[,i]))
   box <-  vioplot(S[,i], add=TRUE,at=ipos+0.3,horizontal=TRUE,wex=0.5, border="blue", 
              h=bandwidth, colMed="black", pchMed=21, col=bvcol)
            }

   if(cplot=="bean" ){
 
  box <-  beanplot(S[,i], add=TRUE,at=ipos+0.3,horizontal=TRUE, maxwidth=0.5,
                   what = c(0,1,0,1),  bw="nrd0", log="", overallline="median",beanlinewd=0, 
      kernel="epanechnikov",  axes=FALSE, border="blue", col=bvcol,method="overplot",
      ll=0.1 , frame.plot=FALSE,  pars = list(boxwex = 0.8, staplewex = 0.5, 
                                              outwex = 0.5,cex=0.5))  
  }
  
  if(cplot=="density") {
  dens <- density(S[,i])
  maxden <- max(dens$y)
  dens$y <- 0.5* dens$y/maxden +ipos
  
  zero <- c(rep(ipos,times=length(dens$x)) ) 
  
  ## do fill-in with handy bit of code grabbed from somewhere! 
  polygon(c(dens$x, rev(dens$x)), c(dens$y, zero),
          col = bvcol, border = NA)
    lines(dens$x,dens$y,col="blue")
  segments(min(dens$x),ipos,max(dens$x),ipos,col="black")
    }

    if(cplot=="boxplot"){
  ## draw box plot for !dummies with > ndisc distinct values

    box <-  boxplot(S[,i], add=TRUE,at=ipos+0.3,horizontal=TRUE,  border="blue",col=bvcol,
             pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.5,cex=0.5))
         }
   
  
  ## spikes, frequencies of unique values
  
  if(cplot=="spikes"){
    ## continuous frequency table - will be ordered by  increasing S 
    
    frq <- as.numeric(table(S[,i]))
    
    maxfreq <- max(frq)
    
      
    dy <- 0.7*frq/maxfreq 
    dx <- dy * aspect
    xleft <-  uniq - dx
    xright <- uniq + dx
    ybottom <- c(rep(ipos,times=luniq)) 
    ytop <- c(rep(ipos,times=luniq)) + dy
    
  ##spikes of continuos  
    for(nspike in 1:length(uniq)){
     
      segments( uniq[nspike], ybottom[nspike], uniq[nspike], ytop[nspike], lwd=3,col="blue")
    }
    ## also add an axis 
    segments( min(uniq),ipos,max(uniq),ipos) 
    }

     } ##if(distribution) 
nmax = round(min(10, 1+10*(max(S[,i]) -min(S[,i]))/lscale ),0)


## now add the variable scale ruler
##  extract beta value for this numeric variable. 
if(!dummies){
  Beta <- betas[which(vact==i)]
}
else
{Beta <- betas[i]
}



  if(!dummies){  mean <- m[which(vact==i)] }
 
  else
  {mean <- m[i]}
    
  ticks     <- pretty(X[,i]+mean,n=nmax) 
  ticks_pos <- ticks*Beta -mean*Beta 
 

## add the scale 

shift <- 0.3
if(!distribution | cplot=="density") {shift <- 0}
shift <- 0
segments(min(ticks_pos),ipos+shift,max(ticks_pos),ipos+shift,col="black",)
segments(ticks_pos,ipos+shift,ticks_pos,ipos+shift -0.05,col="black",)
subticks(ticks_pos,ticks_pos,ipos+shift)
text(x=ticks_pos,y= ipos+shift - 0.15,paste(ticks),cex=0.65,col="black",)


     } ##if(luniq >  n_do_a_box)
##------------------------------------------------------------------------
   else  ##if(luniq >ndisc )
     
 ## draw scaled boxes for each discrete numeric or factor data
   { 
     
    factr[i] <- TRUE 
## frequency table - will be ordered by  increasing S 

     frq <- as.numeric(table(S[,i]))
   
      tot <- sum(frq)  
    

    if(!dummies){
      
      if(isafactor[i])
    
    {
      ## is a factor, so get corresponding beta values of each level
      ## in order of the model
      Xuniq <- unlist(reg$xlevels[variable_names[i]])
      
 #browser()
 
 
      B  <- betas[which(vact_names==variable_names[i])]
      
      
      
      ## augment  B with baseline value zero
      if(nointercept & firstfactor==i){
       
          beta_values <- B}
      
      else
       
       {
      beta_values <- c(0,B)
      ## artificial beta=0 for ref category 
      ## pick up baseline category name
      names(beta_values[1]) <- names(Xuniq[1])
      }
      
      
      ##  Need to order with 0 added to bring in line with order of frq table
      obetas <- order(beta_values)
     
     
     ##browser()
      Xuniq <- Xuniq[obetas]
      
      luniq <- length(Xuniq)
      
      ##val <- Xuniq[o]
      
      val <- Xuniq
      
  ##  but uniq also needs to be ordered by  magnitude and X uniq too    
      
      ofact <- order(uniq)
      uniq <- uniq[ofact]
      
    }
    
    else
    {
      
 ## !dummies requested but not a factor, a discrete numeric variable
     
    
     val <- signif(Xuniq+m[which(vact==i)], 2)
     
    } 
      
    }  ##if(!dummies)
    
    if(dummies){
     
    
        val <- signif(Xuniq+m[i], 2)
     
      
    } 
   
  
o <- order(frq,decreasing=TRUE)
frq <- frq[o]
uniq  <- uniq[o]
val <- val[o]

    if(distribution){
      ## draw squares so that  lagests is first, smaller superimposed
      if(dplot=="spikes"){
              dy <- 0.7*frq/tot 
              dx <- 0.001* diff
              xleft <-  uniq - dx
              xright <- uniq + dx
              ybottom <- c(rep(ipos,times=luniq)) 
              ytop <- c(rep(ipos,times=luniq)) + dy
             ##browser()
              
      }
      else
      {
        
   dy <- 0.35*sqrt(frq/tot) 
   dx <- dy * aspect
   xleft <-  uniq - dx
   xright <- uniq + dx
   ybottom <- c(rep(ipos,times=luniq)) - dy
   ytop <- c(rep(ipos,times=luniq)) + dy
   
   
      }
    }
   else
    {dx <- 0.01
    dy <- 0.01
    xleft <-  uniq - dx
    xright <- uniq + dx
    ybottom <- c(rep(ipos,times=luniq)) - dy
    ytop <- c(rep(ipos,times=luniq)) + dy
   }
   
    
   ## is possible that there is no such combination in interactions
   if(length(uniq) >0) {
   

            ## also add an axis ??
       segments( min(uniq),ipos,max(uniq),ipos,col="gray") 
      
   if(distribution){
     
       colfill <- sqcol
     
     if(dplot=="spikes"){
       for(nspike in 1:length(uniq)){
        
       segments( uniq[nspike], ybottom[nspike], uniq[nspike], ytop[nspike], lwd=3,col="blue")
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
     text(uniq,ipos-dy-0.05,paste(val),cex=0.7,col="black")
     }
     
   }  ##if(distribution) 
   
    ## add labels to squares,  
   else  ##if(distribution) 
     
   {
  
   
   text(uniq,ipos-dy-0.15,paste(val),cex=0.7,col="black")
   size <- 0.3
   
   points(0.5*(xleft+xright),(ybottom+ytop)*0.5,cex=1.3, col="black")
   
   }
   ##  put in point at center of box
 
 
 
}  ##if(length(uniq) >0)
   
}  ##else  of if(luniq

if(person){
    
  if (!is.na(newX[i]) ){
#  add in  RED points  of observation 
    
    if(dummies){
      
        addX <- (newX[i]-m[i])*betas[i]
        
        if(luniq >  n_do_a_box){ text(addX,ipos+0.15,paste(signif(newX[i],3)),cex=0.7,col=obscol,font=2)}
        
        
    }  ##if(!!dummies)
    
    else
      
    {
      
      dummys <-  which(vact_names==variable_names[i])
      
      
      #browser()
      if(isafactor[i]){
        
           
          
          
          value_one <- which(newX[dummys] ==1)
          
          
          
          
          levels <- unlist(reg$xlevels[variable_names[i]])
          
          
          
          
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
            
            variable_names[i]
            addX <-  beta_values[names(value_one)]           
          }  
          j <- which(Xuniq[o]==level)
          posx <- uniq[j]
         
          
          ##write label text in bold (font=2) 
       ##   text(posx,ipos+0.15,paste(level),cex=0.7,col=obscol,font=2)
         #browser()
      }  ##if(isafactor
    else
    {
      
        #!dummies & !isafactor
      addX <-  betas[dummys]*(newX[dummys]-m[dummys])
     

      ## output  actual value in red. Necessary only for continuous
      if(luniq >  n_do_a_box){ text(addX,ipos+0.15,paste(signif(newX[dummys],3)),cex=0.7,col=obscol,font=2)}
          
    }
    
    } ##if(!dummies  
        
## add RED dot point at observation value

   #browser()
    cex <- 1
    col <- obscol
    if(isint){
      cex <- 0.7
      col <- obscol}
  points(addX,ipos, col=col,cex=cex, pch=21, bg=obscol)

## also add red point to top beta X scale
  npos <- nrows+1 + make_space
  
  points(addX,npos+0.4, col=col,cex=cex, pch=21, bg=obscol)
#join by vertical line??
if(droplines){
  segments(addX,ipos+0.02,addX,npos+0.38, col="pink") 
}
if(!is.na(addX)){

if(pointscale){
  
  sumcheck <- sumcheck + 100*( addX-L )/( M-L )
}
else
{ sumcheck <- sumcheck +addX}
}
} ##if (!is.na(newX[i]) )

} ##if(person)


 }
## end of main loop for (i in 1:(nrows))
##-----------------------------------------------------

## note if interactions?

if(isint){text(x=L-0.25*(M-L),y= ipos+0.6,paste("Interactions:"),cex=0.9,col="#888888")} 

#-----------------------------------------------------
# add beta X  or points scale  contributions at top
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
npos <- nrows+1 + make_space

##npos <- npars+1 + make_space
segments(min(tickp),npos+0.4,max(tickp),npos+0.4)
segments(tickp,npos+0.4,tickp,npos+0.33)
cexv <- 0.9
if(npars >15) {cexv <- 0.7}
text(x=tickp,y= npos+0.4+delta,paste(signif(tickval,3)),cex=cexv) 

# subticks?

subticks(tickp, tickp,npos+0.4)



if(pointscale) {ex <- "          Points"}
  else
{ex <- expression(italic(paste(beta,     "X contributions")))}

text((L+M)/2,npos +0.4 -delta , ex ,cex=0.8, font=3) 

#write titleheading above  upper scale 
titlepos <- min(npos+0.4 + 2.5*delta,npars+2+ make_space)

text((L+M)/2,titlepos, title ,cex=1.05 ,font=4) 


  if(logistic){    
## trim ridiculous logit scores >|8|
##  tot_score <- tot_score[which(abs(tot_score + intercept) < 8)]
}

# -------------------------------------------------------- 
# add  nomogram at the bottom 
# -------------------------------------------------------- 

  ##  add  Total score axis 
  ## font=4 is bold italic =3 is italic
  if(pointscale){
  text(L, make_space+0.5,paste("Total-points-to-outcome nomogram:"),font=4)}
  else
  {text(L, make_space+0.5,paste("Total-score-to-outcome nomogram:"),font=4)}
  
  
## now force the position of the box plot to coincide with L, M at the end points

Max_tot_score <- max(tot_score)
Min_tot_score <- min(tot_score)
Range <- Max_tot_score-Min_tot_score



scale_score <- (M-L)* (tot_score -Min_tot_score)/Range + L


npos <- 1

##--------------------------------------------------------------
##------------------Distribution of total score-----------------
if(distribution){
  
  score_uniq <- unique(scale_score)
  luniq <- length(score_uniq)
  
## allow possibility for boxes of total score if few distinct values (< ndisc)
  if(dplot=="boxes" & luniq < ndisc){
    frq <- as.numeric(table(scale_score))
    
    tot <- sum(frq)
    
    dy <- 0.35*sqrt(frq/tot) 
    #dx <- dy*scale
    dx <- dy * aspect
    xleft <-  score_uniq - dx
    xright <- score_uniq + dx
    ybottom <- c(rep(npos+1,times=luniq)) - dy
    ytop <- c(rep(npos+1,times=luniq)) + dy
    
    
    
    rect(xleft, ybottom, xleft+2*dx, ytop,border="blue",col="white")
    
    
    points(0.5*(xleft+xright),(ybottom+ytop)*0.5,cex=0.3, col="black")
     
    
  } ##if(dplot=="boxes" & length(score_uniq) < ndisc)  
  
  
  else
    
  {
  
    if(cplot=="violin"){
    
    bandwidth <- 0.01*(max(scale_score)-min(scale_score))
    box <- vioplot(scale_score, add=TRUE,at=npos+1,horizontal=TRUE,  wex=0.7, h=bandwidth, 
                   colMed="black", pchMed=21, border="blue",col=bvcol)
  }
   if(cplot=="bean") {
     

     box <-  beanplot(scale_score, add=TRUE,log="", at=npos+1,horizontal=TRUE,wex=0.8, bw="nrd0", 
                    what = c(0,1,0,1),maxwidth=0.5,ll=0.1,overallline="median",beanlinewd=0, method="overplot",
                     kernel="epanechnikov",  axes=FALSE, border="blue", col=bvcol)
             pars = list(boxwex = 0.4, staplewex = 0.5, outwex = 0.5,cex=0.5) }                  
  if(cplot=="boxplot"){

  box <-  boxplot(scale_score, add=TRUE,at=npos+1,horizontal=TRUE,  border="blue",col=bvcol,
                  pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.5,cex=0.5))
      
  }
  
  if(cplot=="density") {
    ## testing doing kernel density
    dens <- density(scale_score)
    dens$y <- dens$y <- 0.6* dens$y/max(dens$y) +npos+0.5
 ##   lines(dens$x,dens$y,col="blue")
    zero <- c(rep(npos+0.5,times=length(dens$x)) ) 
    
  ## fill in with handy bit of code  
    polygon(c(dens$x, rev(dens$x)), c(dens$y, zero),
            col = bvcol, border = NA)
    lines(dens$x,dens$y,col="blue")
    
   
  }
  
  
  if(cplot=="spikes"){
    ## frequency table - will be ordered by  increasing S 
    
    frq <- as.numeric(table(scale_score))
    score_uniq <- sort(  unique(scale_score)  )
    
    maxfreq <- max(frq)
    
    
    dy <- 0.7*frq/maxfreq 
    dx <- 0.001* diff
    
    
    #dx <- dy*scale
    dx <- dy * aspect
    xleft <-  score_uniq - dx
    xright <- score_uniq + dx
    luniq <- length(score_uniq)
    ybottom <- c(rep(npos+0.5,times=luniq)) 
    ytop <- c(rep(npos+0.5,times=luniq)) + dy
    
    ## spikes of continuous  
    for(nspike in 1:luniq){
     
      segments( score_uniq[nspike], ybottom[nspike], score_uniq[nspike], ytop[nspike], lwd=3,col="blue")
    }
  
  }  ##if(cplot=="spikes")
 
  }
    
  
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
  
  ## convert to points, but need to account for there being sum of npars parameters
  tickv <- 100*(tickv -L*npars)/(M-L)
  tickv <- pretty(tickv,n=8)
  ## convert back into actual total scores
  pts  <- ((M-L)*tickv)/100 +L*npars
  tickp <- (M-L)* (pts -Min_tot_score)/Range + L 
} 


segments(min(tickp),npos+0.5,max(tickp),npos+0.5)
segments(tickp,npos+0.5,tickp,npos+0.43)
  text(x=tickp,y= npos+0.5 -1.*delta ,paste(tickv),cex=0.7)

# subticks

subticks(tickp, tickv,npos+0.5)

tickvalue <- tickv

if(pointscale) {ex <- "Total points"}
else
{ex <- expression(italic(paste("Total ", beta,"X ")))}

text(L-0.2*(M-L),npos+0.5+delta, ex,cex=0.8 ,font=3 ) 
# draw box plot for variables with > 5 distinct values



if(person){

 
  
  if(length(newX) != length(betas)){
    
 

    message(" \"observation\" inconsistency with regression data. There maybe NAs, not allowed")}

  else
    
    {
   totS <- sum((newX-m) *betas, na.rm=TRUE) 

    
    }
    
    totSpos <- (M-L)* (totS -Min_tot_score)/Range + L
 ## drawing pointer arrow scale in red
 ## filled diamond pch=23
    sz <- 1.3
    if(npars > 15) {sz <- 1}
  
  points(totSpos,npos+0.5, col=obscol,cex=sz, pch=23, bg=obscol)
 
  ## draw downward arrow cursor  at just above lower scale at 0.5 +0.35 
  segments(totSpos,npos+0.5, totSpos, 0.5 + 0.45, col=obscol,cex=1.)
  ## arrow tip is at 0.5 + 0.35, because this is the position of lower scale 
  points(totSpos, 0.5 + 0.45 , col=obscol,cex=sz, pch=25, bg=obscol)
 
  
   ## total in red
  
  exval <- paste(signif(totS,3),sep="")
  
  if(pointscale) {
    ## need to frig around to get the outputted points value
    
    exval <- paste(signif(sumcheck,3), sep="")
    

  }
  ## output total in red of red diamond  on the total score axis   
 
 text(totSpos,npos+.7, exval,cex=0.9 ,font=3,col=obscol ) 
 
}
  
## (font=3 doesn't seem to work on an expression)
## add score-to-probability nomograms at the bottom 

# ------logistic  nomogram-----------------------------------------  

if(logistic){
npos <- 0.5
E <- exp(tot_score+intercept)
prettyprobs <- pretty( E/(1+E), n=10)


# fudge if end points are 0 and 1.


nticks <- length(prettyprobs)

if(prettyprobs[1] < 0.005 ) {prettyprobs[1]  <- 0.005}

if(prettyprobs[nticks] > 0.995 ){prettyprobs[nticks] <- 0.995 }
# take account of intercept, by subtracting from logit

prettyscores <- log(prettyprobs/(1-prettyprobs)) -intercept

tickpos <- (M-L)* (prettyscores  - Min_tot_score)/Range + L
  
limits <- which(tickpos > L-0.3*(M-L) & tickpos < M+0.1*(M-L) )
tickpos <-  tickpos[limits]
prettyprobs  <-  prettyprobs[limits]


segments(min(tickpos),npos+0.35,max(tickpos),npos+0.35)
segments(tickpos,npos+0.35,tickpos,npos+0.27)

## use sieve function to avoid text overlay on axis

sievePr <- sieve(tickpos,prettyprobs, 0.04*(M-L))
#logistic
delta <- 0.05 + (0.2-0.05)*(npars-1)/(12-1)
## line at npos +0.35, text a little above 
text(x=unlist(sievePr[1]),y= npos+0.35 + delta,paste(signif(unlist(sievePr[2]),3)),cex=0.7)


text((max(tickpos)+min(tickpos))/2,npos+0.0, paste("Pr(",yvar,")") ,cex=0.8 ,font=3)
##browser()
subticks( tickpos,prettyprobs,npos+0.35,func="expit", decide=TRUE)

if(person){
  
  pred_mean <- exp(totS+intercept)/(1 + exp(totS+intercept))
  text(totSpos,npos+0.25, paste(signif(pred_mean,3)),cex=0.9,font=3,col=obscol  ) 
  

}


}
# ------lm  nomogram-----------------------------------------  

if(lm | ols ){
  npos <- 0.5
  tickmean <- tot_score+intercept
  prettyscores <- pretty(tickmean)
  
  
  ## tick positions  correspond to scale of the boxplot and its position
  ## after subtracting out the  intercept
  
  tickpos <- (M-L)* (prettyscores -intercept - Min_tot_score)/Range + L
    
  
  
  ## draw line to fit L-M boundary
   
  segments(min(tickpos),npos+0.35,max(tickpos),npos+0.35)

  segments(tickpos,npos+0.35,tickpos,npos+0.27)
  
 #logistic
delta <- 0.05 + (0.2-0.05)*(npars-1)/(12-1)
## line at npos +0.35, text a little above 
  
  text(x=tickpos,y= npos+0.35 +delta,paste(signif(prettyscores,3)),cex=0.7)                                     
  text((max(tickpos)+min(tickpos))/2,npos+0.0, paste("mean",yvar),cex=0.8 ,font=3 ) 
  ##subticks(tickpos, tickmean,npos+0.45)
  subticks(tickpos, prettyscores,npos+0.35)
 
  if(person){
   
    pred_mean <- totS+intercept
    text(totSpos,npos+0.25, paste(signif(pred_mean,3)),cex=0.9,font=3,col=obscol  ) 
    

  }
  
  }

# ------poisson  nomogram-----------------------------------------  

if(poisson){
  npos <- 0.5
  
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
  
  
  ##browser()
  
  limits <- which(tickpos > L-0.3*(M-L) & tickpos < M+0.1*(M-L) )
  
  tickpos <-  tickpos[limits]  
  prettymean <- prettymean[limits]
  
## also limit by upper value of counts 
   mxcount <- exp(intercept+max(tot_score))
   tickpos <-  tickpos[which(prettymean< mxcount)]  
   prettymean <- prettymean[which(prettymean< mxcount)]
  segments(min(tickpos),npos+0.35,max(tickpos),npos+0.35)
  segments(tickpos,npos+0.35,tickpos,npos+0.27)


## use sieve function to aviod text overlay on axis
  
sievePr <- sieve(tickpos,prettymean, 0.05*(M-L))

## logistic
delta <- 0.05 + (0.2-0.05)*(npars-1)/(12-1)
## line at npos +0.35, text a little above 
text(x=unlist(sievePr[1]),y= npos+0.35 +delta,paste(signif(unlist(sievePr[2]),3)),cex=0.7)

 ## text(x=tickpos,y= npos+0.6,paste(signif(prettymean,3)),cex=0.7)  
  
  subticks( tickpos,prettymean,npos+0.35,func="exp")
  text((max(tickpos)+min(tickpos))/2,npos+.05, paste("Poisson mean",yvar),cex=0.8,font=3  ) 

   if(person){
     
     pred_mean <- exp(totS+intercept)
     text(totSpos,npos+0.25, paste(signif(pred_mean,3)),cex=0.9,font=3,col=obscol  ) 

   }

}

# ------cox nomogram-----------------------------------------  

if(cox){
  ## use basehaz() to get hazards at X=0. Only do FIRSTRUN as it
  ## is time consuming
  if(FIRSTRUN) {
  h <- basehaz(reg, centered=FALSE)
  }
  ## possibility that hazard not computable - NAs in betas? 
  if(!is.na( h$hazard[1] ) ){
    
 
    
   sumexp <- exp( sum(betas*m,na.rm=TRUE) )
    hz <- h$hazard*sumexp 
 
  omit_nomo <- FALSE
 ## no specified failtime  - use median 

  if(is.na(tcut)){
   tmed <-  quantile(  h$time, probs=c(0.49,0.5,0.51), na.rm=FALSE )
   h5 <- hz[which(h$time<tmed[3] & h$time>tmed[1])]
   tcut <- tmed[2]
  }
  
  else
  { 
## user specified failtime. Select hazard "close" to tcut
    st <- sort(h$time)
    ibelow <- max(which(st<=tcut ))         
    h5 <- hz[c(ibelow, ibelow+1)]
  
    if(tcut > max(st) | tcut < min(st)){
    
     
     return( message("failtime parameter mis-specified. Out of range of observed") )
     omit_nomo <- TRUE
    }
  }
 
 if(!omit_nomo){

 ##  take default mean-time unit baseline 
 s5 <- mean(exp(-h5))
 ## get pretty survival scores
 
 prettyprobs <- pretty( s5 ^ exp(tot_score), n=8)
 ## fudge add extra fill in on survival scale if  starts at =0.9
 nticks <- length(prettyprobs)
 
 
 if(prettyprobs[1] <0.005 ) {prettyprobs[1]  <- 0.005}
 
 if(prettyprobs[nticks]>0.995 ) {prettyprobs[nticks] <- 0.995 }
 
 
 
 ## fix position of lower nomogram scale 
 npos <- 0.5
 
 
 ## use S(t)^exp(BX)=1-P to get scores corresponding to probabilities
 exp_score <- log(prettyprobs)/log(s5)
 
 
 ## corresponding  betaX scores are 
 prettyscore <- log(exp_score)
  
 
 tickpos <- (M-L)* (prettyscore  - Min_tot_score)/Range + L
 
 
 limits <- which(tickpos > L-0.3*(M-L) & tickpos < M+0.1*(M-L) )
 tickpos <-  tickpos[limits]
 
 
 prettyprobs  <-  prettyprobs[limits]

 
 ## main tick scale
 segments(min(tickpos),npos+0.35,max(tickpos),npos+0.35)
 segments(tickpos,npos+0.35,tickpos,npos+0.27)
 
  
 subticks(tickpos,prettyprobs,npos+0.35,func="cox",s5,decide=TRUE)
 

## use sieve function to avoid text overlay on axis

sievePr <- sieve(tickpos,prettyprobs, 0.05*(M-L))

Xpos <- unlist(sievePr[1])
Pr <-   unlist(sievePr[2])

 ## put values on the axis 

delta <- 0.05 + (0.2-0.05)*(npars-1)/(12-1)
 if(fail) {
 text(x=Xpos,y= npos+ 0.35 +delta,paste(signif(1-Pr,3)),cex=0.7)
 }
 else
 {
  text(x=Xpos,y= npos+0.35 + delta,paste(signif(Pr,3)),cex=0.7)}
 
 ## extract time variable name  from  "Surv(time,fail)" structure
 stripped  <- sub(pattern="Surv\\(", replacement="", x=yvar )
 stripped  <- sub(pattern="\\)",     replacement="",x=stripped)
 stripped  <- sub(pattern=" ",       replacement="",x=stripped)
 
 outcome_names <- strsplit(stripped,",")
 yvar <- unlist(outcome_names)[1]
 deadvar <-  unlist(outcome_names)[2] 


newY <- observation[which(colnames(observation)==yvar)]
deadval <- observation[which(colnames(observation)==deadvar)]
# 
 if(fail){ chr <- "<" } else { chr <- ">" }
   
  text((max(tickpos)+min(tickpos))/2, npos+0.0, paste("Pr(",yvar,chr, signif(tcut,4),")" ), cex=0.8,font=3) 


 }


if(person){
  ##S(t)^exp(BX)=1-P  
  pred_mean <- s5^exp(totS)
  if(fail) { pred_mean <- 1-pred_mean }
  text(totSpos,npos+.25, paste(signif(pred_mean,3)),cex=0.9,font=3,col=obscol  ) 
  
  if(FIRSTRUN){
  ## pick up observation status indicator of "fail" from last row of appended data
  ## column 1 has 2 elements time, status. So pick 2nd. 
  died <- both[nrow(both),1][2]
  cdied <- "(no fail)"
  if(died == 1) {cdied <- "(failed)"}
  ## suppress actual outcome
  ## text( totSpos,npos+0.15, 
  ##      paste(paste(yvar,"=",signif(newY,3) ,sep=""),cdied ),
  ##      cex=0.8,font=3,col=obscol)
}
}
} ##if(!is.na(h$hazard[1])

else
{message("Hazard and survival  not computable. Possibly due to NA betas")
person <- FALSE
}

}
# end of ------cox nomogram-----------------------------------------  

 ##   ANIMATE <- !is.null(observation)

   if(ANIMATE){
   if(FIRSTRUN) message('Click is expected.  Click ESc or Press Esc to quit')
     
     FIRSTRUN <- FALSE
      
     gap <- (npars+2 + make_space)/100
     ## position just down from top
  
     dialogpos <- npars+2+ make_space + 1.5*gap

     xgap <- (M-L)/10
     rect( L, dialogpos-gap, L+xgap*8.5, dialogpos+gap,border="blue",col="white")
     
     text( L+1*xgap,dialogpos,  paste("none"),cex=0.7,col="blue")
     text( L+2*xgap,dialogpos,  paste("boxplot"),cex=0.7,col="blue")
     text( L+3*xgap,dialogpos,  paste("bean"),cex=0.7,col="blue")
     text( L+4*xgap,dialogpos,  paste("violin"),cex=0.7,col="blue")
     text( L+5*xgap,dialogpos,  paste("density"),cex=0.7,col="blue")
     text( L+6*xgap,dialogpos,  paste("spikes"),cex=0.7,col="blue")
     text( L+7*xgap,dialogpos,  paste("boxes"),cex=0.7,col="blue")
     text( L+8*xgap,dialogpos,  paste("Esc"),cex=0.7,col="blue")                          
     
    ## click on graphic to locate one (n=1) x,y coordinate
    
    XY <- locator(n=1)
    
    row <- floor(XY$y -  make_space ) 
   
  esc <- FALSE

  
  #Press Esc to escape
   if(length(row)==0){
     esc <- TRUE
                    }
  
  else

  {
   
 ## plot_select <- floor( (dialogpos-XY$y)/gap + 1.5)
    plot_select <- floor((XY$x-L)/xgap +0.5)
    
  
 ## click on dialog area
 
  if(XY$y > dialogpos - gap & plot_select <=8 & plot_select >= 1 ){  
  if(plot_select == 1) cplot <- "none"
  if(plot_select == 2) cplot <- "boxplot"
  if(plot_select == 3) cplot <- "bean"
  if(plot_select == 4) cplot <- "violin"
  if(plot_select == 5) cplot <- "density"
  if(plot_select == 6) {dplot <- "spikes"
   cplot <- "spikes"                     
 }
  
  if(plot_select == 7) {dplot <- "boxes"
                        if(!distribution ) {cplot <- "boxplot"}
  }
  if(plot_select == 8)   {  esc <- TRUE }
     
  }
  
  else
  ## click within rows of plot if person 
    
   {
   if(!person){
     message("observation is NULL. Specify an observation argument")
   }
    if(person){ 
      
     if(row >= 1 &  row <= npars ){
    
     #browser()
    if(!isinteraction_row[row]){
  
   
  
#------------------------------------------------------------
  if(!dummies){
    
    
    
    if(isafactor[row]) {
    
     
    
    Xuniq <- unlist(reg$xlevels[variable_names[row]])
    
    if(nointercept & firstfactor == row){
      beta_values <- betas[which(vact_names==variable_names[row])]
      
      eps <- abs(XY$x-beta_values )
     
      
    }
    else
    {
    
    beta_values <- betas[which(vact_names==variable_names[row])]
    ##  Need to order with 0 added
     
    ## find closest point clicked
    eps <- abs(XY$x-c(0,beta_values) )
    
    
    } 
    oeps <- order(eps)
    
  ## update adddata with clicked value   
    
        adddata[variable_names[row] ]  <- Xuniq[oeps[1]]
    
    
    }
    
    else  ##if(isafactor[row])
    
      { nfact <- which(vact==row)
      value <- XY$x/betas[nfact]
   ##update adddata with clicked value   
       var <- variable_names[row]
      value <- value +m[nfact]
      adddata[var] <- value
       
  ##    
   } 
    
    ##browser()
}

else
 ##if(!dummies)................... 
{
  ##do for dummies=TRUE
  
    value <- XY$x/betas[row]

    value <- value +m[row]
    
    
    
    
    if(factr[row]){
      value <- round(value)}
      var <- colnames(S)[row]
  
if(isadummy[row]){
  var_name <- variable_names[vact[row]]
  level <- substring(var,nchar(var_name)+1,nchar(var))
  
  ## if clicked on 0, the reference level has been chosen, which is first item of unlisted xlevels
  
 if(value == 0) level <- unlist(reg$xlevels[var_name])[1]
 ## observation[var_name] <- level
 adddata[var_name] <- level
}
else
{
  
  
  if(var =="(offset)") {var <- offset_varname}
adddata[var] <- value
}


} 

}  ##if(!interaction

else
  
{
  
  message("Cannot click on an interaction") 
  

}

}  ##if(row>=1 & <=npars
else
  
{
  
  esc <- TRUE
  
  
}


} ##if(person)

}  ##else of XY$x < L & plot_select <=5 & plot_select >= 1


}  ##if(length(row)==0 i.e. pressed Esc

  

if(esc){
## temp suppress this for plot with header menu shown  
## erase dialog box with a white-over box slightly bigger

  rect( L-0.1*xgap, dialogpos-gap, L+10*xgap, dialogpos+gap,border="white",col="white")
  text((L+M)/2,titlepos, title ,cex=1.05 ,font=4) 
return("Exiting regplot")
}

   }  ##if(ANIMATE)


}  ##WHILE(ANIMATE)
}  ##end of regplot function
##------------------------------------------------------------------------

subticks <- function(maintickpos, maintickvalue,npos,func=NA,s5, decide=FALSE){
l <- length(maintickvalue)
## add 4 subtick marks between major ones, unless decide = TRUE
n4 <-  c(1,2,3,4)
mtv <- maintickvalue
if( !is.na(func) ) {
  if(func == "exp")  {mtv <- log( maintickvalue )}
  if(func == "expit"){mtv <- log( maintickvalue / (1-maintickvalue) )}
  if(func == "cox")  {mtv <- log( log(maintickvalue) / log(s5) )}
                 
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
    stv <- log( log(stv) / log(s5) )}        
  }
  
  
    R <- (maintickpos[i+1] - maintickpos[i])/(mtv[i+1]-mtv[i])
    subtickpos <-  R*(stv - mtv[i]) + maintickpos[i] 

  segments( subtickpos,npos,subtickpos,npos-0.035 )
}
  
} 
##-----------------------------------------------------------------------------------
## Utility functions to  extract variable from  function e.g  "age" from "log(age)"
##  instr grabbed from a web page and adapted  find positions of "(" and ")"
instr <- function(str1, str2, startpos=1, n=1){
  
  aa <- unlist(strsplit(substring(str1,startpos),str2))
 
  if(length(aa) < n ) return(0);
  return(sum(nchar(aa[1:n])) + startpos+(n-1)*nchar(str2) )
}

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
  

  if( right > left  & left >1) {str[i] <- substr(str[i], left+1, right-1)
  func[i] <-  substring(in_str[i],1,left-1) }
  
  }
  
  outgetX <- c(str,func)
  
  return(outgetX)
}


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
      