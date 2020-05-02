#Transfrom a list of variable name denoted selection results to 0-1 matrix result
f01 <- function(object, full.var, p){

  matrix.01 <- data.frame(matrix(0, length(object), p))
  colnames(matrix.01) <- full.var
  for(i in 1:length(object))
  {
    matrix.01[i, (full.var %in% object[[i]])] <- 1
  }
  return(matrix.01)
}

#Calculate the Model confidences bounds and their corresponding freq rate
CI <- function(var.list, var.matrix, p, B)
{
  full.var <- colnames(var.matrix)
  colsum <- apply(var.matrix, 2, sum)
  order <- order(colsum, decreasing = T)
  freq <- vector(length = p+1)
  freq[1] <- 0
  lower <- vector(mode="list", length = p+1)
  upper <- vector(mode="list", length = p+1)
  for(i in 0:p)
  {
    cap <- vector(length=p+1)
    cap[1] <- 0
    for(j in 0:(p-i))
    {
      if (j==0 & i!=0)
      {
        uppertest <- full.var[order[1:i]]
        for(r in 1:length(var.list))
        {
          if(all(var.list[[r]] %in% uppertest)) cap[j+1]<-cap[j+1]+1
        }
      }else{
        if(j!=0){
          lowtest <- full.var[order[1:j]]
          uppertest <- full.var[order[1:(j+i)]]
          for(r in 1:length(var.list))
          {
            if(all(all(lowtest %in% var.list[[r]]),all(var.list[[r]] %in% uppertest))) cap[j+1]<-cap[j+1]+1
          }
        }
        if (j == 0){
          for (r in 1:length(var.list)){
            if (identical(var.list[[r]], character(0))) cap[j+1] <- cap[j+1] + 1
          }
        }
      }
    }
    freq[i+1] <- max(cap)/B
    maxlocation <- which.max(cap)
    if(maxlocation==1)
    {
      if (i != 0){
        lower[[i+1]] <- ''
        upper[[i+1]] <- full.var[order[1:i]]
      } else if (i == 0){
        lower[[1]] <- ''
        upper[[1]] <- ''
      }
    }else{
      lower[[i+1]] <- full.var[order[1:(maxlocation-1)]]
      upper[[i+1]] <- full.var[order[1:(maxlocation-1+i)]]
    }
  }
  result <- list(freq=freq,lower=lower,upper=upper)
  return(result)
}

# Get boostrap models from the modified lasso
BOOT.MODI.LASSO <- function(x, dep.var.index, r, lmbd, seed, full.var){
  var_lasso <- vector(mode="list", length=r)
  set.seed(seed=seed)
  for (i in 1:r){
    ind <- sample(1:nrow(x),nrow(x),replace=T)
    boot.data <- x[ind,]
    opt.lambda <- cv.glmnet(x=as.matrix(boot.data[,-dep.var.index]),y=boot.data[,dep.var.index],alpha=1)$lambda.min
    lasso.fit <- glmnet(x=as.matrix(boot.data[,-dep.var.index]),y=boot.data[,dep.var.index],family='gaussian',alpha=1)
    beta <- coef(lasso.fit,s=opt.lambda)[,1][which(abs(coef(lasso.fit,s=opt.lambda)[,1]) > lmbd)]
    var_lasso[[i]] <- full.var[full.var%in%names(beta)[beta!=0]]
  }
  return(var_lasso)
}

# GET the residual bootsrap variable selection results by using adaptive lasso and
# you can give a specified lambda
RES.BOOT.CI2<-function(x, r, lmbd, seed, p, full.var){
  var.instances <- vector(mode="list",length=r)
  ## = adaLASSO = ##
  tau <- 1
  lasso_init <- glmnet(as.matrix(x[,1:p]),x$y) #nlambda = 100
  first.step.coef <- lasso_init$beta[,which.min(abs(lasso_init$lambda-lmbd))]
  penalty.factor <- abs(first.step.coef+1/sqrt(nrow(x)))^(-tau)
  adalasso <- glmnet(as.matrix(x[,1:p]),x$y,penalty.factor=penalty.factor)
  beta_est <- adalasso$beta[,which.min(abs(adalasso$lambda-lmbd))] #lmbd找最近的
  res_original <- x$y - as.matrix(x[,1:p]) %*% beta_est
  res_after_center <- res_original - mean(res_original)
  constant <- as.matrix(x[,1:p]) %*% beta_est
  set.seed(seed=seed)
  for(j in 1:r) {
    ind=sample(1:nrow(x),nrow(x),replace=T)
    new_response <- constant + res_after_center[ind]
    boot.data <- cbind(x[,1:p], new_response)
    colnames(boot.data)[p+1] <- "y"
    lasso_init_boot=glmnet(as.matrix(boot.data[,1:p]),boot.data$y)
    first.step.coef_boot=lasso_init_boot$beta[,which.min(abs(lasso_init_boot$lambda-lmbd))]
    penalty.factor_boot=abs(first.step.coef_boot+1/sqrt(nrow(x)))^(-tau)
    adalasso_boot=glmnet(as.matrix(boot.data[,1:p]),boot.data$y,penalty.factor=penalty.factor_boot)
    beta_est_boot=adalasso_boot$beta[,which.min(abs(adalasso_boot$lambda-lmbd))]
    var.instances[[j]]<-full.var[beta_est_boot!=0]
  }
  return(var.instances)
}

# Use LAD and Sqrt method to get residual bootstrap variable models
RES.BOOT.CI4 <- function(x, r, q, lmbd,seed, p, full.var){
  var<-vector(mode="list",length=r)
  fit<-flare::slim(X=as.matrix(x[,1:p]),Y=x$y,method="lq",q=q)
  if(lmbd==''){
    beta <- fit$beta[,which.min(abs(fit$lambda-1))]
  }else{
    beta <- fit$beta[,which.min(abs(fit$lambda-lmbd))]
  }

  res_original <- x$y - as.matrix(x[,1:p]) %*% beta
  res_after_center <- res_original - mean(res_original)
  constant <- as.matrix(x[,1:p]) %*% beta
  set.seed(seed=seed)
  for(j in 1:r) {
    ind=sample(1:nrow(x),nrow(x),replace=T)
    new_response <- constant + res_after_center[ind]
    boot.data <- cbind(x[,1:p], new_response)
    colnames(boot.data)[p+1] <- "y"
    tem <- flare::slim(X=as.matrix(boot.data[,1:p]),Y=boot.data$y,method = 'lq', q = q,verbose = FALSE)
    if(lmbd==''){
      var[[j]]<-full.var[tem$beta[,which.min(abs(fit$lambda-1))]!=0]
    }else{
      var[[j]]<-full.var[tem$beta[,which.min(abs(fit$lambda-lmbd))]!=0]
    }

  }
  return(var)
}



# get residual bootstrap variable selection models from SCAD and MCP
RES.BOOT.CI3 <- function(x, r, pnlt, lmbd, seed, p, full.var){
  var<-vector(mode="list",length=r)
  fit<-cv.ncvreg(X=x[,1:p],y=x$y,penalty=pnlt)
  fit2<-fit$fit
  if(lmbd == ''){
    beta<-fit2$beta[,fit$min]
  }else{
    beta<-fit2$beta[,which.min(abs(fit$lambda-lmbd))]
  }
  beta<-beta[2:(p+1)]
  res_original <- x$y - as.matrix(x[,1:p]) %*% beta
  res_after_center <- res_original - mean(res_original)
  constant <- as.matrix(x[,1:p]) %*% beta
  set.seed(seed=seed)
  for(j in 1:r) {
    ind=sample(1:nrow(x),nrow(x),replace=T)
    new_response <- constant + res_after_center[ind]
    boot.data <- cbind(x[,1:p], new_response)
    colnames(boot.data)[p+1] <- "y"
    tem <- cv.ncvreg(X=as.matrix(boot.data[,1:p]),y=boot.data$y,penalty=pnlt)
    fit2_tem<-tem$fit
    if(lmbd==''){
      beta_tem<-fit2_tem$beta[,tem$min]
    }else{
      beta_tem<-fit2_tem$beta[,which.min(abs(tem$lambda-lmbd))]
    }
    beta_tem <- beta_tem[2:(p+1)]
    var[[j]]<-names(beta_tem)[which(beta_tem!=0)]
  }
  return(var)
}

# GET residual boostrap variable selection models with stepwise BIC
RES.BOOT.CI5 <- function(x, p, r, lmbd, seed, full.var){
  if(lmbd==''){
    var<-vector(mode="list",length=r)
    fit<-regsubsets(y~.,data=x,method="seqrep",nvmax=p)
    if (names(coef(fit,which.min(summary(fit)$bic)))[1] == '(Intercept)'){
      beta <- vector(mode = 'numeric', length = p)
      beta[full.var%in%names(coef(fit,which.min(summary(fit)$bic)))] <- coef(fit,which.min(summary(fit)$bic))[-1]
    } else {
      beta <- vector(mode = 'numeric', length = p)
      beta[full.var%in%names(coef(fit,which.min(summary(fit)$bic)))] <- coef(fit,which.min(summary(fit)$bic))
    }
    res_original <- x$y - as.matrix(x[,1:p]) %*% beta
    res_after_center <- res_original - mean(res_original)
    constant <- as.matrix(x[,1:p]) %*% beta
    set.seed(seed=seed)
    for(j in 1:r) {
      ind=sample(1:nrow(x),nrow(x),replace=T)
      new_response <- constant + res_after_center[ind]
      boot.data <- cbind(x[,1:p], new_response)
      colnames(boot.data)[p+1] <- "y"
      tem <- regsubsets(y~.,data=boot.data,method="seqrep",nvmax=p)
      if (names(coef(tem,which.min(summary(tem)$bic)))[1] == '(Intercept)'){
        var[[j]]<-names(coef(tem,which.min(summary(tem)$bic)))[-1]
      } else {
        var[[j]]<-names(coef(tem,which.min(summary(tem)$bic)))
      }
    }
    return(var)
  }else{
    var<-vector(mode="list",length=r)
    fit<-regsubsets(y~.,data=x,method="seqrep",nvmax=p)
    if (names(coef(fit,which.min(abs(summary(fit)$bic - lmbd))))[1] == '(Intercept)'){
      beta <- vector(mode = 'numeric', length = p)
      beta[full.var%in%names(coef(fit,which.min(abs(summary(fit)$bic - lmbd))))] <- coef(fit,which.min(summary(fit)$bic))[-1]
    } else {
      beta <- vector(mode = 'numeric', length = p)
      beta[full.var%in%names(coef(fit,which.min(abs(summary(fit)$bic - lmbd))))] <- coef(fit,which.min(summary(fit)$bic))
    }
    res_original <- x$y - as.matrix(x[,1:p]) %*% beta
    res_after_center <- res_original - mean(res_original)
    constant <- as.matrix(x[,1:p]) %*% beta
    for(j in 1:r) {
      ind=sample(1:nrow(x),nrow(x),replace=T)
      new_response <- constant + res_after_center[ind]
      boot.data <- cbind(x[,1:p], new_response)
      colnames(boot.data)[p+1] <- "y"
      tem <- regsubsets(y~.,data=boot.data,method="seqrep",nvmax=p)
      if (names(coef(tem,abs(which.min(summary(tem)$bic - lmbd))))[1] == '(Intercept)'){
        var[[j]]<-names(coef(tem,abs(which.min(summary(tem)$bic - lmbd))))[-1]
      } else {
        var[[j]]<-names(coef(tem,abs(which.min(summary(tem)$bic - lmbd))))
      }
    }
    return(var)
  }
}

BOOT.MODI.ELASTIC <- function(x, dep.var.index, r, lmbd, seed, full.var){
  var_lasso<-vector(mode="list",length=r)
  set.seed(seed=seed)
  for (i in 1:r){
    ind=sample(1:nrow(x),nrow(x),replace=T)
    boot.data<-x[ind,]
    opt.lambda<-cv.glmnet(x=as.matrix(boot.data[,-dep.var.index]),y=boot.data[,dep.var.index],alpha=1)$lambda.min
    elastic.fit<-glmnet(x=as.matrix(boot.data[,-dep.var.index]),y=boot.data[,dep.var.index],family='gaussian',alpha=0.5)
    beta<-coef(elastic.fit,s=opt.lambda)[,1][which(abs(coef(elastic.fit,s=opt.lambda)[,1]) > lmbd)]
    var_lasso[[i]]<-full.var[full.var%in%names(beta)[beta!=0]]
  }
  return(var_lasso)
}

getmcb <- function(result, c){
  n = length(result$freq)
  fit_freq = c(result$freq)[result$freq - c >= 0]
  best_fit = which.min(fit_freq - c) + n - sum(result$freq - c >= 0)
  mc = list()
  mc$lbm <- result$lower[[best_fit]]
  mc$ubm <- result$upper[[best_fit]]
  mc$bcr <- result$freq[best_fit]
  return(mc)
}

getmcbt <- function(result){
  mcbf <- data.frame(lbm = matrix(result$lower), bcr = result$freq, ubm = matrix(result$upper))
  n = length(mcbf$lbm)
  mcbf$width <- c(0:(n-1))
  mcbf <- mcbf[,c('width','lbm','bcr','ubm')]
  return(mcbf)
}

check_flare <- function() {
  if (!requireNamespace("flare", quietly = TRUE)) {
    stop("Package flare needed for SQRT or MCP method to work. Please install it.",
         call. = FALSE)
  }
}

# main function
mcb.compare <- function(x, y, B=200, lambdas=NA, methods = NA, level=0.95,seed=122){

  # methods = c('adalasso', 'lasso', 'SCAD', 'MCP','stepwise','LAD','SQRT', NA)

  x = as.data.frame(x)
  n = length(x)
  for (i in 1:n) {
    name = paste('x',i,sep='')
    colnames(x)[i] = name
  }
  y = as.data.frame(y)

  data <- cbind(x,y)

  full.var <- paste(colnames(x),sep="")
  colnames(data)[dim(data)[2]] = "y"
  p <- dim(data)[2]-1
  r <- B


  final_result = ''

  if(is.na(methods)){
    methods <- c('aLasso', 'Lasso','Elastic', 'SCAD', 'MCP','stepwise','LAD','SQRT')
  }
  if(is.na(lambdas)){
    lambdas = ''
  }

  mcbfit = list()
  mcbframe = list()

  # adaptivelasso
  if('aLasso' %in% methods){
    if(lambdas == ''){
      var_adaptivelasso <- RES.BOOT.CI2(data, r, lmbd=1, seed=seed, p, full.var)
    }else{
      var_adaptivelasso <- RES.BOOT.CI2(data, r, lmbd=lambdas[which(methods=='adaptiveLasso')], seed=seed, p, full.var)
    }
    var_01_ada_lasso<-f01(var_adaptivelasso, full.var, p)
    result_aLasso <- CI(var_adaptivelasso, var_01_ada_lasso, p, r)
    mcbfit$aLasso <- getmcb(result_aLasso, c=level)
    mcbframe$aLasso <- getmcbt(result_aLasso)
  }

  # lasso
  if('Lasso' %in% methods){
    if(lambdas == ''){
      var_lasso <- BOOT.MODI.LASSO(data, p + 1, r, lmbd=0.05, seed=seed, full.var)
    }else{
      var_lasso <- BOOT.MODI.LASSO(data, p + 1, r, lmbd=lambdas[which(methods=='Lasso')],seed=seed)
    }
    var_01_lasso<-f01(var_lasso, full.var, p)
    result_Lasso <- CI(var_lasso, var_01_lasso, p, r)
    mcbfit$Lasso<-getmcb(result_Lasso, c=level)
    mcbframe$Lasso<-getmcbt(result_Lasso)
  }

  # Elastic
  if('Elastic' %in% methods){
    if(lambdas == ''){
      var_elastic <- BOOT.MODI.ELASTIC(data, p + 1, r, lmbd=0.05, seed=seed, full.var)
    }else{
      var_elastic <- BOOT.MODI.ELASTIC(data, p + 1, r, lmbd=lambdas[which(methods=='Lasso')], seed=seed, full.var)
    }
    var_01_elastic<-f01(var_elastic, full.var, p)
    result_Elastic <- CI(var_elastic, var_01_elastic, p, r)
    mcbfit$Elastic<-getmcb(result_Elastic, c=level)
    mcbframe$Elastic<-getmcbt(result_Elastic)
  }

  # SCAD
  if('SCAD' %in% methods){
    if(lambdas==''){
      var_SCAD <- RES.BOOT.CI3(data, r, pnlt='SCAD',lmbd=lambdas, seed=seed, p, full.var)
    }else{
      var_SCAD <- RES.BOOT.CI3(data, r, pnlt='SCAD',lmbd=lambdas[which(methods=='SCAD')], seed=seed, p, full.var)
    }
    var_01_SCAD<-f01(var_SCAD, full.var, p)
    result_SCAD <- CI(var_SCAD, var_01_SCAD, p, r)
    mcbfit$SCAD <- getmcb(result_SCAD, c=level)
    mcbframe$SCAD <- getmcbt(result_SCAD)
  }

  # MCP
  if('MCP' %in% methods){
    if (lambdas==''){
      var_MCP <- RES.BOOT.CI3(data, r, pnlt='MCP',lmbd=lambdas, seed=seed, p, full.var)
    }else{
      var_MCP <- RES.BOOT.CI3(data, r, pnlt='MCP',lmbd=lambdas[which(methods=='MCP')], seed=seed, p, full.var)
    }
    var_01_MCP<-f01(var_MCP, full.var, p)
    result_MCP <- CI(var_MCP, var_01_MCP, p, r)
    mcbfit$MCP <- getmcb(result_MCP, c=level)
    mcbframe$MCP <- getmcbt(result_MCP)
  }

  # stepwise
  if('stepwise' %in% methods){
    if(lambdas==''){
      var_stepwise <- RES.BOOT.CI5(data, p, r,lmbd=lambdas, seed=seed, full.var)
    }else{
      var_stepwise <- RES.BOOT.CI5(data, p, r,lmbd=lambdas[which(methods=='stepwise')], seed=seed, full.var)
    }
    var_01_stepwise <- f01(var_stepwise, full.var, p)
    result_stepwise <- CI(var_stepwise, var_01_stepwise, p, r)
    mcbfit$stepwise <- getmcb(result_stepwise, c=level)
    mcbframe$stepwise <- getmcbt(result_stepwise)
  }

  # LAD
  if('LAD' %in% methods){
    check_flare()
    if(lambdas==''){
      var_LAD <- RES.BOOT.CI4(data, r, q = 1, lmbd = lambdas,seed=seed, p, full.var)
    }else{
      var_LAD <- RES.BOOT.CI4(data, r, q = 1, lmbd = lambdas[which(methods=='LAD')], seed=seed, p, full.var)
    }
    var_01_LAD <- f01(var_LAD, full.var, p)
    result_LAD <- CI(var_LAD, var_01_LAD, p, r)
    mcbfit$LAD <- getmcb(result_LAD, c=level)
    mcbframe$LAD <- getmcbt(result_LAD)
  }

  # SQRT
  if('SQRT' %in% methods){
    check_flare()
    if(lambdas==''){
      var_SQRT <- RES.BOOT.CI4(data, r, q = 2, lmbd = lambdas, seed=seed, p, full.var)
    }else{
      var_SQRT <- RES.BOOT.CI4(data, r, q = 2, lmbd = lambdas[which(methods=='SQRT')], seed=seed, p, full.var)
    }
    var_01_SQRT <- f01(var_SQRT, full.var, p)
    result_SQRT <- CI(var_SQRT, var_01_SQRT, p, r)
    mcbfit$SQRT <- getmcb(result_SQRT,c=level)
    mcbframe$SQRT <- getmcbt(result_SQRT)
  }

  final_result = ''
  select_name = paste('result_',methods,sep = '')
  for(i in select_name){
    final_result <- cbind(final_result,get(i)$freq)
  }

  final_result <- final_result[,c(2:(length(methods)+1))]
  colnames(final_result) <- methods

  all_result = list()

  # ggplot2 - mucplot
  df <- data.frame(final_result)
  df$x <- seq(0,1,length.out = n+1)
  df$x <- round(df$x,digits = 3)
  df <- melt(df, id=c('x'))
  df$value <- as.numeric(df$value)
  muc <- ggplot(df, aes(x=df$x, y=df$value, color=df$variable, shape=df$variable)) + geom_line() + labs(x = "w/p") + labs(y = "freq") + labs(colour = "Method") + ylim(0,1) + scale_x_continuous(breaks = df$x)
  all_result$mucplot <- muc
  all_result$mcb <- mcbfit
  all_result$mcbframe <- mcbframe

  return(all_result)
}

mcb <- function(x, y, B=200, lambda=NA, method = 'Lasso', level=0.95, seed = 122){
  # methods = c('adalasso', 'lasso', 'SCAD', 'MCP','stepwise','LAD','SQRT', NA)

  x = as.data.frame(x)
  n = length(x)
  for (i in 1:n) {
    name = paste('x',i,sep='')
    colnames(x)[i] = name
  }
  y = as.data.frame(y)

  data <- cbind(x,y)
  full.var <- paste(colnames(x),sep="")
  colnames(data)[dim(data)[2]] = "y"
  p <- dim(data)[2]-1
  r = B

  if(is.na(lambda)){
    lambda = ''
  }

  # adaptivelasso
  if(method == 'aLasso'){
    if(lambda == ''){
      var_adaptivelasso <- RES.BOOT.CI2(data, r, lmbd=1, seed=seed, p, full.var)
    }else{
      var_adaptivelasso <- RES.BOOT.CI2(data, r, lmbd=lambda, seed=seed, p, full.var)
    }
    var_01_ada_lasso<-f01(var_adaptivelasso, full.var, p)
    result <- CI(var_adaptivelasso, var_01_ada_lasso, p, r)
  }

  # lasso
  if(method == 'Lasso'){
    if(lambda == ''){
      var_lasso <- BOOT.MODI.LASSO(data, p + 1, r, lmbd=0.05, seed=seed, full.var)
    }else{
      var_lasso <- BOOT.MODI.LASSO(data, p + 1, r, lmbd=lambda,seed=seed)
    }
    var_01_lasso<-f01(var_lasso, full.var, p)
    result <- CI(var_lasso, var_01_lasso, p, r)
  }

  # Elastic
  if(method == 'Elastic'){
    if(lambda == ''){
      var_elastic <- BOOT.MODI.ELASTIC(data, p + 1, r, lmbd=0.05, seed=seed, full.var)
    }else{
      var_elastic <- BOOT.MODI.ELASTIC(data, p + 1, r, lmbd=lambda, seed=seed, full.var)
    }
    var_01_elastic <- f01(var_elastic, full.var, p)
    result <- CI(var_elastic, var_01_elastic, p, r)
  }

  # SCAD
  if(method == 'SCAD'){
    var_SCAD <- RES.BOOT.CI3(data, r, pnlt='SCAD',lmbd=lambda, seed=seed, p, full.var)
    var_01_SCAD <- f01(var_SCAD, full.var, p)
    result <- CI(var_SCAD, var_01_SCAD, p, r)
  }

  # MCP
  if(method=='MCP'){
    var_MCP <- RES.BOOT.CI3(data, r, pnlt='MCP',lmbd=lambda, seed=seed, p, full.var)
    var_01_MCP <- f01(var_MCP, full.var, p)
    result <- CI(var_MCP, var_01_MCP, p, r)
  }

  # stepwise
  if(method=='stepwise'){
    var_stepwise <- RES.BOOT.CI5(data, p, r,lmbd=lambda, seed=seed, full.var)
    var_01_stepwise <- f01(var_stepwise, full.var, p)
    result <- CI(var_stepwise, var_01_stepwise, p, r)
  }

  # LAD
  if(method=='LAD'){
    check_flare()
    var_LAD <- RES.BOOT.CI4(data, r, q = 1, lmbd = lambda,seed=seed, p, full.var)
    var_01_LAD <- f01(var_LAD, full.var, p)
    result <- CI(var_LAD, var_01_LAD, p, r)
  }

  # SQRT
  if(method=='SQRT'){
    check_flare()
    var_SQRT <- RES.BOOT.CI4(data, r, q = 2, lmbd = lambda, seed=seed, p, full.var)
    var_01_SQRT <- f01(var_SQRT, full.var, p)
    result <- CI(var_SQRT, var_01_SQRT, p, r)
  }

  all_result <- list()

  # ggplot2 - mucplot
  df <- data.frame(result$freq)
  df$x <- seq(0,1,length.out = n+1)
  df$x <- round(df$x,digits = 3)
  muc <- ggplot(df, aes(x=df$x, y=result$freq)) + geom_line() + labs(x = "w/p") + labs(y = "freq") + labs(colour = "Method") + ylim(0,1) + scale_x_continuous(breaks = df$x)
  all_result$mucplot <- muc

  fit_freq = c(result$freq)[result$freq - level >= 0]
  best_fit = which.min(fit_freq - level) + n+1 - sum(result$freq - level >= 0)
  mcb = list()
  mcb$lbm <- result$lower[[best_fit]]
  mcb$ubm <- result$upper[[best_fit]]
  mcb$bcr <- result$freq[best_fit]
  all_result$mcb <- mcb

  mcbframe <- data.frame(lbm = matrix(result$lower), bcr = result$freq, ubm = matrix(result$upper))
  mcbframe$width <- c(0:n)
  mcbframe <- mcbframe[,c('width','lbm','bcr','ubm')]
  all_result$mcbframe <- mcbframe

  return(all_result)
}


