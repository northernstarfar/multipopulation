install.packages("Rsolnp")
install.packages("nloptr")
install.packages("mvtnorm")
install.packages("MBESS")
library(Rsolnp)
library(nloptr)
library(mvtnorm)
library(MBESS)

#Set the seed number.
set.seed(21)
#This nofsimu is the number of simulation we would like to run.
nofsimu=100
#This s1.3 is the list to store all simulation results.
s1.2 <- vector("list",nofsimu)
#This s1.3y_m is the list to store all y from the main model in all nofsimu simulations.
s1.2y_m <- vector("list",nofsimu)
#This s1.3y_1st is the list to store all y from the 1st sets of the helper models in all nofsimu simulations.
s1.2y_1st <- vector("list",nofsimu)
#This s1.3y_2nd is the list to store all y from this  2nd sets of the helper models in all nofsimu simulations.
s1.2y_2nd <- vector("list",nofsimu)
#This s1.3y_3rd is the list to store all y from the 3rd sets of the helper models in all nofsimu simulations.
s1.2y_3rd <- vector("list",nofsimu)
#This s1.3y_all is the list to store all y from the 3rd sets of the helper models in all nofsimu simulations.
s1.2y_all <- vector("list",nofsimu)
#This s1.3estpar is the list to store all parameter estimates in all nofsimu simulations.
s1.2estpar <- vector("list",nofsimu)
s1.2y_ob_out <- vector("list",nofsimu)
s1.2testdata <- vector("list",nofsimu)
s1.2estpar <- vector("list",nofsimu)

for (indexi in 1:nofsimu) {
  #Let us know which simulation we are running.
  print(indexi)
  
  #Data generation part
  
  #This function yijdata_m is used to generate Y_ij data for the logistic model.
  #This can be used to generate Y_ij data for the main model only.
  #The input x should be p dimensional, z should be q dimensional.
  #tbeta should be p dimensional vector.
  #talpha should be q dimensional vector.
  #tgamma should be r dimensional vector.
  #Note that even tough we said tgamma should be r dimensional. We actually can only put 1 dimensional tgamma here,as gamma is an intercept.
  
  yijdata_m <- function(x,z,tbeta,talpha,tgamma,tdelta,zlast,eps) {
    
    pyij=sum(x*tbeta)+sum(z*talpha)+sum(tgamma)+sum(zlast*tdelta)+eps
    
    return (pyij)
  }
  
  #This function yijdata_m_exp is used to generate the expectation of Y_ij data for the logistic model.
  #This can be used to generate the expectation of Y_ij data for the main model only.
  #The input x should be p dimensional, z should be q dimensional.
  #tbeta should be p dimensional vector.
  #talpha should be q dimensional vector.
  #tgamma should be r dimensional vector.
  #Note that even tough we said tgamma should be r dimensional. We actually can only put 1 dimensional tgamma here,as gamma is an intercept.
  
  yijdata_m_exp <- function(x,z,tbeta,talpha,tgamma,tdelta,zlast) {
    lengthofz=length(z)
    pyij=sum(x*tbeta)+sum(z*talpha)+sum(tgamma)+sum(zlast*tdelta)
    
    return (pyij)
  }
  
  #This function yijdata_1st is used to generate Y_ij data for the logistic model.
  #This can be used to generate Y_ij data for the 1st sets of helper models only.
  #The input x should be p dimensional, z should be q dimensional, u should be r dimensional.
  #tbeta should be p dimensional vector.
  #talpha should be q dimensional vector..
  #tgamma should be r dimensional vector..
  #Note that unlike tgamma in yijdata, the dimension of tgamma here could be greater than 1, as here gamma is no longer an intercept in this model.
  #Note that in the paper, the 1st sets of the helper models don't have u in the linear part. This is because in the paper we combined u
  #into the Z part and thus the gamma associated with u is now integrated into alpha associated with z.
  
  yijdata_1st <- function(x,z,u,tbeta,talpha,tgamma,tdelta,eps) {
    lengthofz=length(z)
    pyij=sum(x*tbeta)+sum(z*talpha)+sum(u*tgamma)+eps
    
    return (pyij)
  }
  
  #This function yijdata_1st_mis is used to generate Y_ij data for the misspecification part of the logistic model.
  #This can be used to generate Y_ij data for the misspecification part of 1st sets of helper models only.
  #The input x should be p dimensional, z should be q dimensional, u should be r dimensional.
  #tbeta should be p dimensional vector.
  #talpha should be q dimensional vector.
  #tgamma should be r dimensional vector.
  #tdelta should be p dimensional as it is the parameter for the quadratic term of X.
  #Note that unlike tgamma in yijdata, the dimension of tgamma here could be greater than 1, as here gamma is no longer an intercept in this model.
  #Note that in the paper, the 1st sets of the helper models don't have u in the linear part. This is because in the paper we combined u
  #into the Z part and thus the gamma associated with u is now integrated into alpha associated with z.
  
  yijdata_1st_mis <- function(x,z,u,tbeta,talpha,tgamma,tdelta,zlast,eps) {
    lengthofz=length(z)
    
    pyij=sum(x*tbeta)+sum(z*talpha)+sum(u*tgamma)+sum(zlast*tdelta)+eps
    
    return (pyij)
  }
  
  
  #This function ydata_m is used to generate Y data for all j,i for the logistic model.
  #This can be used to genrate Y data for the main model only.
  #When using ydata to generate Y data for the main model, make sure the input x be a list containing only 1 m*p matrix,
  #and the input z be a list containing only 1 m*q matrix.
  #nhelper1 is the number of populations from the first sets of the helper models.
  #tbeta is a p dimensional vector.
  #talpha is a q dimensional vector.
  #tgamma is a r dimensional vector.
  
  ydata_m <- function(x,z,tbeta,talpha,tgamma,tdelta,zlast,eps) {
    i=length(x)
    j=nrow(x[[1]])
    yall=matrix(0,i,j)
    
    for (k in 1:i) {
      for (l in 1:j) {
        yall[k,l]=yijdata_m(x[[k]][l,],z[[k]][l,],tbeta,talpha,tgamma,tdelta,zlast[[k]][l],eps[l])
      }
    }
    return (yall)
  }
  
  #This function ydata_m_exp is used to generate the expectation of Y data for all j,i for the logistic model.
  #This can be used to genrate the expectation of Y data for the main model only.
  #When using ydata to generate the expectation of Y data for the main model, make sure the input x be a list containing only 1 m*p matrix,
  #and the input z be a list containing only 1 m*q matrix.
  #nhelper1 is the number of populations from the first sets of the helper models.
  #tbeta is a p dimensional vector.
  #talpha is a q dimensional vector.
  #tgamma is a r dimensional vector.
  
  ydata_m_exp <- function(x,z,tbeta,talpha,tgamma,tdelta,zlast) {
    i=length(x)
    j=nrow(x[[1]])
    yall=matrix(0,i,j)
    
    for (k in 1:i) {
      for (l in 1:j) {
        yall[k,l]=yijdata_m_exp(x[[k]][l,],z[[k]][l,],tbeta,talpha,tgamma,tdelta,zlast[[k]][l])
      }
    }
    return (yall)
  }
  
  #This function ydata_1st is used to generate Y data for all j,i for the logistic model.
  #This can be used to genrate Y data for the 1st sets of helper model only.
  #When using ydata1sthelper to gnerate Y data for the 1st sets of helper models, make sure the input x be a list containing nhelper1 m*p matrix,
  #the input z be a list containing nhelper1 m*q matrix, and the input u be a list containing nhelper1 m*r matrix.
  #nhelper1 is the number of populations from the first sets of the helper models.
  #tbeta is a p dimensional vector.
  #talpha is a q dimensional vector.
  #tgamma is a r dimensional vector.
  #Note that in the paper, the 1st sets of the helper models don't have u in the linear part. This is because in the paper we combined u
  #into the Z part and thus the gamma associated with u is now integrated into alpha associated with z.
  
  ydata_1st <- function(x,z,u,tbeta,talpha,tgamma,tdelta,eps) {
    i=length(x)
    j=nrow(x[[1]])
    yall=matrix(0,i,j)
    
    for (k in 1:i) {
      for (l in 1:j) {
        yall[k,l]=yijdata_1st(x[[k]][l,],z[[k]][l,],u[[k]][l,],tbeta,talpha,tgamma,tdelta,eps[l])
      }
    }
    return (yall)
  }
  
  #This function ydata_1st_mis is used to generate Y data for all j,i for the misspecification part of the logistic model.
  #This can be used to genrate Y data for the misspecification part of the 1st sets of helper model only.
  #When using ydata1sthelper to gnerate Y data for the misspecification part of the 1st sets of helper models, make sure the input x be a list containing nhelper1 m*p matrix,
  #the input z be a list containing nhelper1 m*q matrix, and the input u be a list containing nhelper1 m*r matrix.
  #nhelper1 is the number of populations from the first sets of the helper models.
  #tbeta is a p dimensional vector.
  #talpha is a q dimensional vector.
  #tgamma is a r dimensional vector.
  #tdelta should be p dimensional as it is the parameter for the quadratic term of X.
  #Note that in the paper, the 1st sets of the helper models don't have u in the linear part. This is because in the paper we combined u
  #into the Z part and thus the gamma associated with u is now integrated into alpha associated with z.
  
  ydata_1st_mis <- function(x,z,u,tbeta,talpha,tgamma,tdelta,zlast,eps) {
    i=length(x)
    j=nrow(x[[1]])
    yall=matrix(0,i,j)
    
    for (k in 1:i) {
      for (l in 1:j) {
        yall[k,l]=yijdata_1st_mis(x[[k]][l,],z[[k]][l,],u[[k]][l,],tbeta,talpha,tgamma,tdelta,zlast[[k]][l],eps[l])
      }
    }
    return (yall)
  }
  
  
  #This function make.data.setforj_m_2nd is used to obtain the suitable data matrix for the glm function for
  #population from main model and populations from 2nd sets of helper models.
  #m is the observation number for each simulation, corresponding to the index i in the paper.
  #p is the length of the true beta, which is also the length of each X_ji in the paper.
  #q is the length of the ture alpha, which is also the length of each Z_ji in the paper.
  #r is the length of the true gamma.
  #y is an m length vector.
  #x is a m*p matrix from a list of matrix containing all m*p matrix for x, this m*p matrix stands for all observations of x from the same poppulation.
  #z is a m*q matrix from a list of matrix containing all m*q matrix for z, this m*q matrix stands for all observations of z from the same poppulation.
  
  make.data.setforj_m_2nd <- function(m,p,q,r,y,x,z) {
    ## vector indicating which family each observation belongs
    j=rep(m,1)
    family <- rep(1:1,times=j)
    
    ## make y into a vector of values. 
    new.y <- as.vector(t(y))
    
    ## make x into a matrix of values.
    x.values <- x
    
    ## make z into a matrix of values.
    z.values <- z
    
    # make a matrix g.values of 1 to serve as the covariate 1 of gamma_j
    g.values=matrix(1,m,r)
    
    data.set.out <- data.frame(cbind(family,new.y,x.values,z.values,g.values))
    colnames(data.set.out) <- c("family","Y",paste("X",1:p,sep=""),paste("Z",1:q,sep=""),paste("G",1:r,sep=""))
    
    return (data.set.out)
  }
  
  
  
  #This function make.data.setforj_1st is used to obtain the suitable data matrix for the glm function for the populations from 1st sets of helper models.
  #m is the observation number for each simulation, corresponding to the index i in the paper.
  #p is the length of the true beta, which is also the length of each X_ji in the paper.
  #q is the length of the ture alpha, which is also the length of each Z_ji in the paper.
  #r is the length of the true gamma, which is also the length of each U_ji in the paper.
  #y is an m length vector.
  #x is a m*p matrix from a list of matrix containing all m*p matrix for x, this m*p matrix stands for all observations of x from the same poppulation.
  #z is a m*q matrix from a list of matrix containing all m*q matrix for z, this m*q matrix stands for all observations of z from the same poppulation.
  #u is a m*r matrix from a list of matrix containing all m*r matrix for u, this m*r matrix stands for all observations of u from the same poppulation.
  #if we set the u to be a m*r matrix of 1, then it could be used as make.data.setforj function.
  
  make.data.setforj_1st <- function(m,p,q,r,y,x,z,u) {
    ## vector indicating which family each observation belongs
    j=rep(m,1)
    family <- rep(1:1,times=j)
    
    ## make y into a vector of values. 
    new.y <- as.vector(t(y))
    
    ## make x into a matrix of values.
    x.values <- x
    
    ## make z into a matrix of values.
    z.values <- z
    
    ## make u into a matrix of values.
    u.values <- u
    
    data.set.out <- data.frame(cbind(family,new.y,x.values,z.values,u.values))
    colnames(data.set.out) <- c("family","Y",paste("X",1:p,sep=""),paste("Z",1:q,sep=""),paste("G",1:r,sep=""))
    
    return (data.set.out)
  }
  
  
  
  ##Estimating Part
  
  #This function mainmodel is to obtain the MLE of the parameters of the Logistic main Model.
  #The input data.set is from result generated by make.data.setforj_m_2nd function.
  #p is the length of the true beta, which is also the length of each X_ij in the paper.
  #q is the length of the ture alpha, which is also the length of each Z_ij in the paper.
  #r is the length of the true gamma.
  #Note that in the main model, the data is just from one population, which is population 1 in the paper,
  #So when making dataset for the mainmodel function, we should set n to be 1.
  #Since n here is just 1, we have to make sure m is large enough so we could have convergence of MLE.
  
  mainmodel <- function(data.set,p,q,r) {
    xnam <- paste("X",1:p,sep="")
    znam <- paste("Z",1:q,sep="")
    gnam <- paste("G",1:r,sep="")
    totalnam <- c(xnam,znam,gnam)
    
    fmla <- as.formula(paste("Y~-1+",paste(totalnam,collapse="+")))
    fm <- lm(fmla,data=data.set)
    
    ## Extract needed information for coefficients.
    Vcov <- vcov(fm, useScale = FALSE)
    v.delta <- Vcov
    betas <- coef(fm)
    se <- sqrt(diag(Vcov))
    
    zval <- betas / se
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
    out.fixed <- cbind(betas, se, zval, pval)
    
    list(out.fixed=out.fixed)
    #return (fm)
  }
  
  #This function firsthelpmodel is used to obtain the MLE of the parameters of the 1st sets of helper Logistic Models.
  #The input data.set is from result generated by make.data.setforj_1st function.
  #p is the length of the true beta, which is also the length of each X_ij in the paper.
  #q is the length of the ture alpha, which is also the length of each Z_ij in the paper.
  #r is the length of the true gamma.
  #Note that in the first sets of helper models, the data can be from several populations.
  #So when making dataset for the firsthelpmodel function, the n could be arbitrary. 
  
  firsthelpmodel <- function(data.set,p,q,r) {
    xnam <- paste("X",1:p,sep="")
    znam <- paste("Z",1:q,sep="")
    gnam <- paste("G",1:r,sep="")
    totalnam <- c(xnam,znam,gnam)
    
    fmla <- as.formula(paste("Y~-1+",paste(totalnam,collapse="+")))
    fm <- lm(fmla,data=data.set)
    
    ## Extract needed information for coefficients.
    Vcov <- vcov(fm, useScale = FALSE)
    v.delta <- Vcov
    betas <- coef(fm)
    se <- sqrt(diag(Vcov))
    
    zval <- betas / se
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
    out.fixed <- cbind(betas, se, zval, pval)
    
    list(out.fixed=out.fixed)
    #return (fm)
  }
  
  
  #Predition Part
  
  #This function yjpre is used to form the j-th prediction of Y using the estimator B_j_hat and alpha_1_hat, gamma_1_hat in step 2 of the paper.
  #Because our main model is Logistic Model. E(Y) is just Pr(Y=1).
  #The input x is a p dimensional vector, same length as tbeta.
  #The input z is a q dimensional vector, same length as talpha.
  #tbeta is the MLE of the j-th beta, which is from the MLE of the j-th model.
  #talpha is the MLE of the 1st alpha, which is from the MLE of the main model.
  #tgamma is the MLE of the 1st gamma, which is from the MLE of the main model.
  
  yjpre <- function(x,z,tbeta,talpha,tgamma) {
    
    yjpre=sum(x*tbeta)+sum(z*talpha)+sum(tgamma)
    
    return (yjpre)
  }
  
  #This function ywhat is used to combine Y_j_hat's to construct a function of w in step 3 of the paper.
  #The input w is an n dimensional vector of weights.
  #The input y is an n dimensional vector where each component is j-th predction of Y in step 2 where j is from 1 to n(N in the paper).
  
  ywhat <- function(w,y){
    n=length(y)
    ywhat=0
    
    for (i in 1:n) {
      ywhat=ywhat+w[i]*y[i]
    }
    return (ywhat)
  }
  
  
  #This function y1iexcludebfw is used to calculate Y_j_hat of Y_1i_hat^(-i)(w) in step 4 of the paper.
  #The input data set should come from make.data.setforj_m_2nd function, when we are making the dataset, we only need to inlcude observations from
  #the 1st population.
  #The input constbeta is an N-1 length list whose component are from vector beta_2_hat to vector beta_N_hat because when calculating Y_1i_hat^(-i)(w), 
  #we don't need to re-calculate beta_2_hat to beta_N_hat, only beta_1_hat varies as we exclude i-th observation from the 1st population,
  #when calculate MLE of beta_1_hat,alpha_1_hat,gamma_1_hat.
  #p is the length of the true beta, which is also the length of each X_ij in the paper.
  #q is the length of the ture alpha, which is also the length of each Z_ij in the paper.
  #r is the length of the true gamma.
  #i is indicating which observation we are excluding from the 1st population when calculating MLE. The range of i is from 1 to N_1,
  #N_1 is the number of observations in 1st population.
  
  
  y1iexcludebfw <- function(data.set,constbeta,p,q,r,i) {
    m=nrow(data.set)
    n=length(constbeta)+1
    tempbeta=vector("list",n)
    y=rep(0,n)
    
    tempdata.set=data.set[-i,]
    tempout=mainmodel(tempdata.set,p,q,r)
    
    beta1=as.vector(tempout[[1]][1:p,1])
    alpha1=as.vector(tempout[[1]][(p+1):(p+q),1])
    gamma1=as.vector(tempout[[1]][(p+q+1):(p+q+r),1])
    
    for (j in 2:n) {
      tempbeta[[j]]=constbeta[[j-1]]
    }
    tempbeta[[1]]=beta1
    
    
    for (j in 1:n) {
      y[j]=yjpre(as.vector(unlist(data.set[i,3:(2+p)])),as.vector(unlist(data.set[i,(2+p+1):(2+p+q)])),tempbeta[[j]],alpha1,gamma1)
      
    }
    
    return (y)
  }
  
  #This function y1iexcludebfwall is used to calculate Y_j_hat of Y_1i_hat^(-i)(w) for all i in step 4 of the paper.
  #The input data set should come from make.data.setforj_m_2nd function, when we are making the dataset, we only need to inlcude observations from
  #the 1st population.
  #The input constbeta is an N-1 length list whose component are from vector beta_2_hat to vector beta_N_hat because when calculating Y_1i_hat^(-i)(w), 
  #we don't need to re-calculate beta_2_hat to beta_N_hat, only beta_1_hat varies as we exclude i-th observation from the 1st population,
  #when calculate MLE of beta_1_hat,alpha_1_hat,gamma_1_hat.
  #p is the length of the true beta, which is also the length of each X_ij in the paper.
  #q is the length of the ture alpha, which is also the length of each Z_ij in the paper.
  #r is the length of the true gamma.
  #i is indicating which observation we are excluding from the 1st population when calculating MLE. The range of i is from 1 to N_1,
  #N_1 is the number of observations in 1st population.
  #The output is a list of length n_1, each element in this list is a vector of length N whose component is Y_ji_hat in formula in Step 3 in the paper.
  
  y1iexcludebfwall <- function(data.set,constbeta,p,q,r) {
    m=nrow(data.set)
    yall=vector("list",m)
    
    for (i in 1:m) {
      yall[[i]]=y1iexcludebfw(data.set,constbeta,p,q,r,i)
    }
    
    return (yall)
  }
  
  
  #This function cvw is used to calculate the crossvalidation criterion in step 4 of the paper.
  #The input data set should come from make.data.setforj_m_2nd function, when we are making the dataset, we only need to inlcude observations from
  #the 1st population.
  #The input constbeta is an N-1 length list whose component are from vector beta_2_hat to vector beta_N_hat because when calculating Y_1i_hat^(-i)(w), 
  #we don't need to re-calculate beta_2_hat to beta_N_hat, only beta_1_hat varies as we exclude i-th observation from the 1st population,
  #when calculate MLE of beta_1_hat,alpha_1_hat,gamma_1_hat.
  #The input w is an n dimensional vector of weights.
  #p is the length of the true beta, which is also the length of each X_ij in the paper.
  #q is the length of the ture alpha, which is also the length of each Z_ij in the paper.
  #r is the length of the true gamma.
  
  cvw <- function(data.set,constbeta,w,p,q,r) {
    n1=nrow(data.set)
    y1=data.set[,2]
    temp2=0
    temp1=0
    
    y1iexcludeall=y1iexcludebfwall(data.set,constbeta,p,q,r)
    
    y1inegihat=rep(0,n1)
    
    for (i in 1:n1) {
      y1inegihat[i]=sum(y1iexcludeall[[i]]*w)
    }
    
    for (i in 1:n1) {
      temp1=(y1inegihat[i]-y1[i])^2
      temp2=temp2+temp1
    }
    
    cvwout=temp2/n1
    
    return (cvwout)
  }
  
  
  
  #This function sim2 is used to combine previous subroutines to compute the crossvalidation criterion of step 4 of the paper.
  #The input sn is the number of simulation we would like to run, usually 1 is enough.
  #nhelper1 is the number of populations from the first sets of the helper models.
  #nhelper2 is the number of populations from the second sets of the helper models.
  #nhelper3 is the number of populations from the third sets of the helper models.
  #Since the main model is from population 1, so the total number of populations, N, is 1+nhelper1+nhelper2+nhelper3.
  #c2 is the number of observations from each population.
  #Because beta is assumed to be the same accross all N populations. 
  #So tbeta is a p dimensional vector, which is the true beta used to generate the data.
  #Because alpha is different for different populations, controlled by the index j.
  #So talpha is a list of length 1+nhelper1+nhelper2, whose elements are q dimensional vector, each q dimensional vector stands for true alpha_j used to generate the data.
  #The first element of this list of talpha is alpha vector for the main model.
  #The next nhelper1 elements of this list of talpha are alpha vectors for the 1st sets of helper models.
  #The next nhelper2 elements of this list of talpha are alpha vectors for the 2nd sets of helper models.
  #Because gamma is different for different populations, controlled by the index j. 
  #So tgamma is a list of length 1+nhelper1+nhelper2+nhelper3, whose elements are r dimensional vectors, each r dimensional vector stands for the true gamma_j used to generate the data.
  #The first element of this list of tgamma is gamma vector for the main model, so the length r of this vector can only be 1. 
  #The next nhelper1 elements of this list of tgamma are gamma vectors for the 1st sets of the helper models, the length r of this vector could be arbitrary.
  #The next nhelper2 elements of this list of tgamma are gamma vectors for the 2nd sets of the helper models, the length r of this vector can only be 1.
  #The next nhelper3 elements of this list of tgamma are gamma vectors for the 3rd sets of the helper models, the length r of this vector can only be 1.
  #Because delta is different for different populations, controlled by the index j.
  #So tdelta is a list of length 1+nhelper1+nhelper2+nhelper3, whose elements are p dimensional vectors, each p dimensional vector stands for the true delta_j used to generate the data.
  #The first element of this list of tdelta is delta vector for the main model, so the length p of this vector could be arbitrary. Actually, there is no delta in the main model.
  #The next nhelper1 elements of this list of tdelta are delta vectors for the 1st sets of the helper models.
  #The next nhelper2 elements of this list of tdelta are delta vectors for the 2nd sets of the helper models.
  #The next nhelper3 elements of this list of tdelta are delta vectors for the 3rd sets of the helper models.
  #each time we run the simulation, we need to modify the distribution part of x_mainm[[b]] and z_mainm[[b]] manually to accomendate the specific simulation settings. 
  
  sim2 <- function(sn,nhelper1,c2,tbeta,talpha,tgamma,tdelta) {
    c1=1+nhelper1
    p=length(tbeta)
    q=length(talpha[[1]])
    r1=length(tgamma[[1]])   #r1 stands for the length of the first element of tgamma list, which should actually always be 1.
    r2=length(tgamma[[2]])   #r2 stands for the length of the second element of tgamma list, which should be the length of u vector.
    lengthofzlast=length(tdelta[[1]])
    
    rho=0.5
    sd_m=rep(2,p+q+lengthofzlast)
    eps_sd_m=0.5
    eps_sd_1st=0.5
    
    
    
    si=1
    sim2outall=vector("list",sn)
    
    while (si<=sn) {
      print (si)
      
      #Generating list to store x,z data for the main model, in the main model we only have 1 population, so the number of element in these list is just 1.
      x_mainm=vector("list",1)
      z_mainm=vector("list",1)
      zlast_main=vector("list",1)
      
      xandz_mainm=vector("list",1)
      eps_mainm=vector("list",1)
      
      OMEGA=matrix(0,(p+q+lengthofzlast),(p+q+lengthofzlast))
      MU=rep(0,(p+q+lengthofzlast))
      
      for (i in 1:(p+q+lengthofzlast)) {
        for (j in 1:(p+q+lengthofzlast)){
          OMEGA[i,j]=rho^(abs(j-i))
        }
      }
      OMEGAcov=cor2cov(OMEGA,sd_m)
      
      for (b in 1:1) {
        xandz_mainm[[b]]=rmvnorm(c2[[1]][b],MU,OMEGAcov)
        x_mainm[[b]]=xandz_mainm[[b]][,1:p]
        z_mainm[[b]]=xandz_mainm[[b]][,(p+1):(p+q)]
        zlast_main[[b]]=xandz_mainm[[b]][,(p+q+1):(p+q+lengthofzlast)]
        
        eps_mainm[[b]]=matrix(rnorm(c2[[1]][b]*1,0,eps_sd_m),c2[[1]][b],1)
      }
      
      y_mainm=ydata_m(x_mainm,z_mainm,tbeta,talpha[[1]],tgamma[[1]],tdelta[[1]],zlast_main,eps_mainm[[1]])
      
      #Generating list to store x,z,u data for the 1st sets of the helper models, in the 1st sets of the helper models we have nhelper1 populations, so the number of elements in these lists
      #is nhelper1.
      x_1sthelperm=vector("list",nhelper1)
      z_1sthelperm=vector("list",nhelper1)
      zlast_1sthelperm=vector("list",nhelper1)
      
      u_1sthelperm=vector("list",nhelper1)
      xandz_1sthelperm=vector("list",nhelper1)
      eps_1sthelperm=vector("list",nhelper1)
      
      for (b in 1:nhelper1) {
        xandz_1sthelperm[[b]]=rmvnorm(c2[[2]][b],MU,OMEGAcov)
        x_1sthelperm[[b]]=xandz_1sthelperm[[b]][,1:p]
        z_1sthelperm[[b]]=xandz_1sthelperm[[b]][,(p+1):(p+q)]
        zlast_1sthelperm[[b]]=xandz_1sthelperm[[b]][,(p+q+1):(p+q+lengthofzlast)]
        
        u_1sthelperm[[b]]=matrix(runif(c2[[2]][b]*r2,1,1),c2[[2]][b],r2)
        eps_1sthelperm[[b]]=matrix(rnorm(c2[[2]][b]*1,0,eps_sd_1st),c2[[2]][b],1)
      }
      
      talpha_1st <- vector("list",nhelper1)
      tgamma_1st <- vector("list",nhelper1)
      tdelta_1st <- vector("list",nhelper1)
      
      for (i in 1:nhelper1) {
        talpha_1st[[i]]=talpha[[1+i]]
        tgamma_1st[[i]]=tgamma[[1+i]]
        tdelta_1st[[i]]=tdelta[[1+i]]
      }
      
      
      y_1sthelperm=vector("list",nhelper1)
      
      
      for (i in 1:(nhelper1/2)) {
        x_1sthelperm_temp=vector("list",1)
        z_1sthelperm_temp=vector("list",1)
        zlast_1sthelperm_temp=vector("list",1)
        u_1sthelperm_temp=vector("list",1)
        x_1sthelperm_temp[[1]]=x_1sthelperm[[i]]
        z_1sthelperm_temp[[1]]=z_1sthelperm[[i]]
        zlast_1sthelperm_temp[[1]]=zlast_1sthelperm[[i]]
        u_1sthelperm_temp[[1]]=u_1sthelperm[[i]]
        y_1sthelperm[[i]]=ydata_1st_mis(x_1sthelperm_temp,z_1sthelperm_temp,u_1sthelperm_temp,tbeta,talpha_1st[[i]],tgamma_1st[[i]],tdelta_1st[[i]],zlast_1sthelperm_temp,eps_1sthelperm[[i]])
      }
      
      for (i in (nhelper1/2+1):nhelper1) {
        x_1sthelperm_temp=vector("list",1)
        z_1sthelperm_temp=vector("list",1)
        u_1sthelperm_temp=vector("list",1)
        x_1sthelperm_temp[[1]]=x_1sthelperm[[i]]
        z_1sthelperm_temp[[1]]=z_1sthelperm[[i]]
        u_1sthelperm_temp[[1]]=u_1sthelperm[[i]]
        y_1sthelperm[[i]]=ydata_1st(x_1sthelperm_temp,z_1sthelperm_temp,u_1sthelperm_temp,tbeta,talpha_1st[[i]],tgamma_1st[[i]],tdelta_1st[[i]],eps_1sthelperm[[i]])
      }
      
      
      mainm_out_initial=vector("list",1)
      
      for (b in 1:1) {
        x_mainm_temp=x_mainm[[b]]
        z_mainm_temp=z_mainm[[b]]
        y_mainm_temp=y_mainm[b,]
        
        dataset_mainm_temp=make.data.setforj_m_2nd(c2[[1]][b],p,q,r1,y_mainm_temp,x_mainm_temp,z_mainm_temp)
        mainm_out_initial[[b]]=mainmodel(dataset_mainm_temp,p,q,r1)
      }
      
      #Extract alpah_1^hat and gamma_1^hat needed for step 2 of the paper.
      
      talpha1_hat=as.vector(mainm_out_initial[[1]]$out.fixed[,1][(p+1):(p+q)])
      tgamma1_hat=as.vector(mainm_out_initial[[1]]$out.fixed[,1][(p+q+1):(p+q+r1)])
      tbeta1_hat=as.vector(mainm_out_initial[[1]]$out.fixed[,1][1:p])
      
      firsthelperm_out_initial=vector("list",nhelper1)
      
      for (b in 1:nhelper1) {
        x_1sthelperm_temp=x_1sthelperm[[b]]
        z_1sthelperm_temp=z_1sthelperm[[b]]
        u_1sthelperm_temp=u_1sthelperm[[b]]
        y_1sthelperm_temp=y_1sthelperm[[b]]
        
        dataset_1sthelperm_temp=make.data.setforj_1st(c2[[2]][b],p,q,r2,y_1sthelperm_temp,x_1sthelperm_temp,z_1sthelperm_temp,u_1sthelperm_temp)
        firsthelperm_out_initial[[b]]=firsthelpmodel(dataset_1sthelperm_temp,p,q,r2)
      }
      
      #Extract alpha_j^hat and gamma_j^hat from the 1st sets of helper models.
      talpha_1st_hat=vector("list",nhelper1)
      tgamma_1st_hat=vector("list",nhelper1)
      
      for (b in 1:nhelper1) {
        talpha_1st_hat[[b]]=as.vector(firsthelperm_out_initial[[b]]$out.fixed[,1][(p+1):(p+q)])
        tgamma_1st_hat[[b]]=as.vector(firsthelperm_out_initial[[b]]$out.fixed[,1][(p+q+1):(p+q+r2)])
      }
      
      
      
      talpha_1st2nd3rd=vector("list",3)
      tgamma_1st2nd3rd=vector("list",3)
      
      talpha_1st2nd3rd[[1]]=talpha_1st_hat
      
      
      tgamma_1st2nd3rd[[1]]=tgamma_1st_hat
      
      
      #Generate constbeta input for cvwmatrix function.
      
      constbeta=vector("list",nhelper1)
      
      for (b in 1:nhelper1) {
        constbeta[[b]]=as.vector(firsthelperm_out_initial[[b]][[1]][,1][1:p])
      }
      
      
      
      constbetaall=vector("list",1+nhelper1)
      
      constbetaall[[1]]=tbeta1_hat
      
      for (b in 2:(1+nhelper1)) {
        constbetaall[[b]]=constbeta[[b-1]]
      }
      
      data_cvw=make.data.setforj_m_2nd(c2[[1]][1],p,q,r1,y_mainm[1,],x_mainm[[1]],z_mainm[[1]])
      
      #This cvwforsolnp is for solnp function to solve the optimal w.
      
      cvwforsolnp <- function(w) {
        cvwforsolnpout=cvw(data_cvw,constbeta,w,p,q,r1)
      }
      
      #constrain function of the weight vector w.
      
      eqcon <- function(w) {
        sumofw=sum(w)
        return (sumofw)
      }
      ebcon <-1
      
      #start values of solnp function.
      w0 <- rep(0.1,c1)
      
      sqp <- solnp(pars=w0,fun=cvwforsolnp,eqfun=eqcon,eqB=ebcon,LB=rep(0,c1),UB=rep(1,c1))
      
      
      
      #This function lltheta is used to calculate the log-likelihood function of all four models.
      #The input theta is a long vector.
      #The first p elements of this vector theta are the p components from beta parameter.
      #The next q*(1+nhelper1+nhelper2) elements of this vector theta are the q components from each of the alpha_j parameters for all models.
      #The next r1 element of this vector theta is the gamma_1 parameter of the main model, here r1 is actually always 1.
      #The next r2*nhelper1 elements of this vector theta are the r2 components from each of the gamma_j parameter of the 1st sets of the helper models.
      #The next r3*nhelper2 elements of this vector theta are the r3 components from each of the gamma_j parameter of the 2nd sets of the helper models.here r3 is actually always 1.
      #The next r4*nhelper3 elements of this vector theta are the r4 components from each of the gamma_j parameter of the 3rd sets of the helper models.here r4 is actually always 1.
      
      lltheta <- function(theta) {
        
        #We first need to restore the structure of beta,alpha,gamma from the input theta, as theta is a long vector.
        #betainll is used to store the beta parameters for all models, since beta doesn't vary across all populations, betainall is just a p dimensional vector.
        #alphainll_m is used to store the alpha_1 parameter for the main model, thus alphainll_m is a q dimensional vector.
        #alphainll_1st is used to store the alpha_j parameters for the 1st sets of the helper models, it is a list of length nhelper1, each element in this list is
        #a q dimensional vector representing alpha_j.
        #alphainll_2nd is used to store the alpha_j parameters for the 2nd sets of the helper models, it is a list of length nhelper2, each element in this list is
        #a q dimensional vector representing alpha_j.
        #gammainll_m is used to store the gamma_1 parameter for the main model, thus gammainll_m is a r1 dimensional vector, r1 is actually always 1.
        #gammainll_1st is used to store the gamma_j parameter for the 1st sets of the helper models, thus gammainll_1st is a list of length nhelper1,
        #each element in this list is a r2 dimensional vector representing gamma_j.
        #gammainll_2nd is used to store the gamma_j parameter for the 2nd sets of the helper models, thus gammainll_2nd is a list of length nhelper2,
        #each element in this list is a r3 dimensional vector representing gamma_j, r3 is actually always 1.
        #gammainll_3rd is used to store the gamma_j parameter for the 3rd sets of the helper models, thus gammainll_3rd is a list of length nhelper3,
        #each element in this list is a r4 dimensional vector representing gamma_j, r4 is actually always 1.
        betainll <- theta[1:p]
        alphainll <- theta[(p+1):(p+q*(1+nhelper1))]
        gammainll <- theta[(p+q*(1+nhelper1)+1):(p+q*(1+nhelper1)+r1*1+r2*nhelper1)]
        
        
        alphainll_m=alphainll[1:q]
        alphainll_1st=vector("list",nhelper1)
        
        
        gammainll_m=gammainll[1:r1]
        gammainll_1st=vector("list",nhelper1)
        
        
        
        for (b in 1:nhelper1) {
          alphainll_1st[[b]]=alphainll[(q*b+1):(q*b+q)]
          gammainll_1st[[b]]=gammainll[(r1*1+r2*b-r2+1):(r1*1+r2*b)]
        }
        
        
        py1ij=0
        for (b in 1:c2[[1]][1]) {
          py1ijtemp=(1/sqrt(2*pi*eps_sd_m^2))*exp(-(y_mainm[1,b]-(sum(x_mainm[[1]][b,]*betainll)+sum(z_mainm[[1]][b,]*alphainll_m)+sum(gammainll_m)))^2/(2*eps_sd_m^2))
          if (abs(py1ijtemp)<0.00001) {
            py1ij=py1ij+(-20)
          }
          else {py1ij=py1ij+log(py1ijtemp)}
          
        }
        
        py2ij=0
        
        for (c in 1:nhelper1) {
          for (b in 1:c2[[2]][c]) {
            py2ijtemp=(1/sqrt(2*pi*eps_sd_1st^2))*exp(-(y_1sthelperm[[c]][b]-(sum(x_1sthelperm[[c]][b,]*betainll)+sum(z_1sthelperm[[c]][b,]*alphainll_1st[[c]])+sum(u_1sthelperm[[c]][b,]*gammainll_1st[[c]])))^2/(2**eps_sd_1st^2))
            if (abs(py2ijtemp)<0.00001) {
              py2ij=py2ij+(-20)
            }
            else {py2ij=py2ij+log(py2ijtemp)}
          }
        }
        
        
        
        llthetaall=py1ij+py2ij
        
        return (-llthetaall)
      }
      
      
      #The following code is to construct the suitale initial value thetainitial for lltheta function. We will use the input tbeta,talpha,tgamma to construct.
      
      thetainitial=tbeta
      
      for (b in 1:(1+nhelper1)) {
        thetainitial=c(thetainitial,talpha[[b]])
      }
      
      for (b in 1:(1+nhelper1)) {
        thetainitial=c(thetainitial,tgamma[[b]])
      }
      
      #Use bobyqa function to minimize -lltheta function, which is equivalent to maximize  
      thetaout=auglag(thetainitial,lltheta)
      
      sim2out=vector("list",12)
      
      sim2out[[1]]=sqp
      sim2out[[2]]=constbetaall
      sim2out[[3]]=constbeta
      sim2out[[4]]=talpha1_hat
      sim2out[[5]]=tgamma1_hat
      sim2out[[6]]=thetaout
      sim2out[[7]]=talpha_1st2nd3rd
      sim2out[[8]]=tgamma_1st2nd3rd
      sim2out[[9]]=y_mainm
      sim2out[[10]]=y_1sthelperm
      #sim2out[[11]]=y_2ndhelperm
      #sim2out[[12]]=y_3rdhelperm
      #sim2out[[13]]=y_all
      
      
      sim2outall[[si]]=sim2out
      si=si+1
      
    }    
    return (sim2outall)
    
  }
  
  #This function testdata_m is used to generate data from the first population to test our estimator, the simple average method and MLE method.
  #The input n is the number of obeservations we want to generate from the 1st populations.
  #The input tbeta is a p dimensional vector stands for true beta in the main model.
  #The input talpha is a list of length 1, the element in this list is a q dimensional vector stands for true alpha in the main model
  #The input tgamma is a list of length 1, the element in this list is an r dimensional vector stands for true gamma in the main model.
  #Note that p,q,r should be consisitent with the data used for sim2 functions.
  #If we want to change the distribution of x and z, we need to modify this function testdata_m.
  #The output is a list of length 3.
  #the 1st element in this list is a n*1 matrix, standing for the n observations of length 1 expectation of outcome y.
  #The 2nd element in this list is a n*p matrix, standing for the n observations of length p vector x.
  #the 3rd element in this list is a n*q matrix, standing for the n observations of length q vector z.
  
  testdata_m <- function(n,tbeta,talpha,tgamma,tdelta) {
    p=length(tbeta)
    q=length(talpha[[1]])
    r=length(tgamma)
    eps_sd_test=0.5
    lengthofzlast=length(tdelta[[1]])
    
    rho=0.5
    sd_m=rep(2,p+q+lengthofzlast)
    OMEGA=matrix(0,p+q+lengthofzlast,p+q+lengthofzlast)
    MU=rep(0,p+q+lengthofzlast)
    
    for (i in 1:(p+q+lengthofzlast)) {
      for (j in 1:(p+q+lengthofzlast)){
        OMEGA[i,j]=rho^(abs(j-i))
      }
    }
    OMEGAcov=cor2cov(OMEGA,sd_m)
    
    x_test=vector("list",1)
    z_test=vector("list",1)
    zlast_test=vector("list",1)
    xandz_test=vector("list",1)
    eps_test=vector("list",1)
    
    for (b in 1:1) {
      xandz_test[[b]]=rmvnorm(n,MU,OMEGAcov)
      x_test[[b]]=xandz_test[[b]][,1:p]
      z_test[[b]]=xandz_test[[b]][,(p+1):(p+q)]
      zlast_test[[b]]=xandz_test[[b]][,(p+q+1):(p+q+lengthofzlast)]
      eps_test[[b]]=matrix(rnorm(n*1,0,eps_sd_test),n,1)
    }
    
    y_test=ydata_m_exp(x_test,z_test,tbeta,talpha[[1]],tgamma[[1]],tdelta[[1]],zlast_test)
    
    testdata_out=vector("list",3)
    
    testdata_out[[1]]=y_test
    testdata_out[[2]]=x_test
    testdata_out[[3]]=z_test
    
    return (testdata_out)
  }
  
  #This function extest is used to extract all necessary information from the output of sim2 function.
  #The input sim2out is from the output of sim2 functions.
  
  extest <- function(sim2out) {
    sn=length(sim2out)  #this is how many simulations we sotre in sim2out
    dataout=vector("list",sn)
    p=length(sim2out[[1]][[2]][[1]])
    
    for (b in 1:sn) {
      dataout[[b]]=vector("list",8)
      dataout[[b]][[1]]=sim2out[[b]][[1]]$pars   #extract the optimal weight w.
      dataout[[b]][[2]]=sim2out[[b]][[2]][[1]]   #extract the MLE of beta just using the main model, which is a p dimensional vector beta_1_hat.
      dataout[[b]][[3]]=sim2out[[b]][[2]]        #extract the MLE of all beta_j using resepective models, which are p dimensional vectors beta_j_hat.
      dataout[[b]][[4]]=sim2out[[b]][[4]]        #extract the MLE of alpha from the main model, which is a q dimensional vector alpha_1_hat.
      dataout[[b]][[5]]=sim2out[[b]][[5]]        #extract the MLE of gamma from the main model, which is a r dimensional vector gamma_1_hat.
      dataout[[b]][[6]]=sim2out[[b]][[6]]$par[1:p] #extract the MLE of beta using all models, which is a p dimensional vector beta_1_hat_all.
      dataout[[b]][[7]]=sim2out[[b]][[7]]
      dataout[[b]][[8]]=sim2out[[b]][[8]]
    }
    
    return (dataout)
  }
  
  #This function ourest_i is used to calculate the mean square loss of Y_hat from our estimator, using the i-th simulation results.
  #The input testdata is from the output of testdata_m function.
  #The input estpar is from the i-th output of extest function.
  
  ourest_i <- function(testdata,estpar){
    n=length(estpar[[1]])  #This n stands for N in the paper, which is the length of w vector.
    m=nrow(testdata[[2]][[1]])  #This m stands for how many observations we have in this test data set.
    allbeta=estpar[[3]]
    alpha_1=estpar[[4]]
    gamma_1=estpar[[5]]
    w1=estpar[[1]]
    ypre=rep(0,m)
    yreal=testdata[[1]]
    
    for (b in 1:m) {
      x_temp=testdata[[2]][[1]][b,]
      z_temp=testdata[[3]][[1]][b,]
      yjpretemp=rep(0,n)
      for (c in 1:n) {
        yjpretemp[c]=yjpre(x_temp,z_temp,allbeta[[c]],alpha_1,gamma_1)
      }
      ypre[b]=sum(w1*yjpretemp)
    }
    
    msl=0
    
    for (b in 1:m) {
      msltemp=(ypre[b]-yreal[1,b])^2
      msl=msl+msltemp
    }
    
    return (msl/m)
  }
  
  
  ourest_i_y_out <- function(testdata,estpar){
    n=length(estpar[[1]])  #This n stands for N in the paper, which is the length of w vector.
    m=nrow(testdata[[2]][[1]])  #This m stands for how many observations we have in this test data set.
    allbeta=estpar[[3]]
    alpha_1=estpar[[4]]
    gamma_1=estpar[[5]]
    w1=estpar[[1]]
    ypre=vector("list",m)
    yreal=testdata[[1]]
    
    for (b in 1:m) {
      x_temp=testdata[[2]][[1]][b,]
      z_temp=testdata[[3]][[1]][b,]
      yjpretemp=rep(0,n)
      for (c in 1:n) {
        yjpretemp[c]=yjpre(x_temp,z_temp,allbeta[[c]],alpha_1,gamma_1)
      }
      ypre[[b]]=yjpretemp
    }
    
    return (ypre)
  }
  
  
  #This function sam_i is used to predit the Y_i_hat from the simple average method, using the i-th simulation results.
  #The input testdata is from the output of testdata_m function.
  #The input estpar is from the i-th output of extest function.
  
  sam_i <- function(testdata,estpar) {
    n=length(estpar[[1]])  #This n stands for N in the paper, which is the length of w vector.
    m=nrow(testdata[[2]][[1]])  #This m stands for how many observations we have in this test data set.
    allbeta=estpar[[3]]
    alpha_1=estpar[[4]]
    gamma_1=estpar[[5]]
    w1=rep(1/n,n)
    ypre=rep(0,m)
    yreal=testdata[[1]]
    
    for (b in 1:m) {
      x_temp=testdata[[2]][[1]][b,]
      z_temp=testdata[[3]][[1]][b,]
      yjpretemp=rep(0,n)
      for (c in 1:n) {
        yjpretemp[c]=yjpre(x_temp,z_temp,allbeta[[c]],alpha_1,gamma_1)
      }
      ypre[b]=sum(w1*yjpretemp)
    }
    msl=0
    for (b in 1:m) {
      msltemp=(ypre[b]-yreal[1,b])^2
      msl=msl+msltemp
    }
    return (msl/m)
  }
  
  #This function mlemain_i is used to predit the Y_i_hat from the MLE method using MLE of beta from the main model, using the i-th simulation results.
  #The input testdata is from the output of testdata_m function.
  #The input estpar is from the i-th output of extest function.
  
  mlemain_i <- function(testdata,estpar) {
    n=length(estpar[[1]])  #This n stands for N in the paper, which is the length of w vector.
    m=nrow(testdata[[2]][[1]])  #This m stands for how many observations we have in this test data set.
    allbeta=estpar[[3]]
    beta_1=estpar[[2]]
    alpha_1=estpar[[4]]
    gamma_1=estpar[[5]]
    w1=estpar[[1]]
    ypre=rep(0,m)
    yreal=testdata[[1]]
    
    for (b in 1:m) {
      x_temp=testdata[[2]][[1]][b,]
      z_temp=testdata[[3]][[1]][b,]
      ypre[b]=yjpre(x_temp,z_temp,beta_1,alpha_1,gamma_1)
    }
    
    msl=0
    
    for (b in 1:m) {
      msltemp=(ypre[b]-yreal[1,b])^2
      msl=msl+msltemp
    }
    
    return (msl/m)
  }
  
  #This function mleall_i is used to predit the Y_i_hat from the MLE method using MLE of beta from all models, using the i-th simulation results.
  #The input testdata is from the output of testdata_m function.
  #The input estpar is from the i-th output of extest function.
  
  mleall_i <- function(testdata,estpar) {
    n=length(estpar[[1]])  #This n stands for N in the paper, which is the length of w vector.
    m=nrow(testdata[[2]][[1]])  #This m stands for how many observations we have in this test data set.
    allbeta=estpar[[3]]
    beta_1=estpar[[6]]
    alpha_1=estpar[[4]]
    gamma_1=estpar[[5]]
    w1=estpar[[1]]
    ypre=rep(0,m)
    yreal=testdata[[1]]
    
    for (b in 1:m) {
      x_temp=testdata[[2]][[1]][b,]
      z_temp=testdata[[3]][[1]][b,]
      ypre[b]=yjpre(x_temp,z_temp,beta_1,alpha_1,gamma_1)
    }
    
    msl=0
    
    for (b in 1:m) {
      msltemp=(ypre[b]-yreal[1,b])^2
      msl=msl+msltemp
    }
    
    return (msl/m)
  }
  
  #This function tauhat is used to obtain the sum of the weight distributed to the correct models.
  #The input estpar is from the first element of the list of the output of extest function.
  
  tauhat <- function(estpar) {
    w1=estpar[[1]]
    
    tauhat_temp=sum(w1[1])+sum(w1[(1+nhelper1/2+1):(1+nhelper1)])
    
    return (tauhat_temp)
  }
  
  #Parameters estimating part.
  
  nhelper1=2                                        #nhelper1 is the number of populations from the 1st sets of the helper models.
  nhelper2=2                                        #nhelper2 is the number of populations from the 2nd sets of the helper models.
  nhelper3=2                                        #nhelper3 is the number of populations from the 3rd sets of the helper models.
  c2=vector("list",2)
  c2[[1]]=c(500)                                    #c2 is the number of observations from each population.
  c2[[2]]=c(400,300)
  
  
  tbeta=c(0.5,0.6,-0.61,-0.48)
  talpha=vector("list",1+nhelper1+nhelper2)          #a list of length 1+nhelper1+nhelper2 to store alpha_j for different models.
  tgamma=vector("list",1+nhelper1+nhelper2+nhelper3) #a list of length 1+nhelper1+nhelper2+nhelper3 to store gamma_j for different models.
  tdelta=vector("list",1+nhelper1+nhelper2+nhelper3) #a list of length 1+nhelper1+nhelper2+nhelper3 to store delta_j for different models.
  
  
  #parameters for the main model.
  talpha[[1]]=c(0.6,0.5,-0.30,-0.25)
  tgamma[[1]]=0.4  
  tdelta[[1]]=0.1
  
  #parameters for the 1 of 1st helper model.
  talpha[[2]]=c(0.08,0.09,-0.04,-0.06) 
  tgamma[[2]]=0.49
  tdelta[[2]]=2.5
  
  #parameters for the 2 of 1st helper model.
  talpha[[3]]=c(0.07,0.1,-0.05,-0.04) 
  tgamma[[3]]=0.51
  tdelta[[3]]=0
  
  
  simu201804225=sim2(1,nhelper1,c2,tbeta,talpha,tgamma,tdelta) #the 1st argument means how many simulations we would like to run.
  
  #Test data generating part.
  
 
  tbeta_real=tbeta                                   #set the beta used to generate the test data, same beta in the main model.
  talpha_real=vector("list",1)
  talpha_real[[1]]=talpha[[1]]                       #set the alpha used to generate the test data, same alpha in the main model.
  tgamma_real=vector("list",1)
  tgamma_real[[1]]=tgamma[[1]]                       #set the gamma used to generate the test data, same gamma in the main model.
  tdelta_real=vector("list",1)
  tdelta_real[[1]]=tdelta[[1]]
  
  testdata_1=testdata_m(500,tbeta_real,talpha_real,tgamma_real,tdelta_real) #the 1st argument means how many observations we would like in this test data from the main model.
  estpar1=extest(simu201804225)                      #extract the estimated parameters used in our method, the simple average method, the MLE using all models, and the MLE using only the main model for our test data, 
  msltest1=ourest_i(testdata_1,estpar1[[1]])         #The number of list indicator in the 2nd argument stands for which simulation result we are using for the test data.
  msltest2=sam_i(testdata_1,estpar1[[1]])            #The number of list indicator in the 2nd argument stands for which simulation result we are using for the test data.
  msltest3=mlemain_i(testdata_1,estpar1[[1]])        #The number of list indicator in the 2nd argument stands for which simulation result we are using for the test data.
  msltest4=mleall_i(testdata_1,estpar1[[1]])         #The number of list indicator in the 2nd argument stands for which simulation result we are using for the test data.
  y_out_observe=ourest_i_y_out(testdata_1,estpar1[[1]])
  tau1=tauhat(estpar1[[1]])
  tau1
  msltest1
  msltest2
  msltest3
  msltest4
  
  msltemp=rep(0,5)
  msltemp[1]=msltest1
  msltemp[2]=msltest2
  msltemp[3]=msltest3
  msltemp[4]=msltest4
  msltemp[5]=tau1
  
  
  s1.2[[indexi]]=msltemp
  s1.2y_m[[indexi]]=simu201804225[[1]][[9]]
  s1.2y_1st[[indexi]]=simu201804225[[1]][[10]]
  s1.2y_2nd[[indexi]]=simu201804225[[1]][[11]]
  s1.2y_3rd[[indexi]]=simu201804225[[1]][[12]]
  s1.2estpar[[indexi]]=estpar1
  s1.2y_ob_out[[indexi]]=y_out_observe
  s1.2testdata[[indexi]]=testdata_1
  s1.2estpar[[indexi]]=estpar1
  indexi=indexi+1
}

s1.2msl1=rep(0,100)

for (i in 1:100) {
  s1.2msl1[i]=s1.2[[i]][1]
}

s1.2msl2=rep(0,100)

for (i in 1:100) {
  s1.2msl2[i]=s1.2[[i]][2]
}

s1.2msl3=rep(0,100)

for (i in 1:100) {
  s1.2msl3[i]=s1.2[[i]][3]
}

s1.2msl4=rep(0,100)

for (i in 1:100) {
  s1.2msl4[i]=s1.2[[i]][4]
}

s1.2msl5=rep(0,100)

for (i in 1:100) {
  s1.2msl5[i]=s1.2[[i]][5]
}


which(is.nan(s1.2msl1),arr.ind = TRUE)
which(is.nan(s1.2msl2),arr.ind = TRUE)
which(is.nan(s1.2msl3),arr.ind = TRUE)
which(is.nan(s1.2msl4),arr.ind = TRUE)
which(is.nan(s1.2msl5),arr.ind = TRUE)

msl1aver=sum(s1.2msl1)/100
msl2aver=sum(s1.2msl2)/100
msl3aver=sum(s1.2msl3)/100
msl4aver=sum(s1.2msl4)/100
msl5aver=sum(s1.2msl5)/100
msl1aver
msl2aver
msl3aver
msl4aver
msl5aver