library(nlme);library(lme4)
lmer(score ~ Machine + (1|Worker) + (1|Worker:Machine),data=Machines,
     REML=FALSE)
lmm = function(form=0,dat,ref=list()){
  LMMsetup = function(){
    if(length(ref)>0){
    z = matrix(nrow = dim(dat)[1])
    facs = c() #collect number of factors for construction of psi
    for (col in ref){
      s = col[1]
      sc = length(levels(dat[,s]))
      if (length(col)>1){
        for (j in 2:length(col)){
          s = paste(s,':',col[j])
          sc = sc*length(levels(dat[,col[j]]))
        }
      }
      z = cbind(z,model.matrix(formula(paste('~',s,'-1')),dat))
      facs = append(facs,sc)
    }
    z = z[,2:ncol(z)]
    x = model.matrix(formula(form),dat)
    return(list(z=z,x=x,facs=facs))
  }
  }

  f = LMMsetup()
  z = f$z
  x = f$x
  facs = f$facs # for list(worker, c(worker,machine)), facs = c(6,18)
  qrz = qr(z)
  y_col = strsplit(format(form),' ~ ')[[1]][1]
  y = dat[,y_col]
  betas = c()
  thetas = c()
  
  LMMprof = function(theta, y, x, qrz, facs){
    sigma_sq = (exp(theta[1]))**2 #first element is log sigma
    
    n = dim(x)[1]
    p = sum(facs)
    psi = c(rep(0,p))
    start = 1
    for (j in 1:length(facs)){
      psi[start:(start+facs[j]-1)] = (exp(theta[1+j]))**2 #filling in variances of random effects
      start = facs[j]+start
    }
    
    r = qr.R(qrz)
    cp = chol(diag(psi)) #psi = (cp^T)(cp)
    t = r%*%cp
    rpr = t%*%t(t) #R(psi)R^T, i think this part looks a little weird
    mid_v = c(rep(sigma_sq,n))
    mid = diag(mid_v)
    mid[1:p,1:p] = mid[1:p,1:p] + rpr #R(psi)R^T + I(sigma_sq), mid is middle part of W without Qs
    lg_zpz = log(diag(mid)[1:p]) + (n-p)*log(sigma_sq)
    
    m = chol(mid)
    mid_inv = backsolve(m,forwardsolve(t(m),diag(1,nrow = n)))
    xtwy = t(x)%*%qr.qy(qrz,mid_inv)%*%qr.qty(qrz,y)
    xtwx = t(x)%*%qr.qy(qrz,mid_inv)%*%qr.qty(qrz,x)
    cr = chol(xtwx)
    beta_hat = backsolve(cr,forwardsolve(t(cr),xtwy))
    betas = append(betas,beta_hat)
    
    y_xb = y-(x%*%beta_hat)
    l = t(y_xb)%*%qr.qy(qrz,chol(mid_inv))
    -0.5*sum(l^2)-sum(lg_zpz)*0.5
  }
  start = c(3,3,3)
  mle_theta = capture.output(optim(start, LMMprof, y=y, x=x, qrz=qrz, facs=facs, control=list(trace=TRUE)))
  theta = mle_theta$par[1]
  print(theta)
  thetas = append(thetas,theta)
  return(list(thetas = thetas, betas = betas))
}


d = lmm(score~Machine,Machines,list("Worker",c("Worker","Machine")))
