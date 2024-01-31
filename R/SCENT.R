#' SCENT
#'
#' @description This function uses the ADMM algorithm to conduct Simultaneous Clustering and Estimation of Networks
#'
#' @param SS A data array of dimension p by p by K of covariance matrices
#' @param r A scalar of the prespecified rank of the tensor decomposition (an upper bound for # clusters)
#' @param nn A vector of length K of sample size for each population
#' @param rho A scalar of augmented parameter in ADMM, default=0.1
#' @param Tol Convergence tolerance for relative residuals, default=1E-4
#' @param Niter Maximum number of iterations, default=200
#'
#' @returns primal, dual, obj, pen: each being a sequence of values corresponding to
#'      the primal residual, dual residual, objective value, and penalty value
#' @returns Omega A data array of dimension p by p by K of precision matrices
#' @returns Theta0 A matrix of size p by p, common network structure
#' @returns Theta A tensor of size p by p by r (or a smaller rank, can be 0), each layer corresponds to a unit-norm, sparse sub-network
#' @returns U A matrix of size K by r (or a smaller rank, can be 0), the nonsparsity of each column corresponds to the cluster index
#'
#' @examples
#' SS=rWishart(6,20,0.01*diag(20))
#' out=SCENT(SS,r=2,nn=rep(100,6),
#'               rho=50)
#' par(mfrow=c(1,2))
#' plot(out$primal)
#' plot(out$dual)
#' out$U
#'
#' @import MASS
#' @keywords SCENT
#' @export



SCENT=function(SS,r,nn,
               rho=1,  Tol=1E-4,  Niter=200){

  ## initialization
  K=dim(SS)[3]
  p=dim(SS)[1]
  #
  Lambda=array(0,dim=dim(SS)) # lagrange multiplier
  Omega=array(0,dim=dim(SS)) # precision matrix
  Theta=array(0,dim=c(p,p,r)) # subgraph matrix
  U=matrix(0,nrow=K,ncol=r) # loading/cluster index
  for(k in 1:K){
    Omega[,,k]=ginv(SS[,,k]+0.1*mean(diag(SS[,,k]))*diag(p) ) # 1. to make sure initial omega is not too off
  }
  #
  Theta0=apply(Omega, c(1,2), mean)
  deOmega=sweep(Omega,c(1,2),Theta0)
  # apply Tucker1 to get initial est of U and Theta
  out=Tucker1(deOmega,r)
  U=out$loading # K*r unit norm loadings
  Theta=out$tensor # p*p*r tensor of subnetworks
  #
  obj.ini=objective(SS,Omega,nn)




  ## ADMM
  primal=dual=obj=pen=NULL
  niter=0
  diff=Inf
  curr.r=r
  while(niter<=Niter & diff>Tol){
    niter=niter+1
    Omega.old=Omega # prev precision tensor
    #print(paste('iter',niter))

    ### step 1 update Omega
    for(k in 1:K){
      Omegak=Theta0+apply(sweep(Theta,MARGIN=3,U[k,],'*'),MARGIN=c(1,2),sum)
      Mat=SS[,,k]-rho/nn[k]*Omegak+Lambda[,,k]/nn[k] ## more stable
      Q=eigen(Mat)$vector
      D=eigen(Mat)$value
      value=2/(D+sqrt(D^2+4*rho/nn[k]))
      if(abs(prod(D))<=0.001){ ## degenerate
        value[abs(D)<=0.1]=0 ## how to choose cutoff value
      }
      Omega[,,k]=Q%*%diag(value)%*%t(Q)
    }
    #print(paste('Omega OK'))

    ### Step 2 update Theta0
    Theta0=apply(Omega+Lambda/rho, c(1,2), mean) # non-sparse est
    # soft thresholding
    lam.min=min(abs(Theta0))
    lam.max=max(abs(Theta0-diag(diag(Theta0))))
    seg=20
    lam.cand=lam.min+((1:seg)/seg)*(lam.max-lam.min)
    Theta0.BIC=array(Inf,dim = c(1,seg))
    for (ind in 1:seg) {
      Theta0.curr=SoftThres(Theta0,lam.cand[ind])
      Theta0.BIC[ind]=log(sum((Theta0-Theta0.curr)^2)/p^2)+
        log(p^2)*sum(Theta0.curr!=0)/p^2
    }
    Theta0=SoftThres(Theta0,lam.cand[which.min(Theta0.BIC)]) # sparse
    #print(paste('Theta0 OK'))


    ### Step 3 update Theta and U
    fulltensor=sweep(Omega+Lambda/rho, MARGIN=c(1,2), Theta0,'-')
    Theta=array(0,dim=dim(Theta)) # reset
    U=array(0,dim=dim(U)) # reset
    keep=1:K
    for (l in 1:r){ # deflation
      subtensor=fulltensor[,,keep]
      out=Tucker1(subtensor,1)
      u.temp=out$loading
      theta.temp=out$tensor[,,1]
      #print(paste('initial OK'))

      # soft thresholding
      lam.min=min(abs(theta.temp))
      lam.max=max(abs(theta.temp))
      seg=100
      lam.cand=lam.min+(1:seg)*(lam.max-lam.min)/(seg+1)
      theta.BIC=array(Inf,dim=c(length(keep),seg))
      theta.BIC[1,]=log( sum(subtensor^2)/prod(dim(subtensor)) ) # all zero est



      for (ind.u in 2:length(keep)) {
        for (ind.theta in 1:seg) {
          # threshold u.temp
          out=TakeTop(u.temp,ind.u)
          u.curr=out$newvec

          # threshold theta.temp
          theta.curr=SoftThres(theta.temp,lam.cand[ind.theta])

          # BIC
          curr.tensor=array(0,dim=c(p,p,length(keep)))
          for (temp.ind in 1:length(keep)) {
            curr.tensor[,,temp.ind]=u.curr[temp.ind]*theta.curr
          }
          theta.BIC[ind.u,ind.theta]=log( sum((subtensor-curr.tensor)^2)/prod(dim(subtensor)) )+
            log(K*p^2)*(ind.u*sum(theta.curr!=0))/(K*p^2)
        }
      }
      #print(paste('thres OK'))

      out=which(theta.BIC==min(theta.BIC),arr.ind = TRUE)
      #print(out)

      ind.u=out[1]
      ind.theta=out[2]
      if (ind.u==1){ # all zero u is best
        curr.r=l-1 # current rank
        break
      }
      # otherwise
      curr.r=l # current rank
      out=TakeTop(u.temp,ind.u)
      u.final=out$newvec
      U[keep,l]=u.final
      keep=keep[out$unsel] # a subvec of 1:K with candidate indices for next layer

      #
      theta.final=SoftThres(theta.temp, lam.cand[ind.theta])
      Theta[,,l]=theta.final

      #
      if (length(keep)<2){ # no need to continue
        break
      }

    }
    #print(paste('Theta and U OK'))
    #print(paste('iter',niter,', Selected rank:',curr.r))

    ### step 4 Dual step
    # construct Tucker1 tensor
    CC=sweep(Product1(U,Theta),MARGIN = c(1,2),Theta0,'+')
    Lambda=Lambda+rho*(Omega-CC)
    #print(paste('Dual OK'))


    ### Residuals and Stopping Criteria
    primal=c(primal,sqrt(sum((Omega-CC)^2))/sqrt(sum(Omega^2)))
    dual=c(dual,sqrt(sum(Omega-Omega.old)^2)/sqrt(sum(Omega.old^2)))
    obj=c(obj,objective(SS,Omega,nn))
    pen=c(pen,sum(Lambda*(Omega-CC))+rho/2*sum((Omega-CC)^2))
    #
    diff=primal[length(primal)] # use relative primal residual as stopping rule
  }


  ### normalize final estimate
  if(curr.r<r){
    print(paste("The final Tucker 1 rank is",curr.r,"(lower than the given rank",r,")!"))
  }
  if(curr.r==0){
    U=array(0,dim = c(K,0))
    Theta=array(0,dim=c(p,p,0))
  } else {
    U=U[,1:curr.r,drop=FALSE]
    Theta=Theta[,,1:curr.r,drop=FALSE]
    for(ind in 1:curr.r){
      norm=sqrt(sum(Theta[,,ind]^2))
      Theta[,,ind]=Theta[,,ind]/norm
      U[,ind]=U[,ind]*norm
    }
  }


  ### convergence
  if(niter<Niter){
    print(paste("SCENT converges in",niter,"steps."))
  } else {
    print(paste("SCENT does NOT converge within",Niter,"steps! Final primal residual is",diff))
  }


  ###
  return(list("primal"=primal,"dual"=dual,"obj"=obj,"pen"=pen,
              "Omega"=Omega,"Theta0"=Theta0,"U"=U,"Theta"=Theta,
              "CC"=CC,"Lambda"=Lambda))
}



Tucker1=function(TT,r){
  # this function applies rank-r Tucker1 decomposition to the 3rd dimension of TT
  mat=matrix(c(TT),nrow=prod(dim(TT)[c(1,2)]),ncol=dim(TT)[3])
  out=svd(mat)
  if(r==1){
    loading=out$v[,1:r,drop=FALSE]
    tensor=array(out$u[,1]*out$d[1],dim=c(dim(TT)[1:2],r))
  } else {
    loading=out$v[,1:r]
    tensor=array(c(out$u[,1:r]%*%diag(out$d[1:r])),dim=c(dim(TT)[1:2],r))
  }
  return(list('loading'=loading,'tensor'=tensor))
}


Product1=function(mat,tensor){ # product 1 of tensor
  # mat: K*r
  # tensor: p*p*r
  # out: p*p*K
  out=array(0,dim = c(dim(tensor)[1:2],dim(mat)[1]))
  for (j in 1:dim(mat)[1]) {
    out[,,j]=apply(sweep(tensor,MARGIN=3,mat[j,],'*'),MARGIN=c(1,2),sum)
  }
  return(out)
}

objective=function(SS,Omega,n){
  K=dim(SS)[3]
  obj=0
  for(k in 1:K){
    ## might not be p.s.d.
    if(det(Omega[,,k])>0){
      obj=obj+n[k]*(sum(SS[,,k]*Omega[,,k])-log(det(Omega[,,k])))
    } else obj=obj+n[k]*sum(SS[,,k]*Omega[,,k])
  }
  return(obj)
}


TakeTop=function(vec,k){
  # take the top k (abs) value of column vector v
  out=sort.int(abs(vec),decreasing = TRUE,index.return = TRUE)
  ind=out$ix[1:k]
  unind=setdiff(1:nrow(vec),ind)
  temp.v=scale(vec[ind],center=TRUE,scale=FALSE)
  temp.v=temp.v/sqrt(sum(temp.v^2)) # unit norm
  newvec=array(0,dim=dim(vec))
  newvec[ind]=temp.v
  return(list('unsel'=unind,'newvec'=newvec))
}

SoftThres=function(Theta,lam){
  pos=ifelse(Theta > lam, Theta-lam, 0)
  neg=ifelse(Theta < (-lam),Theta+lam, 0)
  return(pos+neg)
}








