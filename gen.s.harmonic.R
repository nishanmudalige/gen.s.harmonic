# Created by Nishan Mudalige and Dr. Peter Kim
# Maintained by Nishan Mudalige. Email: nishan.mudalige.1@ulaval.ca
# Copyright Nishan Mudalige and Dr. Peter Kim, 2020

# This research us partially based upon work supported by a research grant 
# from Google Cloud.

library(orthopolynom)
library(pracma)
library(matrixcalc)
library(mefa)
library(glmnet)
library(Directional)
library(hypergeo)

# Run function in simdd
library(matrixcalc)
library(stringr)
library(assertr)

# For simulations
library(Matrix)
library(pracma)
library(assertr)


rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}




b.q = function(q){
  p = length(q)
  ifelse(q[p-1] + q[p] > 0, 2, 1)
}



g.x = function(x, q){
  p = length(q)
  q.p = q[p]
  q.p_1 = q[p-1]
  
  if(q.p == 0) {
    result = Re(( complex(real=x[2], imaginary=x[1]) )^q.p_1)
    return(Re(result))
  } else if(q.p == 1) {
    result = Im(( complex(real=x[2], imaginary=x[1]) )^q.p_1)
    return(Re(result))
  } else {
    return(NULL) 
  }
  
}



q.norm = function(q, j){
  p_1 = length(q)-1
  return(sum(q[j:p_1]))
}



lambda = function(q, j){
  p = length(q)
  q.norm(q, j+1) + (p-j-1)/2
}


gegenbauer = function(n, lambda, x){
  
  res = ( (orthopolynom::pochhammer(2*lambda, n))/factorial(n) ) * 
    ( hypergeo::hypergeo(tol=0.1, A=-n, B=n+(2*lambda), C=lambda+0.5, z=(1-x)/2) )
  
  return(Re(res))
}


# # gegenbauer.old
# gegenbauer = function(n, lambda, x){
# 
#   if(n == 0){
#     result = 1
#     return(result)
#   } else if(n == 1){
#     result = 2*lambda*x
#     return(result)
#   } else if(n == 2){
#     result = -lambda + 2*lambda*(1 + lambda)*x^2
#     return(result)
#   } else {
#     warning("This function only supports up to second order terms")
#     return()
#   }
# 
# }


h.q.sqr = function(q){
  p = length(q)
  bq = b.q(q)
  result = NULL
  
  for(j in 1:(p-2)){
    lambda.j = lambda(q,j)
    num = factorial(q[j]) * pochhammer( (p-j+1)/2 , q.norm(q, j+1) ) * ( q[j] + lambda.j )
    denom = pochhammer( (2*lambda.j), q[j] ) * pochhammer( (p-j)/2 , q.norm(q, j+1) ) * lambda.j
    
    result[j] = num/denom
  }
  
  return(bq* prod(unlist(result)) )  
  
}


volume.hyper.sphere = function(p, r=1){
  # volume = ( pi^(p/2) / gamma(p/2 + 1) )*(r^p)
  volume = 2*pi^(p/2) / gamma(p/2) * (r^p)
  return(volume)
}




y.q.x = function(x, q, normalize.by.volume = T){
  
  p = as.numeric(length(x))
  
  gx = g.x(x, q)
  h.q = sqrt(h.q.sqr(q))
  
  result = NULL
  
  for(j in 1:(p-2)){
    #j = 2
    x.ss = sum( ( x[1:(p-j+1)] )^2 )
    lambda.j = lambda(q,j)
    
    g.poly = NULL
    
    if(q[j] == 0){
      g.poly = 1  
    } else {
      g.poly = gegenbauer( n=q[j], lambda = lambda.j, x =  as.numeric(x[p-j+1]/sqrt(x.ss)) )
    }
    
    result[j] =   x.ss^(q[j]/2) * g.poly
    
    
  }
  
  if(isTRUE(normalize.by.volume)){
    
    rad = as.numeric(sum(x^2))
    
    volume = volume.hyper.sphere(p, rad)
    
    return( (1/sqrt(volume)) * h.q * gx * prod(unlist(result)) )
  } else{
    return( h.q * gx * prod(unlist(result)) )
  }
  
}



y.q.x(c(0,0,0,1), c(1,0,0,0,0))



combinations.1 = function(size, choose) {
  
  d = do.call("expand.grid", rep(list(0:1), size))
  d[rowSums(d) == choose,]
  
}



# Create q-vectors

########################################


combinations1 = function(size, choose) {
  
  d = do.call("expand.grid", rep(list(0:1), size))
  d[rowSums(d) == choose,]
  
}


combinations2 = function(size, choose) {
  
  d = do.call("expand.grid", rep(list(0:2), size))
  d[rowSums(d) == choose,]
  
}



combinations3 = function(size, choose) {
  
  d = do.call("expand.grid", rep(list(0:3), size))
  d[rowSums(d) == choose,]
  
}


########################################


# Check first order
create.q.1 = function(p){
  
  q1 = combinations1(p, 1)
  
  remove.rows = which( rowSums(q1[,1:p-1]) < 1)
  q1 = q1[-remove.rows,]
  
  colnames(q1) = NA
  
  add.rows = q1[which(q1[, p-1] != 0), ]
  add.rows[, p] = 1
  colnames(add.rows) = NA
  
  # q1 = matrix(as.matrix(rbind(q1, add.rows)), ncol=p)
  q1 = as.matrix(rbind(q1, add.rows), rownames.force = FALSE)
  colnames(q1) = NULL
  
  return(q1)
  
}



# Create Second Order
create.q.2 = function(p){
  
  q12 = combinations2(p, 2)
  
  remove.rows = which( rowSums(q12[,1:p-1]) < 2)
  q12 = q12[-remove.rows,]
  
  colnames(q12) = NA
  
  add.rows = q12[which(q12[, p-1] != 0), ]
  add.rows[, p] = 1
  colnames(add.rows) = NA
  
  # q12 = matrix(as.matrix(rbind(q12, add.rows)), ncol=p)
  q12 = as.matrix(rbind(q12, add.rows), rownames.force = FALSE)
  colnames(q12) = NULL
  
  return(q12)
  
}




# Create Third Order
create.q.3 = function(p){
  
  q12 = combinations3(p, 3)
  
  remove.rows = which( rowSums(q12[,1:p-1]) < 3)
  q12 = q12[-remove.rows,]
  
  colnames(q12) = NA
  
  add.rows = q12[which(q12[, p-1] != 0), ]
  add.rows[, p] = 1
  colnames(add.rows) = NA
  
  # q12 = matrix(as.matrix(rbind(q12, add.rows)), ncol=p)
  q12 = as.matrix(rbind(q12, add.rows), rownames.force = FALSE)
  colnames(q12) = NULL
  
  return(q12)
  
}


####################################


gen.spherical.design1 = function(x, normalize.by.volume = T){
  
  x = as.matrix(x)
  
  p = ncol(x)
  n = nrow(x)
  
  q1 = create.q.1(p)
  
  m1 = nrow(q1)
  temp1 = matrix(0, n, m1)
  
  for(i in 1:n){
    for(j in 1:m1){
      temp1[i,j] = y.q.x(x[i,], q1[j,], normalize.by.volume)
    }
  }
  
  
  out1 = cbind(temp1)
  
  q1.names = apply(q1 , 1 , paste , collapse = "" )
  q1.names = paste("q=", q1.names, sep="")
  
  colnames(out1) = c(q1.names)
  
  return(out1)  
  
}



gen.spherical.design2 = function(x, normalize.by.volume = T){
  
  p = ncol(x)
  n = nrow(x)
  
  # x = as.matrix(x, nrow=x, ncol=p)
  
  q1 = create.q.1(p)
  q2 = create.q.2(p)
  
  m1 = nrow(q1)
  temp1 = matrix(0, n, m1)
  
  for(i in 1:n){
    for(j in 1:m1){
      temp1[i,j] = y.q.x(x[i,], q1[j,], normalize.by.volume)
    }
  }
  
  
  m2 = nrow(q2)
  temp2 = matrix(0, n, m2)
  
  for(i in 1:n){
    for(j in 1:m2){
      temp2[i,j] = y.q.x(x[i,], q2[j,], normalize.by.volume)
    }
  }
  
  
  out1 = cbind(temp1, temp2)
  
  q1.names = apply(q1 , 1 , paste , collapse = "" )
  q1.names = paste("q=", q1.names, sep="")
  
  q2.names = apply(q2 , 1 , paste , collapse = "" )
  q2.names = paste("q=", q2.names, sep="")
  
  colnames(out1) = c(q1.names, q2.names)
  
  return(out1)  
  
}
