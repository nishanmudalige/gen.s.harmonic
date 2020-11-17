# Created by Nishan Mudalige and Dr. Peter Kim
# Maintained by Nishan Mudalige. Email: nishan.mudalige.1@ulaval.ca
# Copyright Nishan Mudalige and Dr. Peter Kim, 2020

# This research us partially based upon work supported by a research grant 
# from Google Cloud.

# We use the general spherical harmonics as defined by Dai and Xu
# in Approximation Theory and Harmonic Analysis on Spheres and Balls

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

# function to duplicate a row times
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

# function to duplicate a column times
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

# The b function as defined by Theorem1.5.1. of Dai and Xu
b.q = function(q){
  p = length(q)
  ifelse(q[p-1] + q[p] > 0, 2, 1)
}


# The g function as defined by Theorem1.5.1. of Dai and Xu
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


# The |q| function as defined by Theorem1.5.1. of Dai and Xu
q.norm = function(q, j){
  p_1 = length(q)-1
  return(sum(q[j:p_1]))
}


# The lambda_j function as defined by Theorem1.5.1. of Dai and Xu
lambda = function(q, j){
  p = length(q)
  q.norm(q, j+1) + (p-j-1)/2
}


# The Gegenbauer polynomial C^{lambda}_{n}(x)
gegenbauer = function(n, lambda, x){
  
  res = ( (orthopolynom::pochhammer(2*lambda, n))/factorial(n) ) * 
    ( hypergeo::hypergeo(tol=0.1, A=-n, B=n+(2*lambda), C=lambda+0.5, z=(1-x)/2) )
  
  return(Re(res))
}


# The h_{\alpha}^{2} function as defined by Theorem1.5.1. of Dai and Xu
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


# The volume of the hyper-sphere
# By default, the unit sphere is considered 
volume.hyper.sphere = function(p, r=1){
  volume = 2*pi^(p/2) / gamma(p/2) * (r^p)
  return(volume)
}


# The spherical harmonic Y_{\alpha}(x) as defined by Theorem1.5.1. of Dai and Xu
# Note the minor error in Theorem1.5.1. of Dai and Xu where 
# (h_{\alpha}^{2})^{-1} should be replaced with h_{\alpha}
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