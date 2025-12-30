###############################################################
# Exact computation of the distribution and moments of R_{v,i}
# ============================================================
# Distribution of S_w
# Based on:
# Gamboa, M. & Lopez-Herrero, M. J. (2022)
# "Measures to assess a warning vaccination level in a
# stochastic SIV model with imperfect vaccine"
# Studies in Applied Mathematics, 148(4), 1411–1438.
# https://doi.org/10.1111/sapm.12479
# ============================================================

# -------------------------------------------------------------
# Parameters: beta, gamma, xi, h, N, w, v0, i0
# Outputs:    Distribution of S_w and moments of R_w
###############################################################


###############################################################
# Auxiliary rate functions
###############################################################

gammai <- function(gamma, i) {
  gamma * i
}

lambdavi <- function(beta, xi, N, v, i) {
  (beta * i / N + xi) * (N - v - i)
}

nuvi <- function(h, beta, xi, N, v, i) {
  h * v * (beta * i / N + xi)
}

qvi <- function(beta, xi, gamma, h, N, v, i) {
  gamma * i +
    (beta * i / N + xi) * (N - v - i) +
    h * v * (beta * i / N + xi)
}


###############################################################
# Model parameters
###############################################################

w<-20
N<-50 
beta<-1.1
gamma<-0.6
xi<-0.3
h<-0.25
v0<-48
i0<-1
negativo<-"F"


###############################################################
# Exact computation of the probability mass function of I_{v,i}
#
# The recursive scheme below is stable for moderate values of N.
# For large N, the computation of S_w is carried out using an
# alternative script specifically designed for large populations.
# and it is also included in
# the repository
###############################################################

xvmenos1ik<-matrix(0,nrow = N-w+1,ncol = N-w) #la i va desde 0 y la k desde 1
for (k in 1:(N-w)) {
  for (i in 0:(N-w)) {
    if (i==k){xvmenos1ik[i+1,k]<-1}
  }
}

for (v in (w+1):v0) { 
  
  cvi<-1
  cvi<-c(cvi,qvi(beta,xi,gamma,h,N,v,0))
  for (i in 1:(N-v)) {
    a<-qvi(beta,xi,gamma,h,N,v,i)*cvi[i+1]-gammai(gamma,i)*lambdavi(beta,xi,N,v,i-1)*cvi[i]
    cvi<-c(cvi,a)
  }
  dvik<-matrix(nrow = N-v+1,ncol = N-w) 
  xvik<-matrix(nrow = N-v+1,ncol = N-w)  
  for (k in 1:(N-w)) {
    dvik[1,k]<-nuvi(h,beta,xi,N,v,0)*xvmenos1ik[2,k]
    for (i in 1:(N-v)) {
      dvik[i+1,k]<-gammai(gamma,i)*dvik[i,k]+nuvi(h,beta,xi,N,v,i)*cvi[i+1]*xvmenos1ik[i+2,k]
    }
    xvik[N-v+1,k]<-dvik[N-v+1,k]/cvi[length(cvi)]
    for (i in (N-v-1):0){
      suma<-0
      for (j in i:(N-v-1)) {
        if (j>(N-v-1)){mm<-0}else{
          mm<-dvik[j+1,k]/cvi[j+2]}
        prod<-1
        for (m in i:(j-1)) {
          if (i>(j-1)){prod<-1}else{
            prod<-prod*lambdavi(beta,xi,N,v,m)*cvi[m+1]/cvi[m+2]
           }
        }
        suma<-suma+mm*prod
         }
      producto<-1
      for (m in i:(N-v-1)) {
        xx<-lambdavi(beta,xi,N,v,m)*cvi[m+1]/cvi[m+2]
        producto<-producto*xx
      }
      producto<-producto*xvik[N-v+1,k]
      xvik[i+1,k]<-suma+producto 
    }
  }
  xvmenos1ik<-xvik
}

fprobav0io1<-xvik[2,]
sum(fprobav0io1)
show(fprobav0io1)


#Momentos de las variables R_vi
##############################################

k<-0
v<-0
i<-0
mvi0<-matrix(nrow = N-v+1,ncol = w+1)
for (v in 0:w) {
  for (i in 0:(N-v)) {
    mvi0[i+1,v+1]<-1
  }
}
mvi1<-matrix(nrow = N+1,ncol = w+1)  
for (v in 0:w) {
  for (i in 0:(N-v0)) {
    mvi1[i+1,v+1]<-0
  }
}
Avi1<-matrix(nrow = N+1,ncol = w+1) 
Bvi1<-matrix(nrow = N+1,ncol = w+1) 
Gvi1<-matrix(nrow = N+1,ncol = w+1) 
for (v in 0:w) {
  Avi1[N-v+1,v+1]<-gammai(gamma,N-v)
  Bvi1[N-v+1,v+1]<-qvi(beta,xi,gamma,h,N,v,N-v)
}
Gvi1[,1]<-1
v<-0
for (i in (N-v-1):(N-v0+1)) {
  
  Avi1[i+1,v+1]<-gammai(gamma,i)* Bvi1[i+2,v+1]
  Bvi1[i+1,v+1]<-qvi(beta,xi,gamma,h,N,v,i)* Bvi1[i+2,v+1]-lambdavi(beta,xi,N,v,i)*Avi1[i+2,v+1]
  Gvi1[i+1,v+1]<-Bvi1[i+2,v+1]+lambdavi(beta,xi,N,v,i)*Gvi1[i+2,v+1]
}
mvi1[N-v0+2,v+1]<-Gvi1[N-v0+2,v+1]/Bvi1[N-v0+2,v+1]
for (i in (N-v0+2):(N-v)) {
  suma<-0
  for (l in (N-v0+2):i) {
    m<-Gvi1[l+1,v+1]/Bvi1[l+1,v+1]
    moment<-1
    if ((l+1)>i){moment<-1}else{
      for (j in (l+1):i) {
        moment<-moment*(Avi1[j+1,v+1]/Bvi1[j+1,v+1])
      }
    }
    moment<-moment*m
    suma<-suma+moment
  }
  producto<-1
  for (l in (N-v0+2):i) {
    producto<-producto*(Avi1[l+1,v+1]/Bvi1[l+1,v+1]) 
  }
  producto<-producto*mvi1[N-v0+2,v+1]
  mvi1[i+1,v+1]<-suma+producto
}

for (v in 1:w) {
  Gvi1[N-v+1,v+1]<-nuvi(h,beta,xi,N,v,N-v)*mvi1[N-v+2,v]+1
  for (i in (N-v-1):(N-v0+1)) {
    Avi1[i+1,v+1]<-gammai(gamma,i)* Bvi1[i+2,v+1]
    Bvi1[i+1,v+1]<-qvi(beta,xi,gamma,h,N,v,i)* Bvi1[i+2,v+1]-lambdavi(beta,xi,N,v,i)*Avi1[i+2,v+1]
    Gvi1[i+1,v+1]<-(nuvi(h,beta,xi,N,v,i)*mvi1[i+2,v]+1)*Bvi1[i+2,v+1]+lambdavi(beta,xi,N,v,i)*Gvi1[i+2,v+1]
  }
  mvi1[N-v0+2,v+1]<-Gvi1[N-v0+2,v+1]/Bvi1[N-v0+2,v+1]
  for (i in (N-v0+2):(N-v)) {
    suma<-0
    for (l in (N-v0+2):i) {
      m<-Gvi1[l+1,v+1]/Bvi1[l+1,v+1]
      moment<-1
      if ((l+1)>i){moment<-1}else{
        for (j in (l+1):i) {
          moment<-moment*(Avi1[j+1,v+1]/Bvi1[j+1,v+1])
        }
      }
      moment<-moment*m
      suma<-suma+moment
    }
    producto<-1
    for (l in (N-v0+2):i) {
      producto<-producto*(Avi1[l+1,v+1]/Bvi1[l+1,v+1]) 
    }
    producto<-producto*mvi1[N-v0+2,v+1]
    mvi1[i+1,v+1]<-suma+producto
  }
}

mvi2<-matrix(nrow = N+1,ncol = w+1)  
for (v in 0:w) {
  for (i in 0:(N-v0)) {
    mvi2[i+1,v+1]<-0
  }
}
Avi2<-matrix(nrow = N+1,ncol = w+1) 
Bvi2<-matrix(nrow = N+1,ncol = w+1) 
Gvi2<-matrix(nrow = N+1,ncol = w+1) 
for (v in 0:w) {
  Avi2[N-v+1,v+1]<-gammai(gamma,N-v)
  Bvi2[N-v+1,v+1]<-qvi(beta,xi,gamma,h,N,v,N-v)
}
Gvi2[N+1,1]<-2*mvi1[N+1,1]
v<-0
for (i in (N-v-1):(N-v0+1)) {
  
  Avi2[i+1,v+1]<-gammai(gamma,i)* Bvi2[i+2,v+1]
  Bvi2[i+1,v+1]<-qvi(beta,xi,gamma,h,N,v,i)* Bvi2[i+2,v+1]-lambdavi(beta,xi,N,v,i)*Avi2[i+2,v+1]
  Gvi2[i+1,v+1]<-2*mvi1[i+1,v+1]*Bvi2[i+2,v+1]+lambdavi(beta,xi,N,v,i)*Gvi2[i+2,v+1]
}
mvi2[N-v0+2,v+1]<-Gvi2[N-v0+2,v+1]/Bvi2[N-v0+2,v+1]
for (i in (N-v0+2):(N-v)) {
  suma<-0
  for (l in (N-v0+2):i) {
    m<-Gvi2[l+1,v+1]/Bvi2[l+1,v+1]
    moment<-1
    if ((l+1)>i){moment<-1}else{
      for (j in (l+1):i) {
        moment<-moment*(Avi2[j+1,v+1]/Bvi2[j+1,v+1])
      }
    }
    moment<-moment*m
    suma<-suma+moment
  }
  producto<-1
  for (l in (N-v0+2):i) {
    producto<-producto*(Avi2[l+1,v+1]/Bvi2[l+1,v+1]) 
  }
  producto<-producto*mvi2[N-v0+2,v+1]
  mvi2[i+1,v+1]<-suma+producto
}

for (v in 1:w) {
  Gvi2[N-v+1,v+1]<-nuvi(h,beta,xi,N,v,N-v)*mvi2[N-v+2,v]+2*mvi1[N-v+1,v+1]
  for (i in (N-v-1):(N-v0+1)) {
    Avi2[i+1,v+1]<-gammai(gamma,i)* Bvi2[i+2,v+1]
    Bvi2[i+1,v+1]<-qvi(beta,xi,gamma,h,N,v,i)* Bvi2[i+2,v+1]-lambdavi(beta,xi,N,v,i)*Avi2[i+2,v+1]
    Gvi2[i+1,v+1]<-(nuvi(h,beta,xi,N,v,i)*mvi2[i+2,v]+2*mvi1[i+1,v+1])*Bvi2[i+2,v+1]+lambdavi(beta,xi,N,v,i)*Gvi2[i+2,v+1]
  }
  mvi2[N-v0+2,v+1]<-Gvi2[N-v0+2,v+1]/Bvi2[N-v0+2,v+1]
  for (i in (N-v0+2):(N-v)) {
    suma<-0
    for (l in (N-v0+2):i) {
      m<-Gvi2[l+1,v+1]/Bvi2[l+1,v+1]
      moment<-1
      if ((l+1)>i){moment<-1}else{
        for (j in (l+1):i) {
          moment<-moment*(Avi2[j+1,v+1]/Bvi2[j+1,v+1])
        }
      }
      moment<-moment*m
      suma<-suma+moment
    }
    producto<-1
    for (l in (N-v0+2):i) {
      producto<-producto*(Avi2[l+1,v+1]/Bvi2[l+1,v+1]) 
    }
    producto<-producto*mvi2[N-v0+2,v+1]
    mvi2[i+1,v+1]<-suma+producto
  }
} 

######################################################################

mwiesperanza<-mvi1[(N-v0+2):(N-w+1),w+1]
mwimm2<-mvi2[(N-v0+2):(N-w+1),w+1]



MRW1<-fprobav0io1[(N-v0+1):(N-w)]*mwiesperanza
MRW1<- sum(MRW1)
cat("El M_{Rw}^1=",MRW1)

MRW2<-fprobav0io1[(N-v0+1):(N-w)]*mwimm2
MRW2<- sum(MRW2)
cat("El M_{Rw}^2=",MRW2)