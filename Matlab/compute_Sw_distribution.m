clear; clc;

%% ============================================================
% Distribution of S_w
% Based on:
% Gamboa, M. & Lopez-Herrero, M. J. (2022)
% "Measures to assess a warning vaccination level in a
% stochastic SIV model with imperfect vaccine"
% Studies in Applied Mathematics, 148(4), 1411–1438.
% https://doi.org/10.1111/sapm.12479
%% ============================================================

%% Parameters
N  = 500;      % Total population
w  = 200;      % Warning vaccination level
v0 = 480;      % Final vaccination level
i0=1;
beta = 6.5;
xi   = 0.01;
gamma = 1;
h    = 0.1;


%% Example of use:

[fprobSw]=fprobabilidadSW(w,N,beta,gamma,xi,h,v0,i0);

%% --------------------------
function[fprobSw]=fprobabilidadSW(w,N,beta,gamma,xi,h,v0,i0)
% ------------------------------------------------------------
% Computes distribution of Sw
%
% Inputs:
%   w     : warning level
%   N     : total population size
%   beta  : internal transmission rate
%   gamma : recovery rate
%   xi    : external transmission rate
%   h     : vaccine failure probability
%   v0    : initial number of vaccinated individuals
%   i0    : initial number of infected individuals
%
% ------------------------------------------------------------
fprobv01=fprobabilidadI0i0Thomas(w,N,beta,gamma,xi,h,v0,i0);
for k = 0:(N-w-1) 
  fprobSw(k+1)=fprobv01(N-w-k);
end
end
function[x]=sistematridiag(a,b,c,d)
     n=length(d);
     w(1)=b(1);
     g(1)=d(1)/w(1);
     for i=2:n
          w(i)=b(i)-a(i)*c(i-1)/w(i-1);
          g(i)=(d(i)-a(i)*g(i-1))/w(i);
     end
     x(n)=g(n);
     for i=n-1:-1:1
         x(i)=g(i)-c(i)*x(i+1)/w(i);
     end
 end    
function[x]=qvi(beta,xi,gamma,h,N,v,i)
  x=gamma*i+(beta*i/N+xi)*(N-v-i)+h*v*(beta*i/N+xi);
end 
function[x]=lambdavi(beta,xi,N,v,i)
  x=(beta*i/N+xi)*(N-v-i);
end
function[x] =gammai(gamma,i)
  x=gamma*i;
end
function[x]=nuvi(h,beta,xi,N,v,i)
  x=h*v*(beta*i/N+xi);
end
function[dnu]=diagnuv(h,beta,xi,N,v)
 dnu=0:(N-v);
 for i=  0:(N-v) 
    dnu(i+1)=nuvi(h,beta,xi,N,v,i);
 end
  dnu=diag(dnu);
end
function[tridiagonal]=tridiag(a,b,c) 
tridiagonal = diag(a, 0) + diag(b, -1) + diag(c, 1);
 end
function[QVV]=matrizQvv(beta,xi,gamma,h,N,v)
 if xi== 0
     prompt = 'NO es un modelo SIVS con reinfeccion externa';
 else
     for k = 1:(N-v)
        b(k)=-k*gamma;
     end
     for k = 0:(N-v)
      a(k+1)=qvi(beta,xi,gamma,h,N,v,k);
     end
    for k = 0:(N-v-1)
      c(k+1)=-lambdavi(beta,xi,N,v,k);
    end
    QVV=tridiag(a, b, c);
 end
end   
function[fprob01]=fprobabilidadI0i0(w,N,beta,gamma,xi,h,v0,i0)
matrizxvmenos1k=eye(N-w);
 for v =(w+1):v0 
    Dv=diagnuv(h,beta,xi,N,v); 
    Qv=matrizQvv(beta,xi,gamma,h,N,v);
     for k = 1:(N-w) 
     xvmenos1k=matrizxvmenos1k(:,k);
     b=Dv * xvmenos1k;
     xvk=Qv\b;
      if xvk<0 
         negativo='T';
      end
     matrizxvk(1:length(xvk),k)=xvk;
    end
  matrizxvmenos1k=matrizxvk; 
  matrizxvmenos1k=matrizxvmenos1k(2:(N-v+1),:);
 end
fprob01=matrizxvk(i0+1,:);
end 
function[a,b,c]=matrizQvvVectores(beta,xi,gamma,h,N,v)
 if xi== 0
     prompt = 'NO es un modelo SIVS con reinfeccion externa';
 else
     a=zeros(N-v+1,1);
     for k = 0:(N-v)
        a(k+1)=-k*gamma;
     end
      b=zeros(N-v+1,1);
     for k = 0:(N-v)
      b(k+1)=qvi(beta,xi,gamma,h,N,v,k);
     end
      c=zeros(N-v+1,1);
    for k = 0:(N-v)
      c(k+1)=-lambdavi(beta,xi,N,v,k);
    end
   
 end
end 
function[fprob01]=fprobabilidadI0i0Thomas(w,N,beta,gamma,xi,h,v0,i0)
matrizxvmenos1k=eye(N-w);
 for v =(w+1):v0 
    Dv=diagnuv(h,beta,xi,N,v); 
    [a,b,c]=matrizQvvVectores(beta,xi,gamma,h,N,v);
    for k = 1:(N-w) 
     xvmenos1k=matrizxvmenos1k(:,k);
     d=Dv * xvmenos1k;
     xvk=sistematridiag(a,b,c,d);
      if xvk<0 
         negativo='T';
      end
     matrizxvk(1:length(xvk),k)=xvk;
    end
  matrizxvmenos1k=matrizxvk; 
  matrizxvmenos1k=matrizxvmenos1k(2:(N-v+1),:);
 end
fprob01=matrizxvk(i0+1,:);
end 

