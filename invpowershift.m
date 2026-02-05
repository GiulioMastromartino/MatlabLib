function [lambda,x,iter] = invpowershift(A,mu,tol,nmax,x0)
% 
% [lambda,x,iter] = invpowershift(A,mu,tol,nmax,x0)
%
% Metodo delle potenze inversecon shift per approssimare l'autovalore di 
% di una matrice A (A x = lambda x) piu' prossimo allo shift mu. Il sistema
% lineare viene risolto applicando il metodo della fattorizzazione LU.
%
% Parametri di ingresso:
% A         Matrice quadrata (n x n)
% mu        Valore di shift
% tol       Tolleranza sul criterio d'arresto (differenza iterate relativa)
%           (valore default 1e-6) 
% nmax      Numero massimo di iterazioni (valore default 100)
% x0        Iterata iniziale per l'autovettore, vettore colonna 
%           (vettore di default (1,1,...,1)^T)
%
% Parametri di uscita:
% lambda    Approssimazione dell'autovalore di A piu' prossimo a mu
% x         Approssimazione dell'autovettore corrispondente a lambda (non
%           normalizzato)
% iter      Numero di iterazioni effettuate
% 
%                                         Politecnico di Milano, 04/04/2024
%

[n,m] = size(A);
if n ~= m, error('Solo per matrici quadrate'); end
if nargin == 2
   tol = 1.e-06;   x0 = ones(n,1);   nmax = 100;
end
M = A - mu*eye(n);
% calcolo la fattorizzazione LU una volta per tutte
[L,U,P]=lu(M);
% iterazione zero fuori dal ciclo while
iter = 0;
y = x0/norm(x0); % y0
lambda = y'*A*y; % lambda0
err = tol*abs(lambda) + 1; % dobbiamo entrare nel ciclo

while (err>tol*abs(lambda)) && (abs(lambda)~=0) && (iter<nmax)
   iter = iter + 1; % iter=1,2,3,...
   % risolvo Mx^{(k)}=y^{(k-1)}
   z=fwsub(L,P*y);
   x=bksub(U,z);
   y= x/norm(x);
   lambdanew = y'*A*y;
   err = abs(lambdanew - lambda);
   lambda = lambdanew; 
end
if (iter < nmax)
     fprintf(['Il metodo delle potenze inverse con shift converge ',...
              'in %d iterazioni all''autovalore \n'], iter);
     lambda
else
     fprintf(['Il metodo delle potenze inverse con shift non converge ',...
              'in %d iterazioni. \n'], iter);
end
return
