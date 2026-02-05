function [lambda,x,iter]=invpower(A,tol,nmax,x0)
%
% [lambda,x,iter] = invpower(A,tol,nmax,x0)
% 
% Metodo delle potenze inverse per approssimare l'autovalore di modulo
% minimo di una matrice A (A x = lambda x). Il sistema lineare viene
% risolto applicando il metodo della fattorizzazione LU.
%
% Parametri di ingresso:
% A         Matrice quadrata (n x n)
% tol       Tolleranza sul criterio d'arresto (differenza iterate relativa)
%           (valore default 1e-6) 
% nmax      Numero massimo di iterazioni (valore default 100)
% x0        Iterata iniziale per l'autovettore, vettore colonna 
%           (vettore di default (1,1,...,1)^T)
%
% Parametri di uscita:
% lambda    Approssimazione dell'autovalore di modulo minimo
% x         Approssimazione dell'autovettore corrispondente a lambda (non
%           normalizzato)
% iter      Numero di iterazioni effettuate
% 
%                                         Politecnico di Milano, 04/04/2024
%

[n,m] = size(A);
if n ~= m, error('Solo per matrici quadrate'); end
if nargin == 1
   tol = 1.e-06;   x0 = ones(n,1);   nmax = 100;
end

% calcolo la fattorizzazione LU una volta per tutte
[L,U,P]=lu(A);

% iterazione zero fuori dal ciclo while
iter = 0;
y = x0/norm(x0); % y0
lambda = y'*A*y; % lambda0
err = tol*abs(lambda) + 1; % dobbiamo entrare nel ciclo

while (err>tol*abs(lambda)) && (abs(lambda)~=0) && (iter<nmax)
   iter = iter + 1; % iter=1,2,3,...
   % risolvo Ax^{(k)}=y^{(k-1)}
   z=fwsub(L,P*y);
   x=bksub(U,z);
   y= x/norm(x);
   lambdanew = y'*A*y;
   err = abs(lambdanew - lambda);
   lambda = lambdanew; 
end

if (err <= tol*abs(lambda))
     fprintf('Il metodo delle potenze inverse converge in %d iterazioni \n', iter);
else
     fprintf('Il metodo delle potenze inverse non converge in %d iterazioni. \n', iter)
end

return
