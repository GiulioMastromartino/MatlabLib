function [lambda,x,iter] = invpowershift(A,mu,tol,nmax,x0)
% 
% [lambda,x,iter] = invpowershift(A,mu,tol,nmax,x0)
%
% Metodo delle potenze inverse con shift per approssimare l'autovalore di 
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

% Gestione parametri opzionali
if nargin < 5
    x0 = ones(n,1);
end
if nargin < 4
    nmax = 100;
end
if nargin < 3
    tol = 1.e-06;
end

% Matrice shiftata M = A - mu*I
M = A - mu*eye(n);

% Fattorizzazione LU una volta per tutte
[L,U,P] = lu(M);

% Inizializzazione
iter = 0;
x0 = x0(:); % Assicura colonna
if norm(x0) == 0, error('Vettore iniziale nullo'); end

y = x0 / norm(x0); 
% Stima iniziale dell'autovalore (quoziente di Rayleigh)
lambda = y' * A * y; 

% Entra nel ciclo almeno una volta
err = tol * abs(lambda) + 1; 

while (err > tol * abs(lambda)) && (iter < nmax)
   iter = iter + 1;
   
   % Risolvo (A - mu*I) * x_k = y_{k-1}
   % Usando la fattorizzazione: L*U*x = P*y
   z = fwsub(L, P*y);
   x_k = bksub(U, z);
   
   % Normalizzo
   y = x_k / norm(x_k);
   
   % Quoziente di Rayleigh per lambda
   lambdanew = y' * A * y;
   
   % Errore relativo sull'autovalore
   if abs(lambda) > eps
       err = abs(lambdanew - lambda) / abs(lambda);
   else
       err = abs(lambdanew - lambda);
   end
   
   lambda = lambdanew; 
end

% Assegna autovettore finale
x = y;

if (iter < nmax)
     fprintf(['Il metodo delle potenze inverse con shift converge ',...
              'in %d iterazioni all''autovalore \\n'], iter);
else
     fprintf(['Il metodo delle potenze inverse con shift non converge ',...
              'in %d iterazioni. \\n'], iter);
end

return
end