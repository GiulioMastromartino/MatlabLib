function [x,k,res_norm] = gs(A,b,x0,toll,nmax)
%
% [x,k,res_norm] = gs(A,b,x0,toll,nmax)
%
% Metodo di Gauss-Seidel per l'approssimazione della soluzione di A x = b
%
% Parametri di ingresso:
% A        Matrice del sistema lineare
% b        Termine noto, vettore colonna
% x0       Iterata iniziale, vettore colonna
% toll     Tolleranza sul criterio d'arresto del residuo normalizzato
% nmax     Massimo numero di iterazioni
%
% Parametri di uscita:
% x        Approssimazione del vettore soluzione, vettore colonna
% k        Numero di iterazioni effettuate
% res_norm Residuo normalizzato finale ||b - Ax_k|| / ||b||
%
%                                         Politecnico di Milano, 04/04/2024
%

n = length(b);
b = b(:);
x0 = x0(:);

if (( size(A,1)~=n) || (size(A,2)~=n) || (length(x0) ~= n) )
  error('dimensioni incompatibili')
end

if (prod(diag(A)) == 0)
    error('errore: elementi diagonali nulli')
end

% Matrice triangolare inferiore (parte di M per lo splitting)
M = tril(A); 
% Parte strettamente superiore (N nello splitting A = M - N)
N = -triu(A, 1);

xv = x0;
r = b - A * x0;
norm_b = norm(b);
if norm_b == 0, norm_b = 1; end % Evita div by zero

err = norm(r) / norm_b;
k = 0;

while ( err > toll && k < nmax )
  k = k + 1;
  
  % Risolvo M * x_{k+1} = N * x_k + b
  % x_{k+1} = M \ (N * x_k + b)
  % Poiché M è tril, usiamo fwsub
  rhs = N * xv + b;
  xn = fwsub(M, rhs);
  
  r = b - A * xn;
  err = norm(r) / norm_b;
  xv = xn;
end

x = xv;
res_norm = err;

return
end