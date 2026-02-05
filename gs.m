function [x,k] = gs(A,b,x0,toll,nmax)
%
% [x,k] = gs(A,b,x0,toll,nmax)
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
%
%                                         Politecnico di Milano, 04/04/2024
%

n = length(b);
xn = zeros( n, 1 );
k = 0;

if (( size(A,1)~=n) || (size(A,2)~=n) || (length(x0) ~= n) )
  error('dimensioni incompatibili')
end

if (prod(diag(A)) == 0)
    error('errore: elementi diagonali nulli')
end

T = tril(A);
xv = x0;
r = b - A * x0;
err = norm(r) / norm(b);

while ( err > toll && k < nmax )
  k = k + 1;
  z = fwsub(T,r);
  xn = xv + z;
  r = b - A*xn;  
  err = norm(r) / norm(b);
  xv = xn;
end

x = xn;

