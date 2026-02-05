function [x,it]=jacobi(A,b,x0,toll,nmax)
%
% [x,it] = jacobi(A,b,x0,toll,nmax)
%
% Metodo di Jacobi per l'approssimazione della soluzione di A x = b
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
% it       Numero di iterazioni effettuate
%
%                                         Politecnico di Milano, 04/04/2024
%

n = length(b);
it = 0;

if ((size(A,1) ~= n) || (size(A,2) ~= n) || (length(x0) ~= n))
  error('Dimensioni incompatibili')
end

if (prod(diag(A)) == 0)
  error('Errore: elementi diagonali nulli')
end

D_inv = diag(1./diag(A));
x = x0;
r = b - A*x;
err = norm(r) / norm(b);

while (err > toll && it < nmax)
    it = it + 1;
    z = D_inv*r;
    x = x + z;
    r = b - A*x;
    err = norm(r)/norm(b);
end



