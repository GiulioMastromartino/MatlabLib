function [mu,x,min_eigenvec_comp] = invpower_mod(A,nmax,x0)
%
% [mu,x,min_eigenvec_comp] = invpower_mod(A,nmax,x0)
% 
% Metodo delle potenze inverse per approssimare l'autovalore di modulo
% minimo di una matrice A (A x = lambda x). Il sistema lineare viene
% risolto applicando il metodo della fattorizzazione LU.
%
% Parametri di ingresso:
% A                  Matrice quadrata (n x n) 
% nmax               Numero massimo di iterazioni
% x0                 Iterata iniziale per l'autovettore, vettore colonna 
%
% Parametri di uscita:
% mu                 Approssimazione dell'autovalore di modulo minimo dopo 
%                    nmax iterazioni
% x                  Approssimazione dell'autovettore corrispondente a mu 
%                    (non normalizzato)
% min_eigenvec_comp  Vettore di lunghezza NMAX+1. Detta y_i 
%                    l'approssimazione dell'autovettore ottenuta 
%                    all'iterazione i-esima min_eigenvec_comp(i) contiene 
%                    la componente di y_i nella direzione definita 
%                    dall'autovettore associato all'autovalore di modulo 
%                    minimo della matrice A.
% 
%                                         Politecnico di Milano, 04/04/2024
%

[V,D] = eig(A);
D = diag(D);
[~,ind] = min(abs(D));
eigvec = V(:,ind);
norm_eigvec2 = norm(eigvec)^2;

[n,m] = size(A);
if n ~= m, error('Solo per matrici quadrate'); end

% calcolo la fattorizzazione LU una volta per tutte
[L,U,P]=lu(A);

% iterazione zero fuori dal ciclo while
iter = 0;
y = x0/norm(x0); % y0
min_eigenvec_comp = abs(y'*eigvec/norm_eigvec2);
mu = y'*A*y; % mu0

while abs(mu) ~= 0 && iter<nmax
   iter = iter + 1;
   % risolvo Ax^{(k)}=y^{(k-1)}
   z=fwsub(L,P*y);
   x=bksub(U,z);
   y= x/norm(x);
   min_eigenvec_comp = [min_eigenvec_comp, abs(y'*eigvec/norm_eigvec2)];
   munew = y'*A*y;
   err = abs(munew - mu);
   mu = munew; 
end


return
