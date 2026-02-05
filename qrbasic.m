function D = qrbasic(A,tol,nmax)
%
% D = qrbasic(A,tol,nmax)
%
% Metodo delle iterazioni QR per l'approssimazione degli n autovalori di 
% una matrice A quadrata a valori reali.
%
% Parametri di ingresso:
% A       Matrice a valori reali quadrata (n x n)
% tol     Tolleranza sul criterio d'arresto
% nmax    Numero massimo di iterazioni
%
% Parametri di uscita:
% D       Vettore colonna (n x 1) contenente gli n autovalori approssimati
%
%                                         Politecnico di Milano, 04/04/2024
%

[n,m]=size(A);
if n ~= m
  error('La matrice deve essere quadrata')
end
T = A; 
niter = 0; 
test = max(max(abs(tril(T,-1))));
while niter < nmax && test > tol
  [Q,R]=qr(T);    
  T = R*Q;
  niter = niter + 1;
  test = max(max(abs(tril(T,-1))));
end
if niter > nmax
 fprintf(['Il metodo non converge nel massimo',...
             ' numero di iterazioni permesso']);
else
 fprintf('Il metodo converge in %d iterazioni\n',niter)
end
D = diag(T);
return
