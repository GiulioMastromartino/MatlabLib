function [x, k] = richardson(A, b, P, x0, tol, nmax, alpha)
%
% [x, k] = richardson(A, b, P, x0, tol, nmax, alpha)
%
% Metodo di Richardson stazionario precondizionato o metodo
% dinamico precondizionato (ovvero metodo del gradiente precondizionato) 
% per l'approssimazione della soluzione di A x = b. Se P = I, il metodo e'
% non precondizionato.
%
% Parametri di ingresso:
% A         Matrice del sistema lineare (n x n)
% b         Termine noto, vettore colonna (n x 1)
% P         Precondizionatore, matrice di precondizionamento (n x n)
% x0        Iterata iniziale, vettore colonna (n x 1)
% tol       Tolleranza criterio d'arresto del residuo normalizzato
% nmax      Numero massimo di iterazioni ammesse
% alpha     Parametro di accelerazione stazionario per l'applicazione del 
%           metodo di Richardson stazionario. Se l'ingresso alpha non e' 
%           assegnato si applica il metodo dinamico, ovvero il metodo del 
%           gradiente precondizionato
%
% Parametri di uscita:
% x         Approssimazione del vettore soluzione, vettore colonna (n x 1)
% k         Numero di iterazioni effettuate
%
%                                         Politecnico di Milano, 04/04/2024
%

n = length(b);
if ((size(A,1) ~= n) || (size(A,2) ~= n) || (length(x0) ~= n))
  error('Dimensioni incompatibili')
end

% E' possibile utilizzare una sola variabile x al posto di xn e xv viste
% nel laboratorio precedente.

x = x0;
k = 0;
r    = b - A * x;
res_normalizzato  = norm(r) / norm(b);

while ((res_normalizzato > tol) && (k < nmax))
     z = P \ r;
     if (nargin == 6)
         alpha = (z' * r) / (z' * A * z); % alpha dinamico
     end
     x = x + alpha * z;
     r = b - A * x; % equivalente a: r = r - alpha * A * z;
     res_normalizzato  = norm(r) / norm(b);
     k = k + 1;
end

if (res_normalizzato < tol)
     fprintf('Richardson converge in %d iterazioni \n', k);
else
     fprintf('Richardson non converge in %d iterazioni. \n', k)
end
