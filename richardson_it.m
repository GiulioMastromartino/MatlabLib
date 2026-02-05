function [xk,k] = richardson_it(A,b,P,x0,toll,nmax,alpha)
%
% [xk,k] = richardson_it(A,b,P,x0,toll,nmax,alpha)
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
% toll      Tolleranza criterio d'arresto del residuo normalizzato
% nmax      Numero massimo di iterazioni ammesse
% alpha     Parametro di accelerazione stazionario per l'applicazione del 
%           metodo di Richardson stazionario. Se l'ingresso alpha non e' 
%           assegnato si applica il metodo dinamico, ovvero il metodo del 
%           gradiente precondizionato
%
% Parametri di uscita:
% xk        Matrice contenente tutte le iterate (vettori approssimanti x)
%           xk = [ x0, x1, x2, ... ] dimensione ( n x ( k + 1 ) )      
% k         Numero di iterazioni effettuate
%
%                                         Politecnico di Milano, 04/04/2024
%

n = length(b);
if ((size(A,1) ~= n) || (size(A,2) ~= n) || (length(x0) ~= n))
  error('Dimensioni incompatibili')
end

x = x0;
k = 0;
r    = b - A * x;
res  = norm(r) / norm(b);
xk = [x0];

while ((res > toll) && (k < nmax))
     z = P \ r;
     if (nargin == 6)
         alpha = (z' * r) / (z' * A * z); % alpha dinamico
     end
     x = x + alpha * z;
     r = b - A * x;
     res  = norm(r) / norm(b);
     k = k + 1;
     xk=[xk x];
end

if (res < toll)
     fprintf('Richardson converge in %d iterazioni \n', k);
else
     fprintf('Richardson non converge in %d iterazioni. \n', k)
end

end