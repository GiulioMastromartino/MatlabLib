function [ L,U ] = lugauss( A )
%
% [L,U] = lugauss( A )
% 
% Metodo di Eliminazione di Gauss senza pivoting per la costruzione dei
% fattori L (matrice triangolare inferiore) ed U (matrice triangolare
% superiore) di una matrice A quadrata. A = L U
%
% Parametri di ingresso:
% A     Matrice quadrata da fattorizzare (n x n)
%
% Parametri di uscita:
% L     Matrice triangolare inferiore (n x n), con elementi pari ad 1 sulla
%       diagonale principale
% U     Matrice triangolare superiore (n x n) 
%
%                                         Politecnico di Milano, 04/04/2024
%

[n, m] = size (A);

if (n ~= m)
    error ('A non e'' una matrice quadrata'); 
end

L = eye (n);
for k = 1:n-1
    if (A (k,k) == 0)
        error ('Un elemento pivot si e'' annullato'); 
    end
    for i = k+1:n % fissa l'indice di riga               
        L (i, k) = A (i, k) / A (k, k);       
        for j = k+1:n % a partire dall'elemento fissato scorre la riga
            A (i, j) = A (i, j) - L (i, k) * A (k, j);            
        end
    end
end

U = triu(A);














