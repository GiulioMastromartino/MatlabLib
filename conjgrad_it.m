function [xk,it] = conjgrad_it(A,b,x0,nmax,toll)
%
% [xk,it] = conjgrad_it(A,b,x0,nmax,toll)
%
% Metodo del gradiente coniugato (non precond.) per risolvere Ax = b
%
% Parametri di ingresso:
% A          Matrice quadrata ( n x n ) del sistema lineare Ax = b
% b          Termine noto, vettore colonna
% x0         Iterata iniziale, vettore colonna
% nmax       Numero massimo di iterazioni
% toll       Tolleranza sul criterio d'arresto (residuo normalizzato)
%
% Parametri di uscita:
% xk         Matrice contenente tutte le iterate (vettori approssimanti x)
%            xk = [ x0, x1, x2, ... ] dimensione ( n x ( it + 1 ) )      
% it         Numero di iterazioni effettuate
%
%                                         Politecnico di Milano, 04/04/2024
%

x = x0;
r = b - A*x0;
p = r;
xk = x0;

res_norm = norm(r) / norm(b);
it = 0;
while it < nmax && res_norm > toll
    it = it + 1;
    
    alpha = (p' * r) / (p' * A * p);
    x = x + alpha * p;
    r = r - alpha * A * p;
    beta = (p' * A * r) / (p' * A * p);
    p = r - beta * p;
    
    res_norm = norm(r) / norm(b);
    xk = [xk, x];
end