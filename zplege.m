function [x,w] = zplege(n,a,b)
%
% [x,w] = zplege(n,a,b)
% 
% Calcola i nodi e i pesi di quadratura della formula di quadratura di
% Gauss-Legendre di ordine n riferiti all'intervallo [a,b]
%
% Parametri di ingresso:
% n     Ordine della formula (n+1 = numero dei nodi di quadratura)
% a,b   Estremi dell'intervallo [a,b]
%
% Parametri di uscita:
% x     Vettore (colonna) contenente le coordinate dei nodi di quadratura
% w     Vettore (colonna) contenente i corrispondenti pesi della formula
%
%                                         Politecnico di Milano, 12/06/2024
%


n = n + 1;
if n <= 1
   z = 0;   p = 2;
   x = (a+b)/2 + (b-a)*z/2;
   w = (b-a)*p/2;
   return
end
jac = zeros(n);
k  = [1:n-1];
v  = k./(sqrt(4*k.^2-1));
jac = jac+diag(v,1)+diag(v,-1);
[p,z] = eig(jac);
norm2 = sqrt(diag(p'*p));    
p = (2*p(1,:)'.^2)./norm2;   
z = diag(z);

x = (a+b)/2 + (b-a)*z/2;
w = (b-a)*p/2;