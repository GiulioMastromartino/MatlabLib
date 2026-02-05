function [succ, it] = ptofis(x0, phi, nmax, toll)
%
% [succ, it] = ptofis(x0, phi, nmax, toll)
%
% Metodo delle iterazioni di punto fisso x_{k+1} = phi( x_{k} ), k>=0 
%
% Parametri di ingresso:
% x0      Iterata iniziale
% phi     Funzione di punto fisso (definita inline o anonimous)
% nmax    Numero massimo di iterazioni
% toll    Tolleranza sul criterio d'arresto (differenza iterate successive)
%
% Parametri di uscita:
% succ    Vettore contenente tutte le iterate calcolate (l'ultima 
%         componente e' l'approssimazione del punto fisso di phi)
% it      Numero di iterazioni effettuate
%
%                                         Politecnico di Milano, 04/04/2024
%

err   = 1 + toll;
it    = 0;
succ  = x0;
xv    = x0;
while (it < nmax && err > toll)
   xn    = phi(xv);
   err   = abs(xn - xv);
   succ = [succ; xn];
   it    = it + 1;
   xv    = xn;
end
%fprintf(' \n Numero di Iterazioni    : %d \n',it);
%fprintf(' Punto fisso calcolato   : %12.13f \n',succ(end));
