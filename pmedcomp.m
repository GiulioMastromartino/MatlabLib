function I = pmedcomp(a,b,N,f)
%
% I = pmedcomp(a,b,N,f)
%
% Formula del punto medio composita
%
% Parametri di ingresso:
% a,b   Estremi di integrazione, intervallo [a,b]
% N     Numero di sottointervalli equispaziali; se N=1 formula di semplice
% f     Funzione da integrare definita come inline o anonimous
%
% Parametri di uscita:
% I     Valore approssimato dell'integrale tramite la formula
%
%                                         Politecnico di Milano, 12/06/2024
%

h = ( b - a ) / N;

x = [ a + h / 2 : h : b - h / 2 ]; % vettore dei punti medi

I = h * sum( f(x) );


