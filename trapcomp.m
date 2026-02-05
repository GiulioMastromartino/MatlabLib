function I = trapcomp(a,b,N,f)
%
% I = trapcomp(a,b,N,f)
%
% Formula dei trapezi composita
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

x = [ a : h : b ]; % vettore dei nodi di integrazione

y = f( x );
   
I = h * ( 0.5 * y( 1 ) + sum( y( 2 : N ) ) + 0.5 * y( N + 1 ) );


