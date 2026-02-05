function I = simpcomp(a,b,N,f)
%
% I = simpcomp(a,b,N,f)
%
% Formula di Simpson composita
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

x = [ a : h / 2 : b ]; % vettore dei nodi di integrazione

y = f( x );
   
I = ( h / 6 ) * ( y( 1 ) + ...
                    2 * sum( y( 3 : 2 : 2 * N - 1 ) ) + ...
                    4 * sum( y( 2 : 2 : 2 * N ) ) + ...
                    y( 2 * N + 1 ) );


