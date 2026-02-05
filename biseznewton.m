function [xvect,it] = biseznewton(a,b,nmax_b,nmax_n,toll,fun,dfun)
%
% [xvect,it] = biseznewton(a,b,nmax_b,nmax_n,toll,fun,dfun) 
%
% Approssimazione dello zero di f(x) nell'intervallo [a,b],
% applicando prima il metodo di bisezione per l'avvicinamento allo zero
% e successivamente il metodo di Newton
%
% Parametri di ingresso:
% a, b       Estremi intervallo di partenza
% nmax_b     Numero di iterazioni da applicare per il metodo di Bisezione
% nmax_n     Numero massimo di iterazioni per il metodo di Newton
% toll_n     Tolleranza sul criterio d'arresto per il metodo di Newton
% fun, dfun  Funzione e la sua derivata definite come inline
%
% Parametri di uscita:
% xvect      Vettore contenente tutte le iterate calcolate
%            (l'ultima componente e' la soluzione approssimata)
% it         Numero di iterazioni totali (bisezione + Newton) effettuate
%
%                                         Politecnico di Milano, 04/04/2024
%

xvect = [];
it = [];

disp('----------------- mi avvicino allo zero con bisezione ---------------')

% Metodo di Bisezione
toll_b = ( b - a ) / ( 2^( nmax_b + 1) );
[ xvect_b,  it_b ] = bisez( a, b, toll_b, fun );
it = it_b;
xvect = [ xvect_b ];

disp('----------------- lancio il metodo di Newton ---------------')


% Metodo di Newton
xv = xvect( end );
[ xvect_n, it_n ] = newton( xv, nmax_n, toll, fun, dfun );  
it = it + it_n;
xvect = [ xvect; xvect_n ]; 
