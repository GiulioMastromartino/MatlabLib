function [t,u] = eulero_avanti_sistemi(fun,tv,y0,Nh)
%
% [t,u] = eulero_avanti_sistemi(fun,tv,y0,Nh)
%
% Approssima il problema di Cauchy vettoriale utilizzando il 
% metodo di Eulero in avanti
%
% y' = f(t,y)   t \in (t_0,t_max)
% y(t_0) = y0
% 
% m >= 1,   y : R -> R^m,     f : R x R^m -> R^m,    y0 \in R^m, 
% Parametri di ingresso:
%
% fun       Function che descrive il problema di Cauchy dichiarata come 
%           inline o anonymous @. Funzione di due argomenti: f=fun(t,y)
% tv        Vettore riga per l'intervallo temporale (t_0,t_max) 
% y0        Dato iniziale, vettore colonna m x 1
% Nh        Numero di passi temporali di ampiezza h
%
% Parametri di uscita:
% t         Vettore degli istanti in cui si calcola la soluzione discreta
% u         Matrice della soluzione approssimata agli istanti temporali t
%           dimensione: m x (Nh+1)
%
%                                         Politecnico di Milano, 12/06/2024
%

h = ( tv( end ) - tv( 1 ) ) / Nh;

u = zeros( size( y0, 1 ), Nh + 1 );
t = linspace( tv( 1 ), tv( end ), Nh + 1 );

u( :, 1 ) = y0;

for n = 1 : Nh
    u( :, n + 1 ) = u( :, n ) + h * fun( t( n ), u( :, n ) );
        
end

return