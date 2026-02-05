function [t_h,u_h] = Heun(f,t_max,y_0,h)
%
% [t_h,u_h] = Heun(f,t_max,y_0,h)
%
% Approssima il problema di Cauchy scalare utilizzando il 
% metodo di Heun
%
% y' = f(t,y)   t \in (0,t_max)
% y(0) = y_0
%
% Parametri di ingresso:
% f         Function che descrive il problema di Cauchy dichiarata come 
%           inline o anonymous @. Funzione di due argomenti: f=f(t,y)
% t_max     Tempo finale dell' intervallo temporale (0,t_max) 
%           (l'istante iniziale e' t_0=0)
% y_0       Dato iniziale
% h         Ampiezza del passo di discretizzazione temporale (h>0)
%
% Parametri di uscita:
% t_h       Vettore degli istanti in cui si calcola la soluzione discreta
% u_h       Vettore della soluzione approssimata agli istanti temporali t_h
%
%                                         Politecnico di Milano, 12/06/2024
%

% vettore degli istanti in cui risolvo la edo

t0=0;

t_h=t0:h:t_max;

% inizializzo il vettore che conterra' la soluzione discreta

N_istanti=length(t_h);

u_h=zeros(1,N_istanti);

u_h(1)=y_0;

for it=2 : N_istanti
    
    u_h(it) = u_h(it-1) + h/2 * ( f( t_h(it-1), u_h(it-1) ) + ...
                                  f( t_h(it), u_h(it-1) + h * f( t_h(it-1), u_h(it-1) ) ) ); 
end

