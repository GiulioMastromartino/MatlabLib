function [t_h,u_h] = Runge_Kutta_4(f,t_max,y_0,h)
%
% [t_h,u_h] = Runge_Kutta_4(f,t_max,y_0,h)
%
% Risolve il problema di Cauchy scalare utilizzando il 
% metodo di Runge-Kutta di ordine 4 esplicito
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

% ciclo iterativo per calcolare la soluzione
for i=2:N_istanti
    k_1 = f (t_h (i-1), u_h (i-1));
    k_2 = f (t_h (i-1) + 0.5 * h, u_h (i-1) + 0.5 * h * k_1);
    k_3 = f (t_h (i-1) + 0.5 * h, u_h (i-1) + 0.5 * h * k_2);
    k_4 = f (t_h (i-1) + h, u_h (i-1) + k_3 * h);

    u_h (i) = u_h (i-1) + 1/6 * (k_1 + 2 * k_2 + 2 * k_3 + k_4) * h;
end
