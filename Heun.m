function [t_h,u_h] = Heun(f,t_max,y_0,h)
%
% [t_h,u_h] = Heun(f,t_max,y_0,h)
%
% Approssima il problema di Cauchy (scalare o vettoriale) utilizzando il 
% metodo di Heun.
%
% y' = f(t,y)   t \in (0,t_max)
% y(0) = y_0
%
% Parametri di ingresso:
% f         Function che descrive il problema di Cauchy dichiarata come 
%           inline o anonymous @. Funzione di due argomenti: f=f(t,y).
%           Deve restituire un vettore colonna se y è un vettore.
% t_max     Tempo finale dell' intervallo temporale (0,t_max) 
%           (l'istante iniziale e' t_0=0)
% y_0       Dato iniziale (scalare o vettore colonna)
% h         Ampiezza del passo di discretizzazione temporale (h>0)
%
% Parametri di uscita:
% t_h       Vettore degli istanti in cui si calcola la soluzione discreta
% u_h       Matrice della soluzione approssimata agli istanti temporali t_h
%           (dimensione neq x N_istanti)
%
%                                         Politecnico di Milano, 12/06/2024
%

% Assicuriamoci che y_0 sia una colonna
y_0 = y_0(:);
neq = length(y_0);

% vettore degli istanti in cui risolvo la edo
t0 = 0;
t_h = t0:h:t_max;

% inizializzo il vettore/matrice che conterra' la soluzione discreta
N_istanti = length(t_h);
u_h = zeros(neq, N_istanti);

u_h(:,1) = y_0;

for it = 2 : N_istanti
    % t_n corrente (step precedente)
    tn = t_h(it-1);
    un = u_h(:,it-1);
    
    % Predittore (Eulero esplicito)
    k1 = f(tn, un);
    
    % Correttore (Trapezi esplicito / Heun)
    u_pred = un + h * k1;
    t_next = t_h(it);
    k2 = f(t_next, u_pred);
    
    u_h(:,it) = un + (h/2) * (k1 + k2);
end

end