function [t_h,u_h,vett_it_newton] = eulero_indietro_newton(f,df_du,t_max,valore_iniziale,delta_t)
%
% [t_h,u_h,vett_it_newton] = eulero_indietro_newton(f,df_du,t_max,valore_iniziale,delta_t)
%
% Approssima il problema di Cauchy scalare utilizzando il 
% metodo di Eulero all'indietro con metodo di Newton
%
% y' = f(t,y)   t \in (0,t_max)
% y(0) = y_0
%
%
% Parametri di ingresso:
% f                 Function che descrive il problema di Cauchy dichiarata
%                   come inline o anonymous @. 
%                   Funzione di due argomenti: f=f(t,y)
% df_du             Derivata di f rispetto a u dichiarata come inline o 
%                   anonymous @. Funzione di due argomenti: 
%                   df_du=df_du(t,y)(è funzione di t e u)
% t_max             Tempo finale dell' intervallo temporale (0,t_max) 
%                   (l'istante iniziale e' t_0=0)
% valore_iniziale   Dato iniziale
% delta_t           Ampiezza del passo di discretizzazione temporale
% 
% Parametri di uscita:
% t_h               Vettore degli istanti in cui si calcola la soluzione 
%                   discreta
% u_h               Vettore della soluzione approssimata agli 
%                   istanti temporali t_h
% vett_it_newton    Vettore che contiene il numero di iterazioni di Newton
%                   effettuate per risolvere l'equazione non lineare ad 
%                   ogni istante temporale (tolleranza prefissata)
%
%                                         Politecnico di Milano, 12/06/2024
%

t0=0;

t_h=t0:delta_t:t_max;

% inizializzo il vettore che conterra' la soluzione discreta

N_istanti=length(t_h);

u_h=zeros(1,N_istanti);

u_h(1)=valore_iniziale;

% parametri per le iterazioni di punto fisso
N_max=100;
toll=1e-5;

vett_it_newton=zeros(1,N_istanti);


for it=2:N_istanti
    
    % preparo le variabili per le sottoiterazioni
    
    u_old=u_h(it-1);
    t_newton=t_h(it);
        
    F=@(u) u_old + delta_t * f( t_newton, u) - u;
    dF = @(u) delta_t * df_du(t_newton, u) - 1;
    
    % sottoiterazioni
    
    [u_newton, it_newton] = newton(u_old, N_max, toll, F, dF);
    
    u_h(it)=u_newton(end);
    
    % tengo traccia dei valori di it_newton per valutare la convergenza 
    % delle iterazioni di punto fisso
    vett_it_newton(it)=it_newton;
    
end

