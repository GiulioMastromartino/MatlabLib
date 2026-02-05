function [t_h,u_h,vett_it_pf] = eulero_indietro_pto_fisso(f,t_max,valore_iniziale,delta_t)
%
% [t_h,u_h,vett_it_pf] = eulero_indietro(f,t_max,valore_iniziale,delta_t)
%
% Approssima il problema di Cauchy scalare utilizzando il 
% metodo di Eulero all'indietro con metodo delle iterazioni di punto fisso
%
% y' = f(t,y)   t \in (0,t_max)
% y(0) = y_0
%
%
% Parametri di ingresso:
% f                  Function che descrive il problema di Cauchy dichiarata
%                    come inline o anonymous @. 
%                    Funzione di due argomenti: f=f(t,y)
% t_max              Tempo finale dell' intervallo temporale (0,t_max) 
%                    (l'istante iniziale e' t_0=0)
% valore_iniziale    Dato iniziale
% delta_t            Ampiezza del passo di discretizzazione temporale
%
% Parametri di uscita:
% t_h         Vettore degli istanti in cui si calcola la soluzione discreta
% u_h         Vettore della soluzione approssimata agli istanti temporali t_h
% vett_it_pf  Vettore che contiene il numero di iterazioni di punto fisso 
%             effettuate per risolvere l'equazione non lineare ad ogni 
%             istante temporale (tolleranza prefissata)
%
%                                         Politecnico di Milano, 12/06/2024
%


t0=0;

t_h=t0:delta_t:t_max;

% inizializzo il vettore che conterra' la soluzione discreta

N_istanti=length(t_h);

u_h=zeros(1,N_istanti);

% ciclo iterativo che calcola u_(n+1)=u_n+h*f_(n+1) . Ad ogni iterazione temporale devo eseguire delle sottoiterazioni
% di punto fisso per il calcolo di u_(n+1): u_(n+1)^(k+1) = u_n + h * f_( t_(n+1) , u_(n+1)^k ). Per garantire la
% convergenza di tale metodo la derivata della funzione di iterazione phi(x) u_n + h * f_( t_(n+1) , x ) deve avere
% derivata minore di 1 in modulo. Questo introdurra' delle condizioni su h.

u_h(1)=valore_iniziale;

% parametri per le iterazioni di punto fisso
N_max=100;
toll=1e-5;

vett_it_pf=zeros(1,N_istanti);



for it=2:N_istanti
    
    % preparo le variabili per le sottoiterazioni
    
    u_old=u_h(it-1);
    t_pf=t_h(it);
        
    phi=@(u) u_old + delta_t * f( t_pf, u );
    
    % sottoiterazioni
    
    [u_pf, it_pf] = ptofis(u_old, phi, N_max, toll);
    
    u_h(it)=u_pf(end);
    
    % tengo traccia dei valori di it_pf per valutare la convergenza delle iterazioni di punto fisso
    vett_it_pf(it)=it_pf;
    
end

