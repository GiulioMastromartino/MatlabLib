function [alpha_k,it] = backtracking(Phi,GradPhi,x_k,d_k,c1,rho,nmax)
%
% [alpha_k,it] = backtracking(Phi,GradPhi,x_k,d_k,c1,rho,nmax)
%
% Algoritmo di backtracking per il calcolo di alpha_k (step-size metodi di
% discesa per ottimizzazione). Solo per funzioni obbiettivo in 2-dimensioni 
% (caso n=2) 
%
% Parametri di ingresso:
% Phi       (handle function) Funzione obiettivo
% GradPhi   (handle function) Gradiente della funzione obiettivo
% x_k       (double, vettore) punto corrente
% d_k       (double, vettore) direzione corrente
% c1        (double) parametro metodo, c1 \in (0,1)
% rho       (double) parametro metodo, rho \in [1/10, 1/2]
% nmax      (int) numero massimo di iterazioni
%
% Parametri di uscita:
% alpha_k   (double) valore step-size
% it        (int) numero di iterazioni
%
%                                         Politecnico di Milano, 03/04/2025
%

% Inizializzo alpha_k
alpha_k = 1;

% Inizializzo parametri metodo iterativo
it        = 0;
cond_wolfe = 0;

while cond_wolfe == 0 && it < nmax
    
    % Aggiorno alpha_k (riduco di un fattore rho)
    alpha_k = rho*alpha_k;
    
    % Aggiorno numero di iterazioni
    it = it + 1;
    
    % Aggiorno condizione di uscita
    cond_wolfe = Phi(x_k(1) + alpha_k*d_k(1), x_k(2) + alpha_k*d_k(2)) ...
        <= Phi(x_k(1),x_k(2)) + c1*alpha_k*d_k'*GradPhi(x_k(1), x_k(2));

end

end