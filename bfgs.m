function [xvect,it] = bfgs(Phi,GradPhi,x0,toll,nmax)
%
% [xvect,it] = bfgs(Phi,GradPhi,x0,toll,nmax)
%
% Algoritmo BFGS per la risoluzione di problemi di ottimizzazione per
% funzioni in 2-dimensioni (n=2)
% 
% Utilizza il metodo di backtracking per aggiornare lo step-size alpha_k
% implementato in backtracking.m con parametri fissati
%
% Parametri di ingresso:
% Phi       (handle function) Funzione obiettivo
% GradPhi   (handle function) Gradiente della funzione obiettivo
% x0        (double, vettore) Guess iniziale
% toll      (double) tolleranza sulla norma del gradiente
% nmax      (int) numero massimo di iterazioni
%
% Parametri di uscita:
% xvect     (double, matrice) matrice che, su ogni colonna, contiene 
%               la soluzione approssimata a ogni iterata
% it        (int) numero di iterazioni
%
%                                         Politecnico di Milano, 03/04/2025
%

% Inizializzo variabili per ciclo iterativo
it  = 0;
err = toll+1;

% Inizializzo la matrice soluzione
xvect = x0;

% Calcolo la direzione iniziale di discesa
I  = eye(length(x0));
Bk = I;
d  = - Bk * GradPhi(x0(1), x0(2));

% Parametri backtracking
c1_bt   = 1e-4;
rho_bt  = 0.5;
nmax_bt = 10;

% Ciclo while per calcolo x^{it+1} (it = 0, 1, ...)
while it < nmax && err > toll
    
    % Calcolo il nuovo valore di alpha_k con backtracking
    [alpha_k,~] = backtracking(Phi, GradPhi, x0, d, c1_bt, rho_bt, nmax_bt);    
    
    % Calcolo il nuovo valore di x
    x_new = x0 + alpha_k*d;

    % Calcolo Bk (approssimazione inversa Hessiano di Phi)
    delta_k = x_new - x0;
    s_k     = GradPhi(x_new(1), x_new(2)) - GradPhi(x0(1), x0(2));  
    rho_k   = 1 / (s_k' * delta_k);    
    Bk      = (I - rho_k*delta_k*s_k')*Bk*(I - rho_k*s_k*delta_k') + (rho_k*delta_k)*delta_k';    
    
    % Aggiorno la direzione di discesa
    d = - Bk*GradPhi(x_new(1), x_new(2));  

    % Aggiorno quantità per metodo iterativo
    it    = it+1;
    err   = norm(GradPhi(x_new(1), x_new(2)));
    x0    = x_new;
    xvect = [xvect, x_new];
    
end

if err <= toll
     fprintf(['\nIl metodo BFGS converge in %d iterazioni \n' ...
         'al punto di minimo x = (%f,%f) \n '], it, xvect(1,end), xvect(2,end));
else
     fprintf('Il metodo BFGS non converge in %d iterazioni. \n', it)
end

end