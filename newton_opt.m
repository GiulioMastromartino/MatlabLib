function [xvect,it] = newton_opt(GradPhi,HessPhi,x0,toll,nmax)
%
% [xvect,it] = newton_opt(GradPhi,HessPhi,x0,toll,nmax)
%
% Algoritmo del metodo di Newton per la risoluzione di problemi di ottimizzazione per
% funzioni in 2-dimensioni (n=2) 
% 
% Metodo di Newton esatto, quindi con step-size \alpha_k = 1
%
% Parametri di ingresso:
% GradPhi     (handle function) Gradiente della funzione obiettivo
% HessPhi     (handle function) Hessiano della funzione obiettivo
% x0          (double, vettore) Guess iniziale
% toll        (double) tolleranza sulla norma del gradiente
% nmax        (int) numero massimo di iterazioni
%
% Parametri di uscita:
% xvect        (double, matrice) matrice che, su ogni colonna, contiene 
%               la soluzione approssimata a ogni iterata
% it           (int) numero di iterazioni
%
%                                         Politecnico di Milano, 03/04/2025
%


% Inizializzo variabili per ciclo iterativo
it  = 0;
err = toll+1;

% Inizializzo la matrice soluzione
xvect = x0;

% Calcolo la direzione iniziale di discesa
d = - HessPhi(x0(1),x0(2)) \ GradPhi(x0(1),x0(2));

% Definiamo lo step-size alpha_k
alpha_k = 1;

% Ciclo while per calcolo x^{it+1} (it = 0, 1, ...)
while it < nmax && err > toll
    
    % Calcolo il nuovo valore di x
    x_new = x0 + alpha_k*d;
    
    % Aggiorno la direzione di discesa
    d = - HessPhi(x_new(1), x_new(2)) \ GradPhi(x_new(1), x_new(2));

    % Aggiorno quantità per metodo iterativo
    it    = it+1;
    err   = norm(GradPhi(x_new(1), x_new(2)));
    x0    = x_new;
    xvect = [xvect, x_new];

end

if err <= toll
     fprintf(['\nIl metodo di Newton converge in %d iterazioni \n' ...
         'al punto di minimo x = (%f,%f)^T \n '], it, xvect(1,end), xvect(2,end));
else
     fprintf('Il metodo di Newton non converge in %d iterazioni. \n', it)
end

end