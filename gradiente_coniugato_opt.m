function [xvect,it] = gradiente_coniugato_opt(Phi,GradPhi,flag_metodo,x0,toll,nmax)
%
% [xvect,it] = gradiente_coniugato_opt(Phi,GradPhi,flag_metodo,x0,toll,nmax)
%
% Algoritmo del gradiente e gradiente coniugato per la risoluzione di 
% problemi di ottimizzazione per funzioni in 2-dimensioni (n=2) 
%
% Tramite il terzo input (flag_metodo) è possibile scegliere se usare il 
% metodo del gradiente o il metodo del gradiente coniugato (con parametro 
% \beta_k di Fletcher-Reeves)
%
% Utilizza il metodo di backtracking per aggiornare lo step-size alpha_k
% implementato in backtracking.m con parametri fissati
%
% Parametri di ingresso:
% Phi           (handle function) Funzione obiettivo
% GradPhi       (handle function) Gradiente della funzione obiettivo
% flag_metodo   (string) flag per scegliere tra gradiente ('G') e gradiente
%                     coniugato con parametro di Fletcher-Reeves ('FR')
% x0            (double, vettore) Guess iniziale
% toll          (double) tolleranza sulla norma del gradiente
% nmax          (int) numero massimo di iterazioni
%
% Parametri di uscita:
% xvect         (double, matrice) matrice che, su ogni colonna, contiene 
%               la soluzione approssimata a ogni iterata
% it            (int) numero di iterazioni
%
%                                         Politecnico di Milano, 03/04/2025
%


% Inizializzo variabili per ciclo iterativo
it  = 0;
err = toll+1;

% Inizializzo la matrice soluzione
xvect = x0;

% Calcolo la direzione iniziale di discesa
d = -GradPhi(x0(1), x0(2));

% Parametri backtracking
c1_bt   = 1e-4;
rho_bt  = 0.3;
nmax_bt = 10;

% Ciclo while per calcolo x^{it+1} (it = 0, 1, ...)
while it < nmax && err > toll

    % Calcolo il nuovo valore di alpha_k
    [alpha_k,~] = backtracking(Phi, GradPhi, x0, d, c1_bt, rho_bt, nmax_bt);   
        
    % Calcolo il nuovo valore di x
    x_new = x0 + alpha_k*d;

    % Swtich per scelta metodo gradiente o gradiente coniugato
    switch flag_metodo
        case 'G' % gradiente
            beta_k = 0;
        case 'FR' % gradiente coniugato (Fletcher-Reeves)
            beta_k = norm(GradPhi(x_new(1),x_new(2)))^2 / norm(GradPhi(x0(1),x0(2)))^2;
        otherwise
            error('metodo CG non definito');
    end

    % Aggiorno la direzione di discesa (se beta_k = 0 -> gradiente)
    d = -GradPhi(x_new(1),x_new(2)) + beta_k*d;  

    % Aggiorno quantità per metodo iterativo
    it    = it+1;
    err   = norm(GradPhi(x_new(1), x_new(2)));
    x0    = x_new;
    xvect = [xvect, x_new];
    
end

switch flag_metodo
    case 'G' % gradiente
        if err <= toll
             fprintf(['\n Il metodo del gradiente converge in %d iterazioni \n' ...
                 'al punto di minimo x = (%f,%f)^T \n '], it, xvect(1,end), xvect(2,end));
        else
             fprintf('\n Il metodo del gradiente non converge in %d iterazioni. \n', it)
        end
        
    case 'FR' % gradiente coniugato (Fletcher-Reeves)
        if err <= toll
             fprintf(['\n Il metodo del gradiente coniugato (Fletcher-Reeves) converge in %d iterazioni \n' ...
                 'al punto di minimo x = (%f,%f)^T \n '], it, xvect(1,end), xvect(2,end));
        else
             fprintf('\n Il metodo del gradiente coniugato (Fletcher-Reeves) non converge in %d iterazioni. \n',  it)
        end    
end

end