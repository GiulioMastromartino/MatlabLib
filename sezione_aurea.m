function [xv,k,errv] = sezione_aurea(f,a,b,tol,nmax)
% 
% [xv,k,errv] = sezione_aurea(f,a,b,tol,nmax)
% 
% Metodo della sezione aurea per approssimare il punto di minimo di una 
% funzione scalare (n=1) f : [a,b] -> R
%
% Parametri di ingresso:
% f      funzione obiettivo
% a      estremo inferiore dell'intervallo
% b      estremo superiore dell'intervallo
% tol    tolleranza criterio d'arresto basato su stimatore dell'errore
% nmax   numero massimo di iterazioni ammesse
%
% Parametri di uscita:
% xv      vettore delle approssimazioni x_k del punto di minimo 
% k       numero di iterazioni eseguite
% errv    vettore delle stime dell'errore e_k
%
%                                         Politecnico di Milano, 03/04/2025
%


k = 0;
phi = (1 + sqrt(5))/2; % rapporto aureo
x = (a + b)/2;
% Vettore delle iterate
xv = x;
err = (b - a)/2;
% Vettore delle stime dell'errore
errv = err;
while (k < nmax && err > tol)
    k = k + 1;
    c = a + (b - a)/(phi + 1);
    d = a + (b - a)/phi;
    if(f(c) > f(d))
        a = c;
    else
        b = d;
    end
    x = (a + b)/2;
    xv = [xv, x];

    err = (b-a)/2;
    errv = [errv, err];
end

if (err <= tol)
     fprintf(['Il metodo della sezione aurea converge in %d iterazioni \n' ...
         ' al punto di minimo x = %f \n '], k, x);
else
     fprintf('Il metodo della sezione aurea non converge in %d iterazioni. \n', k)
end
end

