function [xvect,it] = newton(fun,dfun,x0,toll,nmax,mol)

% [xvect,it] = newton(fun,dfun,x0,toll,nmax,mol)
%
% Metodo di Newton per l'approssimazione di uno zero (radice) della 
% funzione fun. Se la molteplicità (mol) dello zero da approssimare è 
% specificata, allora viene applicato il metodo di Newton modificato. 
% Criterio d'arresto basato sulla differenza tra due iterate successive.
%
% Parametri di ingresso:
% fun, dfun   Inline functions contenenti la funzione e la sua derivata
% x0          Iterata iniziale
% toll        Tolleranza sul criterio d'arresto (diff. iter. succ.)
% nmax        Numero massimo di iterazioni
% mol         Molteplicità assegnata allo zero. Se assegnato, permette di 
%             applicare il metodo di Newton modificato. Default: mol = 1
%
% Parametri di uscita:
% xvect      Vettore contenente tutte le iterate calcolate (l'ultima 
%            componente e' lo zero approssimato)
% it         Numero di iterazioni effettuate
%
%                                         Politecnico di Milano, 04/04/2024
%

if nargin < 6
    mol = 1;
end
if nargin < 5
    nmax = 100;
end
if nargin < 4
    toll = 1e-6;
end

err = toll + 1;
it = 0;
xvect = x0;
xv = x0;

while (it < nmax && err > toll)
   dfx = dfun(xv);
   if dfx == 0
      error(' Arresto per azzeramento di dfun');
   else
      xn = xv - mol*fun(xv)/dfx;
      err = abs(xn-xv);
      xvect = [xvect; xn];
      it = it+1;
      xv = xn;
   end
end

if (it < nmax)
    fprintf(' Convergenza al passo k : %d \n',it);
else
    fprintf(' E` stato raggiunto il numero massimo di passi k : %d \n',it);
end
fprintf(' Radice calcolata       : %-12.8f \n', xvect(end));
end