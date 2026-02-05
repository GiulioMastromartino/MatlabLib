function [y_i] = cubicspline(x_nod,y_nod,x_i)
%
% [y_i] = cubicspline(x_nod, y_nod, x_i)
% 
% Interpolazione con Spline Cubica Naturale
%
% Parametri di ingresso:
% x_nod, y_nod    Ascisse e ordinate dei nodi di interpolazione
% x_i             Ascisse in cui si vuole valutare la spline
%
% Parametri di uscita:
% y_i             Valori (ordinate) della spline in corrispondenza di x_i
%
%                                         Politecnico di Milano, 12/06/2024
%
    

SP_nat = csape(x_nod, y_nod, 'variational');

y_i    = fnval(SP_nat, x_i);
