function [xmin, xmax] = floating_point_range(beta, t, L, U)
% FLOATING_POINT_RANGE Calcola il range di numeri rappresentabili in F(beta, t, L, U)
%
%   [xmin, xmax] = floating_point_range(beta, t, L, U)
%
%   Input:
%       beta: base del sistema numerico (intero >= 2)
%       t: numero di cifre della mantissa (precisione)
%       L: limite inferiore dell'esponente
%       U: limite superiore dell'esponente
%
%   Output:
%       xmin: il più piccolo numero positivo rappresentabile normalizzato
%       xmax: il più grande numero positivo rappresentabile
%
%   Esempio:
%       [xmin, xmax] = floating_point_range(2, 53, -1022, 1023) % IEEE 754 Double

    % Verifica input minimi
    if nargin < 4
        error('Richiesti 4 argomenti: beta, t, L, U');
    end

    % Calcolo xmin (normalizzato)
    % La mantissa normalizzata minima è 1.00...0 = 1 * beta^0 = 1 (in base beta, interpretato come m * beta^-t no?)
    % Definizione standard: x = +/- m * beta^(e-t) con beta^(t-1) <= m <= beta^t - 1
    % Oppure x = +/- .d1 d2 ... dt * beta^e con d1 != 0.
    % La definizione più comune nei corsi di calcolo numerico (e usata nelle soluzioni precedenti) è:
    % x_min = beta^(L-1)
    
    xmin = beta^(L-1);
    
    % Calcolo xmax
    % La mantissa massima è 0. (beta-1) (beta-1) ... (beta-1)  = 1 - beta^(-t)
    % L'esponente massimo è U.
    % x_max = (1 - beta^(-t)) * beta^U
    
    xmax = (1 - beta^(-t)) * beta^U;

end
