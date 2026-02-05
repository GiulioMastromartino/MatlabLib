function [xmin, xmax, mach_eps] = floating_point_range(varargin)
% FLOATING_POINT_RANGE Calcola range e epsilon macchina per sistemi floating point
%
% Utilizzo:
%   [xmin, xmax, eps] = floating_point_range(beta, t, L, U)
%   [xmin, xmax, eps] = floating_point_range(type_string)
%   [xmin, xmax, eps] = floating_point_range(struct_config)
%
% Input:
%   Caso 1 (4 argomenti):
%       beta: base (es. 2, 10)
%       t: cifre significative (precisione)
%       L: esponente minimo
%       U: esponente massimo
%
%   Caso 2 (1 argomento stringa):
%       'double': IEEE 754 Double Precision (binary64)
%       'single': IEEE 754 Single Precision (binary32)
%       'half'  : IEEE 754 Half Precision (binary16)
%
%   Caso 3 (1 argomento struct):
%       struttura con campi .beta, .t, .L, .U
%
% Output:
%   xmin: Minimo numero positivo normalizzato (realmin)
%   xmax: Massimo numero positivo finito (realmax)
%   mach_eps: Epsilon macchina (distanza tra 1 e il successivo numero floating point)
%
% Esempi:
%   [min_d, max_d, eps_d] = floating_point_range('double');
%   [min_c, max_c, eps_c] = floating_point_range(2, 8, -20, 20);

    % Gestione Input
    if nargin == 4
        beta = varargin{1};
        t = varargin{2};
        L = varargin{3};
        U = varargin{4};
    elseif nargin == 1
        arg = varargin{1};
        if ischar(arg) || isstring(arg)
            switch lower(arg)
                case 'double'
                    % IEEE 754 Double (64 bit)
                    % Mantissa: 53 bits (52 espliciti + 1 implicito)
                    % Exp: 11 bits. Bias 1023. e_min = -1022, e_max = 1023
                    % Modello frazionario (0.m x beta^e):
                    % L = -1021, U = 1024
                    beta = 2; t = 53; L = -1021; U = 1024;
                case 'single'
                    % IEEE 754 Single (32 bit)
                    % Mantissa: 24 bits (23 espliciti + 1 implicito)
                    % Exp: 8 bits. Bias 127. e_min = -126, e_max = 127
                    % Modello frazionario: L = -125, U = 128
                    beta = 2; t = 24; L = -125; U = 128;
                case 'half'
                    % IEEE 754 Half (16 bit)
                    % Mantissa: 11 bits (10 espliciti + 1 implicito)
                    % Exp: 5 bits. Bias 15. e_min = -14, e_max = 15
                    % Modello frazionario: L = -13, U = 16
                    beta = 2; t = 11; L = -13; U = 16;
                case 'quad'
                    % IEEE 754 Quad (128 bit)
                    % Mantissa: 113 bits
                    % Exp: 15 bits. Bias 16383. e_min = -16382, e_max = 16383
                    % Modello frazionario: L = -16381, U = 16384
                    beta = 2; t = 113; L = -16381; U = 16384;
                otherwise
                    error('Tipo non supportato. Usare "double", "single", "half" o "quad".');
            end
        elseif isstruct(arg)
            beta = arg.beta;
            t = arg.t;
            L = arg.L;
            U = arg.U;
        else
            error('Input non valido. Richiesti 4 scalari, una stringa o una struct.');
        end
    else
        error('Numero di argomenti non valido.');
    end

    % Formule basate sul modello a mantissa frazionaria normalizzata
    % x = +/- 0.d1...dt * beta^e
    % con d1 != 0 (quindi m in [1/beta, 1))
    % e in [L, U]
    
    % Minimo numero positivo (normalizzato)
    % m_min = 1/beta
    % e_min = L
    % xmin = (1/beta) * beta^L = beta^(L-1)
    xmin = beta^(L-1);
    
    % Massimo numero positivo
    % m_max = 1 - beta^(-t)  (es. 0.99...9)
    % e_max = U
    % xmax = (1 - beta^(-t)) * beta^U
    xmax = (1 - beta^(-t)) * beta^U;
    
    % Epsilon macchina
    % Distanza tra 1 e il numero successivo.
    % 1 nel modello frazionario è rappresentato come (1/beta) * beta^1  (se beta=2, 0.5 * 2^1)
    % Oppure come 1.0 * beta^0 nel modello a mantissa intera.
    % La definizione standard di eps è beta^(1-t).
    % Esempio double (t=53): 2^(1-53) = 2^-52 = 2.22e-16. Corretto.
    mach_eps = beta^(1-t);

end
