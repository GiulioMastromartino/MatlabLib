function [fl_x, abs_err, rel_err, u] = floating_point_representation(x, beta, t, L, U)
% FLOATING_POINT_REPRESENTATION Simula la rappresentazione in un sistema F(beta, t, L, U)
% e calcola l'errore di arrotondamento.
%
%   [fl_x, abs_err, rel_err, u] = floating_point_representation(x, beta, t, L, U)
%
%   Input:
%       x:      Numero reale da rappresentare
%       beta:   Base del sistema (es. 2, 10)
%       t:      Numero di cifre della mantissa
%       L:      Estremo inferiore dell'esponente
%       U:      Estremo superiore dell'esponente
%
%   Output:
%       fl_x:    Rappresentazione floating point di x (arrotondamento al più vicino)
%       abs_err: Errore assoluto |x - fl(x)|
%       rel_err: Errore relativo |x - fl(x)| / |x|
%       u:       Unità di arrotondamento (machine epsilon teorico) = 1/2 * beta^(1-t)
%
%   Esempio:
%       [fl, err] = floating_point_representation(exp(pi), 2, 5, -6, 6)

    if x == 0
        fl_x = 0;
        abs_err = 0;
        rel_err = 0;
        u = 0.5 * beta^(1-t);
        return;
    end

    % 1. Calcolo esponente 'e' tale che beta^(-1) <= m < 1
    % x = m * beta^e
    % log_beta(x) = log_beta(m) + e
    % Poiché -1 <= log_beta(m) < 0, allora floor(log_beta(x)) = e - 1
    % Quindi e = floor(log_beta(|x|)) + 1
    
    e = floor(log(abs(x)) / log(beta)) + 1;
    
    % 2. Gestione Overflow/Underflow
    if e < L
        warning('Underflow: x troppo piccolo per F(%d, %d, %d, %d)', beta, t, L, U);
        fl_x = 0; % Flush to zero
        abs_err = abs(x);
        rel_err = 1;
        u = 0.5 * beta^(1-t);
        return;
    elseif e > U
        warning('Overflow: x troppo grande per F(%d, %d, %d, %d)', beta, t, L, U);
        fl_x = sign(x) * inf;
        abs_err = inf;
        rel_err = inf;
        u = 0.5 * beta^(1-t);
        return;
    end
    
    % 3. Calcolo mantissa normalizzata
    m = abs(x) / (beta^e);
    
    % 4. Arrotondamento mantissa a t cifre
    % Moltiplichiamo per beta^t per avere un intero
    m_scaled = m * (beta^t);
    m_rounded = round(m_scaled);
    
    % Recuperiamo mantissa arrotondata
    m_final = m_rounded / (beta^t);
    
    % Caso speciale: se l'arrotondamento ha portato m a 1 (es. 0.99... -> 1.0)
    % allora la mantissa diventa 1/beta e l'esponente aumenta di 1.
    % Esempio base 10: 0.999 -> 1.0. Normalizzato è 0.1 * 10^1.
    if m_final >= 1
        m_final = m_final / beta;
        e = e + 1;
        if e > U
             warning('Overflow dopo arrotondamento.');
             fl_x = sign(x) * inf;
             return;
        end
    end
    
    % 5. Ricostruzione
    fl_x = sign(x) * m_final * (beta^e);
    
    % 6. Errori
    abs_err = abs(x - fl_x);
    rel_err = abs_err / abs(x);
    
    % 7. Unità di arrotondamento (bound teorico errore relativo)
    u = 0.5 * beta^(1-t);

end
