function u = bvp_mini_solver(N, coeffs, bc, method)
%BVP_MINI_SOLVER Risolve -u'' + p u' + q u = f su [0,1]
%   coeffs: [p, q, f] (costanti)
%   bc: [u(0), u(1)] (Dirichlet)
%   method: 'centered' (default) o 'upwind'

    if nargin < 4, method = 'centered'; end
    p = coeffs(1); q = coeffs(2); f = coeffs(3);
    ua = bc(1); ub = bc(2);
    
    h = 1/N; n = N-1;
    
    % Costruzione matrice A (tridiagonale)
    main = zeros(n,1); 
    low = zeros(n,1); % Lunghezza n per spdiags
    up = zeros(n,1);  % Lunghezza n per spdiags
    b = f * ones(n,1);
    
    % Diff Seconda centrata: (-1, 2, -1)/h^2
    % Attenzione: il problema è -u'', quindi i segni sono opposti a u''
    % -u'' approx -(u_{i+1} - 2u_i + u_{i-1})/h^2 = (-u_{i+1} + 2u_i - u_{i-1})/h^2
    c2_main = 2/h^2;
    c2_side = -1/h^2;
    
    % Diff Prima
    if strcmp(method, 'upwind')
        % Upwind per u'
        if p >= 0
            % Backward diff: (u_i - u_{i-1})/h
            c1_main = 1/h;
            c1_low  = -1/h;
            c1_up   = 0;
        else
            % Forward diff: (u_{i+1} - u_i)/h
            c1_main = -1/h;
            c1_low  = 0;
            c1_up   = 1/h;
        end
    else
        % Centered: (u_{i+1} - u_{i-1})/(2h)
        c1_main = 0;
        c1_low  = -1/(2*h);
        c1_up   = 1/(2*h);
    end
    
    for i=1:n
        % Indici per A
        % Diagonale principale
        val_main = c2_main + p*c1_main + q;
        main(i) = val_main;
        
        % Sotto-diagonale (coeff di u_{i-1})
        val_low = c2_side + p*c1_low;
        if i > 1
            low(i) = val_low; % spdiags vuole il valore alla riga i
        else
            % Condizione al bordo u_0 = ua
            b(i) = b(i) - val_low * ua;
        end
        
        % Sopra-diagonale (coeff di u_{i+1})
        val_up = c2_side + p*c1_up;
        if i < n
            up(i) = val_up; % spdiags vuole il valore alla riga i
        else
            % Condizione al bordo u_N = ub
            b(i) = b(i) - val_up * ub;
        end
    end
    
    % Assemblaggio con spdiags
    % low è la diagonale -1. spdiags si aspetta che low(i) sia l'elemento (i, i-1)
    % Ma spdiags shiftati:
    % diag -1: prende gli elementi da indice 2 a n. L'elemento in (2,1) è in pos 2 del vettore.
    % diag +1: prende gli elementi da indice 1 a n-1. L'elemento in (1,2) è in pos 1 del vettore.
    
    % Verifichiamo la convenzione di spdiags:
    % B = [low, main, up]
    % d = [-1, 0, 1]
    % La colonna 'low' corrisponde alla diagonale -1. L'elemento B(i,1) va in (i, i-1).
    % Quindi B(2,1) va in (2,1). B(1,1) viene ignorato.
    % Il nostro ciclo riempie low(i) pensando alla riga i.
    % Quindi low(2) è il coeff di u_1 nella riga 2 -> corretto.
    
    % La colonna 'up' corrisponde alla diagonale +1. L'elemento B(i,3) va in (i, i+1).
    % Quindi B(1,3) va in (1,2).
    % Il nostro ciclo riempie up(i) pensando alla riga i.
    % Quindi up(1) è il coeff di u_2 nella riga 1 -> corretto.
    
    A = spdiags([low, main, up], -1:1, n, n);
    
    u_int = A\b;
    u = [ua; u_int; ub];
end