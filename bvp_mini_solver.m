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
    low = zeros(n-1,1); 
    up = zeros(n-1,1);
    b = f * ones(n,1);
    
    for i=1:n
        % Diff Seconda: (-1, 2, -1)/h^2
        c2 = [-1, 2, -1]/h^2;
        
        % Diff Prima
        if strcmp(method, 'upwind')
            if p>0, c1 = [-1, 1, 0]/h; else, c1 = [0, -1, 1]/h; end
        else
            c1 = [-0.5, 0, 0.5]/h;
        end
        
        row = c2 + p*c1;
        row(2) = row(2) + q;
        
        % Riempimento vettori diagonali
        % main: elemento diagonale A(i,i)
        main(i) = row(2);
        
        % low: elemento sottodiagonale A(i, i-1) -> corrisponde a coeff di u_{i-1} (row(1))
        if i>1, low(i-1) = row(1); end
        
        % up: elemento sopradiagonale A(i, i+1) -> corrisponde a coeff di u_{i+1} (row(3))
        if i<n, up(i) = row(3); end
        
        % Termini noti (Condizioni al Bordo)
        % Se i=1, u_{i-1} è u_0 (ua)
        if i==1, b(i) = b(i) - row(1)*ua; end
        % Se i=n, u_{i+1} è u_N (ub)
        if i==n, b(i) = b(i) - row(3)*ub; end
    end
    
    % Correzione per spdiags:
    % Le colonne devono essere lunghe n.
    % La colonna -1 (sottodiagonale) deve avere A(i, i-1) in posizione i.
    % La colonna +1 (sopradiagonale) deve avere A(i, i+1) in posizione i.
    
    low_col = [0; low]; % Pad all'inizio, così low(1) finisce in indice 2
    up_col = [up; 0];   % Pad alla fine
    
    A = spdiags([low_col main up_col], -1:1, n, n);
    
    u_int = A\b;
    u = [ua; u_int; ub];
end