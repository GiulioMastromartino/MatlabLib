function u = bvp_mini_solver(N, coeffs, bc, method)
%BVP_MINI_SOLVER Risolve -u'' + p u' + q u = f su [0,1]
%   coeffs: [p, q, f] (costanti)
%   bc: [u(0), u(1)] (Dirichlet)
%   method: 'centered' (default) o 'upwind'

    if nargin < 4, method = 'centered'; end
    p = coeffs(1); q = coeffs(2); f = coeffs(3);
    ua = bc(1); ub = bc(2);
    
    h = 1/N; 
    n = N-1; % Nodi interni
    
    % Diff Seconda centrata: (-1, 2, -1)/h^2 -> -u'' approx (-u_{i-1} + 2u_i - u_{i+1})/h^2
    c2_main = 2/h^2;
    c2_side = -1/h^2;
    
    % Diff Prima
    if strcmp(method, 'upwind')
        if p >= 0
            % Backward: (u_i - u_{i-1})/h
            c1_main = 1/h;
            c1_low  = -1/h;
            c1_up   = 0;
        else
            % Forward: (u_{i+1} - u_i)/h
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
    
    % Assemblaggio con sparse triplet (più sicuro di spdiags)
    I = zeros(3*n, 1);
    J = zeros(3*n, 1);
    V = zeros(3*n, 1);
    idx = 0;
    
    b = f * ones(n,1);
    
    for i=1:n
        % Diagonale principale (u_i)
        idx = idx + 1;
        I(idx) = i; J(idx) = i;
        V(idx) = c2_main + p*c1_main + q;
        
        % Sotto-diagonale (u_{i-1})
        coeff_low = c2_side + p*c1_low;
        if i > 1
            idx = idx + 1;
            I(idx) = i; J(idx) = i-1;
            V(idx) = coeff_low;
        else
            % Condizione al bordo sx
            b(i) = b(i) - coeff_low * ua;
        end
        
        % Sopra-diagonale (u_{i+1})
        coeff_up = c2_side + p*c1_up;
        if i < n
            idx = idx + 1;
            I(idx) = i; J(idx) = i+1;
            V(idx) = coeff_up;
        else
            % Condizione al bordo dx
            b(i) = b(i) - coeff_up * ub;
        end
    end
    
    % Trim vectors
    I = I(1:idx); J = J(1:idx); V = V(1:idx);
    
    A = sparse(I, J, V, n, n);
    
    u_int = A\b;
    u = [ua; u_int; ub];
end