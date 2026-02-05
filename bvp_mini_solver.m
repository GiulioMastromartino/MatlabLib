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
    main = zeros(n,1); low = zeros(n-1,1); up = zeros(n-1,1);
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
        
        main(i) = row(2);
        if i>1, low(i-1) = row(1); end
        if i<n, up(i) = row(3); end
        
        % BC
        if i==1, b(i) = b(i) - row(1)*ua; end
        if i==n, b(i) = b(i) - row(3)*ub; end
    end
    
    A = spdiags([low main up], -1:1, n, n);
    u_int = A\b;
    u = [ua; u_int; ub];
end
