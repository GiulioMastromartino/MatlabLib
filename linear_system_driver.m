function [x, info] = linear_system_driver(A, b, cfg)
%LINEAR_SYSTEM_DRIVER Risolve Ax=b con metodi iterativi
%   cfg.type: 'jacobi', 'gs', 'sor', 'richardson', 'grad'
%   cfg.tol, cfg.max_it, cfg.x0
%   cfg.mu (per SOR), cfg.P (per Richardson)

    n = length(b);
    if ~isfield(cfg, 'x0'), x = zeros(n,1); else, x = cfg.x0; end
    tol = 1e-6; if isfield(cfg, 'tol'), tol = cfg.tol; end
    max_it = 1000; if isfield(cfg, 'max_it'), max_it = cfg.max_it; end
    
    norm_b = norm(b); if norm_b==0, norm_b=1; end
    res_hist = zeros(max_it, 1);
    
    % Setup
    method = lower(cfg.type);
    if strcmp(method, 'jacobi'), D_inv = 1./diag(A); end
    if strcmp(method, 'richardson') 
        if isfield(cfg, 'P'), [Lp, Up] = lu(cfg.P); else, Lp=[]; Up=[]; end
    end
    
    for iter = 1:max_it
        r = b - A*x;
        res_norm = norm(r)/norm_b;
        res_hist(iter) = res_norm;
        
        if res_norm < tol, break; end
        
        x_old = x;
        
        switch method
            case 'jacobi'
                x = x_old + D_inv .* r;
            case 'gs'
                for i=1:n
                    sigma = A(i, [1:i-1, i+1:n]) * x([1:i-1, i+1:n]);
                    x(i) = (b(i) - sigma)/A(i,i);
                end
            case 'sor'
                mu = cfg.mu;
                for i=1:n
                    sigma_new = A(i, 1:i-1)*x(1:i-1);
                    sigma_old = A(i, i+1:n)*x_old(i+1:n);
                    x_gs = (b(i) - sigma_new - sigma_old)/A(i,i);
                    x(i) = (1-mu)*x_old(i) + mu*x_gs;
                end
            case 'richardson'
                if ~isempty(Lp), z = Up\(Lp\r); else, z = r; end
                alpha = 1; % Default dinamico o fisso
                x = x_old + alpha*z;
            case 'grad' % Steepest Descent
                Ar = A*r;
                alpha = (r'*r)/(r'*Ar);
                x = x_old + alpha*r;
        end
    end
    info.iters = iter;
    info.res = res_norm;
    info.hist = res_hist(1:iter);
end
