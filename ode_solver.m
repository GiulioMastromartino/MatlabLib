function [t, u, info] = ode_solver(f_fun, t_span, y0, method, opts)
%ODE_SOLVER Solutore universale per ODE y' = f(t,y)
%   Supporta metodi Espliciti e Impliciti (tramite Newton interno).
%
%   INPUT:
%   f_fun   - Function handle @(t,y).
%   t_span  - [t0, tf]
%   y0      - Condizione iniziale (vettore colonna)
%   method  - 'fe' (Eulero Avanti), 'be' (Eulero Indietro), 'heun', 'cn' (Crank-Nicolson), 'rk4'
%   opts    - Struct con: .h (passo), .tol, .max_it, .Jf (Jacobiano per metodi impliciti)

    t0 = t_span(1); tf = t_span(2);
    if isfield(opts, 'h'), h = opts.h; else, h = (tf-t0)/100; end
    N = ceil((tf - t0) / h);
    
    y0 = y0(:);
    neq = length(y0);
    
    t = linspace(t0, tf, N+1);
    u = zeros(neq, N+1);
    u(:, 1) = y0;
    
    % Parametri Newton
    tol = 1e-6; if isfield(opts, 'tol'), tol = opts.tol; end
    max_it = 20; if isfield(opts, 'max_it'), max_it = opts.max_it; end
    
    tot_newton = 0;
    
    for n = 1:N
        tn = t(n); tn1 = t(n+1);
        yn = u(:, n);
        
        switch lower(method)
            case {'fe', 'eulero_avanti'}
                u(:, n+1) = yn + h * f_fun(tn, yn);
                
            case 'heun'
                k1 = f_fun(tn, yn);
                y_pred = yn + h * k1;
                k2 = f_fun(tn1, y_pred);
                u(:, n+1) = yn + (h/2) * (k1 + k2);
                
            case 'rk4'
                k1 = f_fun(tn, yn);
                k2 = f_fun(tn + h/2, yn + 0.5*h*k1);
                k3 = f_fun(tn + h/2, yn + 0.5*h*k2);
                k4 = f_fun(tn + h, yn + h*k3);
                u(:, n+1) = yn + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
                
            case {'be', 'eulero_indietro'}
                % F(y) = y - yn - h*f(tn1, y) = 0
                eq = @(y) y - yn - h * f_fun(tn1, y);
                Jf_eq = [];
                if isfield(opts, 'Jf'), Jf_eq = @(y) eye(neq) - h * opts.Jf(tn1, y); end
                [u(:, n+1), it] = local_newton(eq, Jf_eq, yn, tol, max_it);
                tot_newton = tot_newton + it;
                
            case {'cn', 'crank_nicolson'}
                % F(y) = y - yn - h/2*(f(tn,yn) + f(tn1,y)) = 0
                fn = f_fun(tn, yn);
                eq = @(y) y - yn - (h/2)*(fn + f_fun(tn1, y));
                Jf_eq = [];
                if isfield(opts, 'Jf'), Jf_eq = @(y) eye(neq) - 0.5*h * opts.Jf(tn1, y); end
                [u(:, n+1), it] = local_newton(eq, Jf_eq, yn, tol, max_it);
                tot_newton = tot_newton + it;
                
            otherwise
                error('Metodo non supportato');
        end
    end
    info.newton_iters = tot_newton;
end

function [y, it] = local_newton(fun, jac, y0, tol, max_it)
    y = y0;
    for it = 1:max_it
        F = fun(y);
        if norm(F, inf) < tol, return; end
        
        if ~isempty(jac)
            J = jac(y);
        else
            % Finite Diff Jacobian
            n = length(y); J = zeros(n);
            step = 1e-8;
            for j=1:n
                y_p = y; y_p(j) = y_p(j)+step;
                J(:,j) = (fun(y_p)-F)/step;
            end
        end
        y = y - J\F;
    end
end
