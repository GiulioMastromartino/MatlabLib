function [t, u, info] = rk_generic(Butcher, f_fun, t_span, y0, opts)
%RK_GENERIC Solutore Runge-Kutta universale (Esplicito e Implicito)
%   Risolve y' = f(t,y) dato un Tableau di Butcher generico.
%
%   INPUT:
%   Butcher - Struct con campi:
%             .A (Matrice stadi s x s)
%             .b (Vettore pesi 1 x s)
%             .c (Vettore nodi 1 x s)
%   f_fun   - Function handle @(t,y) -> colonna
%   t_span  - [t0, tf]
%   y0      - Condizione iniziale
%   opts    - Struct opzionale:
%             .h (passo costante) o .N (numero passi)
%             .tol, .max_it (per Newton negli stadi impliciti)
%             .Jf (Jacobiano @(t,y) per Newton veloce)
%
%   OUTPUT:
%   t       - Vettore tempi
%   u       - Matrice soluzione (dim x N+1)
%   info    - Statistiche (iterazioni Newton)

    % --- 1. Setup ---
    A_rk = Butcher.A;
    b_rk = Butcher.b(:).'; % Riga
    c_rk = Butcher.c(:);   % Colonna
    s = length(c_rk);      % Numero stadi
    
    t0 = t_span(1); tf = t_span(2);
    if isfield(opts, 'h'), h = opts.h; else, h = (tf-t0)/100; end
    N = ceil((tf - t0) / h);
    
    y0 = y0(:);
    neq = length(y0);
    
    t = linspace(t0, tf, N+1);
    u = zeros(neq, N+1);
    u(:, 1) = y0;
    
    % Configurazione Newton (per stadi impliciti)
    tol = 1e-6; if isfield(opts, 'tol'), tol = opts.tol; end
    max_it = 20; if isfield(opts, 'max_it'), max_it = opts.max_it; end
    
    tot_newton = 0;
    
    % Analisi Metodo: Esplicito se A è triangolare stretta inferiore
    is_explicit = isequal(triu(A_rk), zeros(s));
    
    % --- 2. Time Loop ---
    for n = 1:N
        tn = t(n);
        yn = u(:, n);
        
        % K matrix: ogni colonna j è il vettore k_j (dim neq)
        K = zeros(neq, s);
        
        for i = 1:s
            % Calcolo termine noto dagli stadi precedenti (j < i)
            % Sum_prev = sum_{j=1}^{i-1} a_{ij} * k_j
            if i > 1
                Sum_prev = K(:, 1:i-1) * A_rk(i, 1:i-1).';
            else
                Sum_prev = zeros(neq, 1);
            end
            
            t_stage = tn + c_rk(i) * h;
            a_ii = A_rk(i, i);
            
            if a_ii == 0
                % --- STADIO ESPLICITO ---
                % k_i = f(tn + ci*h, yn + h * Sum_prev)
                Y_i = yn + h * Sum_prev;
                K(:, i) = f_fun(t_stage, Y_i);
                
            else
                % --- STADIO IMPLICITO ---
                % Risolvere per k_i: k_i = f(t_stage, yn + h*Sum_prev + h*a_ii*k_i)
                % Residuo: R(k) = k - f(..., Y_base + h*a_ii*k) = 0
                
                Y_base = yn + h * Sum_prev;
                
                % Guess iniziale per Newton (uso k_{i-1} o 0)
                if i > 1, k_guess = K(:, i-1); else, k_guess = zeros(neq, 1); end
                
                % Equazione residuo
                eq_fun = @(k) k - f_fun(t_stage, Y_base + h * a_ii * k);
                
                % Jacobiano residuo: I - h * a_ii * Jf
                if isfield(opts, 'Jf') && ~isempty(opts.Jf)
                    jac_fun = @(k) eye(neq) - h * a_ii * opts.Jf(t_stage, Y_base + h * a_ii * k);
                else
                    jac_fun = []; % Newton userà differenze finite
                end
                
                [k_sol, it] = local_newton_solver(eq_fun, jac_fun, k_guess, tol, max_it);
                K(:, i) = k_sol;
                tot_newton = tot_newton + it;
            end
        end
        
        % --- 3. Assemblaggio Finale ---
        % y_{n+1} = y_n + h * sum(b_i * k_i)
        u(:, n+1) = yn + h * (K * b_rk.');
    end
    
    info.newton_iters = tot_newton;
    info.is_explicit = is_explicit;
end

% --- Helper: Newton Solver Locale ---
function [x, it] = local_newton_solver(fun, jac, x0, tol, max_it)
    x = x0;
    for it = 1:max_it
        F = fun(x);
        if norm(F, inf) < tol, return; end
        
        if ~isempty(jac)
            J = jac(x);
        else
            % Jacobiano numerico semplice
            n = length(x); J = zeros(n);
            step = 1e-8;
            for j=1:n
                xp = x; xp(j) = xp(j) + step;
                J(:,j) = (fun(xp) - F) / step;
            end
        end
        
        delta = -J \ F;
        x = x + delta;
    end
    warning('RK Newton stage non converge');
end
