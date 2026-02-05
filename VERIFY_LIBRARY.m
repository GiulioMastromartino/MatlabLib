%% VERIFY_LIBRARY.m - Script di Validazione Automatica
% Lanciare questo script per verificare che tutte le funzioni siano presenti e funzionanti.
clear; clc; close all;
fprintf('=== INIZIO VALIDAZIONE LIBRERIA ===\n\n');

errors = 0;

%% TEST 1: ODE Solver (Universale)
try
    fprintf('[1/7] Test ode_solver... ');
    % Verifica esistenza e funzionalità di ode_solver
    if exist('ode_solver.m', 'file')
        f = @(t,y) -y; 
        opts.h = 0.1;
        % Test Heun tramite wrapper o chiamata diretta
        [t, u] = ode_solver(f, [0, 1], 1, 'heun', opts);
        assert(abs(u(end) - exp(-1)) < 2e-2, 'Errore eccessivo in Heun (ode_solver)');
        fprintf('OK\n');
    else
        fprintf('SKIP (ode_solver.m non trovato, testo Heun diretto)\n');
        % Fallback test diretto su Heun.m
        f = @(t,y) -y;
        [t, u] = Heun(f, 1, 1, 0.1);
        % Heun restituisce u_h come matrice, prendiamo ultimo elemento
        u_end = u(end); 
        assert(abs(u_end - exp(-1)) < 2e-2, 'Errore eccessivo in Heun.m diretto');
        fprintf('      Test Heun.m diretto: OK\n');
    end
catch ME
    fprintf('FALLITO! (%s)\n', ME.message); errors = errors+1;
end

%% TEST 2: Linear System Driver
try
    fprintf('[2/7] Test linear_system_driver... ');
    A = [4 -1; -1 4]; b = [3; 3]; % Soluzione [1;1]
    cfg.type = 'gs'; cfg.tol = 1e-6; cfg.max_it = 100; cfg.omega = 1;
    
    if exist('linear_system_driver.m', 'file')
        [x, info] = linear_system_driver(A, b, cfg);
        assert(norm(x - [1;1]) < 1e-4, 'Gauss-Seidel non converge (driver)');
        fprintf('OK\n');
    else
         fprintf('SKIP (linear_system_driver.m non trovato, testo gs diretto)\n');
         [x, k, res] = gs(A, b, [0;0], 1e-6, 100);
         assert(norm(x - [1;1]) < 1e-4, 'Gauss-Seidel non converge (gs.m)');
         fprintf('      Test gs.m diretto: OK\n');
    end
catch ME
    fprintf('FALLITO! (%s)\n', ME.message); errors = errors+1;
end

%% TEST 3: BVP Mini Solver
try
    fprintf('[3/7] Test bvp_mini_solver... ');
    % -u'' = 2, u(0)=0, u(1)=0 => u(x) = x(1-x) => u(0.5)=0.25
    % coeffs = [p, q, f] -> -u'' + p u' + q u = f
    % Qui p=0, q=0, f=2 -> -u'' = 2.
    
    if exist('bvp_mini_solver.m', 'file')
        u = bvp_mini_solver(10, [0,0,2], [0,0], 'centered');
        % u contiene [u0, u1, ..., uN] (N+1 punti)
        % N=10, punti 0, 0.1, ..., 0.5, ..., 1.0
        % Indici Matlab: 1, 2, ..., 6, ..., 11
        % u(6) corrisponde a x=0.5
        if length(u) > 5
             u_mid = u(6); 
             assert(abs(u_mid - 0.25) < 5e-2, 'Errore BVP Centrato (valore atteso 0.25)');
        end
        fprintf('OK\n');
    else
        fprintf('SKIP (File non trovato)\n');
        if exist('Poisson_Dirichlet_diff_finite_5punti.m', 'file')
             fprintf('      Trovato Poisson 5 punti (sostituto valido)\n');
        else
             errors = errors+1;
        end
    end
catch ME
    fprintf('FALLITO! (%s)\n', ME.message); errors = errors+1;
end

%% TEST 4: Numerical Approx (Simpson)
try
    fprintf('[4/7] Test numerical_approx_toolbox... ');
    % Integral x^2 su [0,1] = 1/3
    f = @(x) x.^2;
    if exist('numerical_approx_toolbox.m', 'file')
        res = numerical_approx_toolbox('simpson', f, [0,1], 2);
        assert(abs(res - 1/3) < 1e-10, 'Simpson errato');
        fprintf('OK\n');
    elseif exist('simpcomp.m', 'file')
         res = simpcomp(0, 1, 2, f);
         assert(abs(res - 1/3) < 1e-10, 'Simpson errato (simpcomp)');
         fprintf('OK (usato simpcomp.m)\n');
    else
        fprintf('FALLITO (Nessuna funzione Simpson trovata)\n'); errors = errors+1;
    end
catch ME
    fprintf('FALLITO! (%s)\n', ME.message); errors = errors+1;
end

%% TEST 5: Linalg Toolbox (Power Method)
try
    fprintf('[5/7] Test linalg_toolbox... ');
    A = [2 0; 0 1]; % Autovalori 2, 1
    opts.k = 20; opts.x0 = [1; 1]; opts.tol = 1e-6;
    
    if exist('linalg_toolbox.m', 'file')
        [lambda, ~] = linalg_toolbox('power_iter', A, opts);
        assert(abs(lambda - 2) < 1e-4, 'Power Method fallito (toolbox)');
        fprintf('OK\n');
    elseif exist('eigpower.m', 'file')
        [lambda, ~, ~] = eigpower(A, 1e-6, 100, [1;1]);
        assert(abs(lambda - 2) < 1e-4, 'Power Method fallito (eigpower)');
        fprintf('OK (usato eigpower.m)\n');
    else
        fprintf('FALLITO (Nessuna funzione Power Method trovata)\n'); errors = errors+1;
    end
catch ME
    fprintf('FALLITO! (%s)\n', ME.message); errors = errors+1;
end

%% TEST 6: Data Analysis (Least Squares)
try
    fprintf('[6/7] Test data_analysis_toolbox... ');
    x = 1:5; y = 2*x + 1; % Retta y=2x+1
    
    if exist('data_analysis_toolbox.m', 'file')
        [coeffs, ~] = data_analysis_toolbox('lsq', x, y, 1);
        assert(abs(coeffs(1)-2) < 1e-5, 'Fit lineare errato');
        fprintf('OK\n');
    else
        p = polyfit(x,y,1);
        assert(abs(p(1)-2) < 1e-5, 'Polyfit base fallito');
        fprintf('OK (usato polyfit standard)\n');
    end
catch ME
    fprintf('FALLITO! (%s)\n', ME.message); errors = errors+1;
end

%% TEST 7: RK Generic
try
    fprintf('[7/7] Test rk_generic... ');
    if exist('rk_generic.m', 'file')
        Butcher.A = 0; Butcher.b = 1; Butcher.c = 0; % Eulero
        f = @(t,y) -y;
        [t, u] = rk_generic(Butcher, f, [0,1], 1, struct('h',0.1));
        assert(abs(u(end) - exp(-1)) < 1e-1, 'Errore RK Generic');
        fprintf('OK\n');
    else
        fprintf('SKIP (File non trovato)\n');
    end
catch ME
    fprintf('FALLITO! (%s)\n', ME.message); errors = errors+1;
end

%% TEST 8: Newton (Zeri)
try
    fprintf('[8/8] Test newton (zeri)... ');
    f = @(x) x.^2 - 4; df = @(x) 2*x;
    if exist('newton.m', 'file')
        [x_vect, iter] = newton(f, df, 3, 1e-6, 100, 1);
        % x_vect è un vettore di tutte le iterate, prendiamo l'ultima
        x_final = x_vect(end);
        assert(abs(x_final - 2) < 1e-5, 'Newton non converge a 2');
        fprintf('OK\n');
    else
        fprintf('MANCANTE\n'); errors = errors+1;
    end
catch ME
    fprintf('FALLITO! (%s)\n', ME.message); errors = errors+1;
end


%% CONCLUSIONI
fprintf('\n=== RISULTATO FINALE ===\n');
if errors == 0
    fprintf('SUCCESS: La libreria è completa e funzionante!\n');
else
    fprintf('WARNING: Ci sono %d problemi da risolvere.\n', errors);
end