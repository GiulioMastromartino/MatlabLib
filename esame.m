%% ESAME_15_02_2024.m
% Soluzione automatica Appello Parte 2 del 15/02/2024
% Richiede la libreria MatlabLib nel path.

clear; clc; close all;
format long e;

%% === TEST - ESERCIZIO 1 (Interpolazione) ===
fprintf('\n=== ESERCIZIO 1 ===\n');
xi = [0, 1, 2, 3];
yi = [2.5, 0.5, 1.5, 0.5];
x_eval = 0.75;

% Costruisco il polinomio interpolante (grado 3 per 4 punti)
c = polyfit(xi, yi, 3); 

% Valutazione in 0.75
val = polyval(c, x_eval);
fprintf('Pi_3(0.75) = %.4f\n', val);

% Derivata: se p(x) = a*x^3 + b*x^2 + c*x + d
% p'(x) = 3a*x^2 + 2b*x + c
der_c = [3*c(1), 2*c(2), c(3)];
val_der = polyval(der_c, x_eval);
fprintf('dPi_3/dx(0.75) = %.4f\n', val_der);


%% === TEST - ESERCIZIO 2 (Integrale Polinomio Lineare a Tratti) ===
fprintf('\n=== ESERCIZIO 2 ===\n');
xi_2 = [0, 0.25, 0.5, 0.75, 1];
yi_2 = [2, 0.5, -0.5, 0.5, 3];

% L'integrale del polinomio lineare a tratti (interpolante) è esattamente
% l'integrale calcolato con la formula dei trapezi sui dati.
I_H1 = trapz(xi_2, yi_2);
fprintf('Integrale Pi_1^H(x) = %.4f\n', I_H1);


%% === TEST - ESERCIZIO 3 (Minimi Quadrati) ===
fprintf('\n=== ESERCIZIO 3 ===\n');
xj = [0, 0.5, 1, 1.5, 2];
yj = [4, 0.25, 0.5, -1.5, 4];

% a) Grado 3, tutti i 5 punti
p3_a = polyfit(xj, yj, 3);
vals_a = polyval(p3_a, xj);
scarto_a = sum((vals_a - yj).^2);
fprintf('Scarto quadratico (n+1=5) = %.4f\n', scarto_a);

% b) Grado 3, solo primi 4 punti (Interpolazione!)
xj_b = xj(1:4);
yj_b = yj(1:4);
p3_b = polyfit(xj_b, yj_b, 3);
vals_b = polyval(p3_b, xj_b);
scarto_b = sum((vals_b - yj_b).^2); % Dovrebbe essere 0 (o quasi macchina)
fprintf('Scarto quadratico (n+1=4) = %.4f\n', scarto_b);


%% === TEST - ESERCIZIO 4 (Punto Medio Composito) ===
fprintf('\n=== ESERCIZIO 4 ===\n');
f4 = @(x) sqrt(3 + 2*abs(x));
a = -1; b = 1; H = 0.5;
% Sottointervalli: [-1, -0.5], [-0.5, 0], [0, 0.5], [0.5, 1]
x_nodes = a:H:b;
mid_points = x_nodes(1:end-1) + H/2;
I_pm = H * sum(f4(mid_points));
fprintf('Integrale Punto Medio = %.4f\n', I_pm);


%% === TEST - ESERCIZIO 5 (Gauss-Legendre) ===
fprintf('\n=== ESERCIZIO 5 ===\n');
f5 = @(x) sin(x).*cos(x);
a5 = 0; b5 = 2;
% Gauss-Legendre ordine 3 su [-1, 1]
nodes_ref = [-sqrt(3/5), 0, sqrt(3/5)];
weights_ref = [5/9, 8/9, 5/9];

% Mappa su [0, 2]
% x = (b-a)/2 * t + (b+a)/2 = t + 1
% dx = dt
nodes_phy = nodes_ref + 1;
weights_phy = weights_ref; % * (b-a)/2 = * 1

I_GL = sum(weights_phy .* f5(nodes_phy));
fprintf('Integrale Gauss-Legendre = %.4f\n', I_GL);
fprintf('Grado di esattezza = %d (2*3 - 1)\n', 2*3 - 1);


%% === TEST - ESERCIZIO 6 (Cauchy Eulero) ===
fprintf('\n=== ESERCIZIO 6 ===\n');
f6 = @(t,y) -t*y^2 - 9*t;
y0 = 1;
h = 0.1;
% u0 = 1 (t=0)
% u1 = u0 + h*f(t0, u0)
u0 = 1; t0 = 0;
u1 = u0 + h * f6(t0, u0);
% u2 = u1 + h*f(t1, u1)
t1 = 0.1;
u2 = u1 + h * f6(t1, u1);
fprintf('u2 (Eulero avanti) = %.4f\n', u2);


%% === TEST - ESERCIZIO 8 (BVP Misto) ===
fprintf('\n=== ESERCIZIO 8 ===\n');
% Problema: -2u'' + 8u' + 4(4-x)u = 0
% u(0)=3, u'(4)=0
% Metodo: diff finite centrate + backward diff ordine 1 su bordo dx
% ATTENZIONE: bvp_mini_solver gestisce solo Dirichlet. 
% Dobbiamo implementare lo schema ad hoc.
L = 4; N = 40; h = L/N;
x = linspace(0, L, N+1)';
% Equazione interna nodi 1...N-1 (indici Matlab 2...N)
% -2 (u_{i+1}-2u_i+u_{i-1})/h^2 + 8 (u_{i+1}-u_{i-1})/(2h) + 4(4-xi)u_i = 0

% Nodo N (x=4): u'(4) = 0 approssimato con differenze finite backward ord 1
% (u_N - u_{N-1})/h = 0 => u_N = u_{N-1}

A = zeros(N+1, N+1);
b = zeros(N+1, 1);

% Nodo 0 (Dirichlet)
A(1,1) = 1; b(1) = 3;

for i = 2:N
    xi_val = x(i);
    % u_{i-1} coeff: -2/h^2 - 8/(2h) = -2/h^2 - 4/h
    c_m1 = - (-2/h^2) - 8/(2*h); % Segno eq: moltiplico tutto per -1 per avere diag positiva?
    % No, scriviamo i coeff diretti:
    % (-2/h^2)*u_{i+1} - (-2)*(-2/h^2)*u_i + (-2/h^2)*u_{i-1} ...
    
    % Riscrivo: -2/h^2 (u_{i+1} - 2u_i + u_{i-1}) + 4/h (u_{i+1} - u_{i-1}) + reaction*u_i = 0
    % u_{i-1}: -2/h^2 - 4/h
    % u_i:     +4/h^2 + 4(4-xi)
    % u_{i+1}: -2/h^2 + 4/h
    
    A(i, i-1) = -2/h^2 - 4/h;
    A(i, i)   =  4/h^2 + 4*(4-xi_val);
    A(i, i+1) = -2/h^2 + 4/h;
end

% Nodo N+1 (Neumann backward 1: u_{N+1} - u_N = 0)
% Indice Matlab N+1 corrisponde a x_N
A(N+1, N+1) = 1;
A(N+1, N)   = -1;
b(N+1) = 0;

u = A\b;
fprintf('u(4) approx = %.4f\n', u(end));


%% === TEST - ESERCIZIO 9 (BVP Upwind) ===
fprintf('\n=== ESERCIZIO 9 ===\n');
% -2u'' + 400u' = 5
% u(0)=0, u(1)=4
% N=10 => h=0.1
% Usa bvp_mini_solver con 'upwind'
coeffs = [400, 0, 5]; % [p, q, f] per eq: -u'' + p u' + q u = f
bc = [0, 4];
N = 10;

% Usa la libreria!
u_up = bvp_mini_solver(N, coeffs, bc, 'upwind');

% Approssimazione u'(1) con diff finite backward
% (u_N - u_{N-1}) / h
h = 1/N;
der_1 = (u_up(end) - u_up(end-1)) / h;
fprintf('u''(1) approx = %.4f\n', der_1);


%% === TEST - ESERCIZIO 10 (Calore) ===
fprintf('\n=== ESERCIZIO 10 ===\n');
% du/dt - 2 d2u/dx2 = 0
% mu = 2
% f = 0
% Dirichlet omogenee
% u0 = 20 sin(pi x)
% h = 0.1, dt = 1e-3
% Eulero avanti => theta = 0

mu = 2;
f_fun = @(x,t) 0*x;
a = 0; b = 1;
us = @(t) 0; ud = @(t) 0;
g0 = @(x) 20*sin(pi*x);
T = 0.003; % Per arrivare a 3 passi (t=0, 0.001, 0.002, 0.003)
h = 0.1;
dt = 1e-3;
theta = 0; % Eulero esplicito

[u_mat, x_vec, t_vec] = EqCalore_DiffFin_Theta(mu, f_fun, a, b, us, ud, g0, T, h, dt, theta);

% u(0.1, 0.003)
% x=0.1 è indice 2 (x=0, 0.1, ...)
% t=0.003 è indice 4 (t=0, 0.001, 0.002, 0.003)
val_calore = u_mat(2, 4);
fprintf('u(0.1, 0.003) = %.4f\n', val_calore);


%% === ESERCIZIO GRANDE (Sistema ODE) ===
fprintf('\n=== ESERCIZIO GRANDE (17pt) ===\n');

% Parametri
A_mat = [-2/3, -4, 0, 4/3; 
          1,    0, 0, 0; 
          0,    4, -3, -4; 
          0,    0, 1, 0];
      
f1 = @(t) -(6349/25)*cos(4*t)*exp(-t/5) - 36*sin(3*t)*exp(-t/5) - (112/5)*sin(4*t)*exp(-t/5);
f2 = @(t) -(28*cos(4*t) - (351/5)*cos(3*t) + (1251/25)*sin(3*t)) * exp(-t/5);

% g(t) = (f1/3, 0, f2, 0)'
g_fun = @(t) [f1(t)/3; 0; f2(t); 0];

y0 = [-7/5; 7; 27; 0];
tf = 50;
h = 0.1;

% Soluzione esatta
sol_exact = @(t) [7*cos(4*t).*exp(-t/5); 9*sin(3*t).*exp(-t/5)]; % x1, x2
% Nota: il vettore y ha 4 comp: [w1, x1, w2, x2]. La sol esatta data è solo per x1, x2.
% Quindi x1 = y(2), x2 = y(4).

% --- Punto 2: Heun ---
% Definisco la funzione per Heun
F_sys = @(t,y) A_mat * y + g_fun(t);

[t_h, u_h] = Heun(F_sys, tf, y0, h);

% Valori u1 e uNh
u1_vec = u_h(:, 2); % t1 = h
uNh_vec = u_h(:, end); % tf

fprintf('u1 (Heun) = \n'); disp(u1_vec);
fprintf('uNh (Heun) = \n'); disp(uNh_vec);

% Tempo tm in cui norma < 0.5
% (un)2 è x1, (un)4 è x2
x1_num = u_h(2,:);
x2_num = u_h(4,:);

mask = (abs(x1_num) < 0.5) & (abs(x2_num) < 0.5);
% Trova il primo indice tale che da lì in poi è sempre vero
% Scorriamo all'indietro
last_false = find(~mask, 1, 'last');
if isempty(last_false)
    tm_idx = 1;
else
    tm_idx = last_false + 1;
end
tm = t_h(tm_idx);
fprintf('tm = %.1f\n', tm);

% --- Punto 3: Errori Heun ---
hs = [1e-1, 5e-2, 2.5e-2, 1.25e-2];
errs = zeros(length(hs), 1);

for i = 1:length(hs)
    hi = hs(i);
    [ti, ui] = Heun(F_sys, tf, y0, hi);
    
    % Soluzione esatta su griglia
    x1_ex = 7 * cos(4*ti) .* exp(-ti/5);
    x1_num = ui(2,:);
    
    errs(i) = max(abs(x1_num - x1_ex));
end
fprintf('Errori Heun: \n'); disp(errs);

% --- Punto 4: Ordine ---
ratio = errs(end-1) / errs(end);
order = log2(ratio);
fprintf('Ordine stimato p = %.4f\n', order);

% --- Punto 5: Stabilità Assoluta ---
% Heun assoluta stabilità: |1 + z + z^2/2| < 1
% z = h * lambda
lam = eig(A_mat);
re_lam = real(lam);
% Poiché gli autovalori hanno parte reale negativa (probabilmente, da verificare),
% cerchiamo h limite.
% Per Heun e autovalori complessi/immaginari, la regione tocca l'asse imm?
% Se sono reali negativi, intervallo (-2, 0).
fprintf('Autovalori A:\n'); disp(lam);

% --- Punto 6: Runge-Kutta ---
% Tableau
Butcher.A = [0 0 0; 1 0 0; 1/4 1/4 0]; % Matrice A
Butcher.b = [1/6 1/6 2/3];             % Vettore b
Butcher.c = [0; 1; 1/2];               % Vettore c (trasposto in colonna)
% La tabella data:
% 0 | 0 0 0
% 1 | 1 0 0
% 1/2 | 1/4 1/4 0
% ----------------
%     | 1/6 1/6 2/3
% Metodo esplicito.

[t_rk, u_rk] = rk_generic(Butcher, F_sys, [0, tf], y0, struct('h', 0.1));

u1_rk = u_rk(:, 2);
uNh_rk = u_rk(:, end);

fprintf('u1 (RK) = \n'); disp(u1_rk);
fprintf('uNh (RK) = \n'); disp(uNh_rk);

%% Finito
fprintf('\n=== ESECUZIONE COMPLETATA ===\n');
