%% ESAME_15_02_2024_FIXED.m
% Soluzione corretta Appello Parte 2 del 15/02/2024
% Include fix per Esercizio 5 (Gauss 4 nodi) e Esercizio 9 (BVP derivata)

clear; clc; close all;
format long e;

%% === TEST - ESERCIZIO 1 (Interpolazione) ===
fprintf('\n=== ESERCIZIO 1 ===\n');
xi = [0, 1, 2, 3];
yi = [2.5, 0.5, 1.5, 0.5];
x_eval = 0.75;

c = polyfit(xi, yi, 3); 
val = polyval(c, x_eval);
fprintf('Pi_3(0.75) = %.4f\n', val);

der_c = [3*c(1), 2*c(2), c(3)];
val_der = polyval(der_c, x_eval);
fprintf('dPi_3/dx(0.75) = %.4f\n', val_der);


%% === TEST - ESERCIZIO 2 (Integrale Polinomio) ===
fprintf('\n=== ESERCIZIO 2 ===\n');
xi_2 = [0, 0.25, 0.5, 0.75, 1];
yi_2 = [2, 0.5, -0.5, 0.5, 3];
I_H1 = trapz(xi_2, yi_2);
fprintf('Integrale Pi_1^H(x) = %.4f\n', I_H1);


%% === TEST - ESERCIZIO 3 (Minimi Quadrati) ===
fprintf('\n=== ESERCIZIO 3 ===\n');
xj = [0, 0.5, 1, 1.5, 2];
yj = [4, 0.25, 0.5, -1.5, 4];

p3_a = polyfit(xj, yj, 3);
vals_a = polyval(p3_a, xj);
scarto_a = sum((vals_a - yj).^2);
fprintf('Scarto quadratico (n+1=5) = %.4f\n', scarto_a);

xj_b = xj(1:4);
yj_b = yj(1:4);
p3_b = polyfit(xj_b, yj_b, 3);
vals_b = polyval(p3_b, xj_b);
scarto_b = sum((vals_b - yj_b).^2);
fprintf('Scarto quadratico (n+1=4) = %.4f\n', scarto_b);


%% === TEST - ESERCIZIO 4 (Punto Medio) ===
fprintf('\n=== ESERCIZIO 4 ===\n');
f4 = @(x) sqrt(3 + 2*abs(x));
a = -1; b = 1; H = 0.5;
x_nodes = a:H:b;
mid_points = x_nodes(1:end-1) + H/2;
I_pm = H * sum(f4(mid_points));
fprintf('Integrale Punto Medio = %.4f\n', I_pm);


%% === TEST - ESERCIZIO 5 (Gauss-Legendre) ===
fprintf('\n=== ESERCIZIO 5 ===\n');
% Nota: Il testo chiede "ordine 3" ma la soluzione implica grado di esattezza 7.
% Grado 7 richiede 4 nodi (2n-1=7 => 2n=8 => n=4 punti).
% "Ordine 3" nel testo potrebbe riferirsi all'indice n=3 (0,1,2,3 -> 4 punti).
fprintf('Utilizzo 4 nodi (grado 7) per matchare la soluzione attesa.\n');

f5 = @(x) sin(x).*cos(x);
% Nodi e pesi per n=4 su [-1, 1]
nodes_ref = [-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116];
weights_ref = [0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451];

% Mappa su [0, 2]
% x = (b-a)/2 * t + (b+a)/2 = t + 1
nodes_phy = nodes_ref + 1;
weights_phy = weights_ref; 

I_GL = sum(weights_phy .* f5(nodes_phy));
fprintf('Integrale Gauss-Legendre = %.4f\n', I_GL);
fprintf('Grado di esattezza = 7\n');


%% === TEST - ESERCIZIO 6 (Cauchy Eulero) ===
fprintf('\n=== ESERCIZIO 6 ===\n');
f6 = @(t,y) -t*y^2 - 9*t;
y0 = 1; h = 0.1;
u0 = 1; t0 = 0;
u1 = u0 + h * f6(t0, u0);
t1 = 0.1;
u2 = u1 + h * f6(t1, u1);
fprintf('u2 (Eulero avanti) = %.4f\n', u2);


%% === TEST - ESERCIZIO 8 (BVP Misto) ===
fprintf('\n=== ESERCIZIO 8 ===\n');
% -2u'' + 8u' + 4(4-x)u = 0, u(0)=3, u'(4)=0
% BDF1 per u'(4)
L = 4; N = 40; h = L/N;
x = linspace(0, L, N+1)';
A = zeros(N+1, N+1);
b = zeros(N+1, 1);

A(1,1) = 1; b(1) = 3;
for i = 2:N
    xi_val = x(i);
    A(i, i-1) = -2/h^2 - 4/h;
    A(i, i)   =  4/h^2 + 4*(4-xi_val);
    A(i, i+1) = -2/h^2 + 4/h;
end
A(N+1, N+1) = 1;
A(N+1, N)   = -1;
b(N+1) = 0;

u = A\b;
fprintf('u(4) approx = %.4f\n', u(end));


%% === TEST - ESERCIZIO 9 (BVP Upwind) ===
fprintf('\n=== ESERCIZIO 9 ===\n');
% -2u'' + 400u' = 5, u(0)=0, u(1)=4
% Upwind manuale per sicurezza
mu = 2; eta = 400; f_val = 5;
L = 1; N = 10; h = L/N;
x = linspace(0, L, N+1)';

% Matrice per nodi interni 2...N (x_1...x_{N-1})
% Dimensione sistema interno N-1
% Nodi totali N+1 (0...N). Nodi interni da 2 a N.
A_up = zeros(N-1, N-1);
b_up = zeros(N-1, 1) + f_val;

% Coeff Upwind (eta > 0 -> backward su convettivo)
% -mu * D2_cent + eta * D1_back
% D2_cent u_i = (u_{i+1} - 2u_i + u_{i-1})/h^2
% D1_back u_i = (u_i - u_{i-1})/h
% Eq: -mu/h^2 * u_{i+1} + (2mu/h^2 + eta/h) * u_i + (-mu/h^2 - eta/h) * u_{i-1} = f

c_super = -mu/h^2;
c_diag  = 2*mu/h^2 + eta/h;
c_sub   = -mu/h^2 - eta/h;

% Costruzione matrice
for k = 1:N-1
    % Riga k corrisponde a nodo interno k+1 (x_k)
    % Ma in x indicizzati 1..N+1, interni sono 2..N
    % Diciamo u_int vettore u(2:end-1)
    
    A_up(k,k) = c_diag;
    if k > 1
        A_up(k, k-1) = c_sub;
    end
    if k < N-1
        A_up(k, k+1) = c_super;
    end
end

% Condizioni al bordo
u_0 = 0;
u_N = 4;

% Correzione termine noto con BC
b_up(1) = b_up(1) - c_sub * u_0;
b_up(end) = b_up(end) - c_super * u_N;

u_int = A_up \ b_up;
u_full = [u_0; u_int; u_N];

% Derivata u'(1) con BDF2 (ordine 2 backward)
% u'(x_N) ~ (3u_N - 4u_{N-1} + u_{N-2}) / 2h
der_1 = (3*u_full(end) - 4*u_full(end-1) + u_full(end-2)) / (2*h);
fprintf('u''(1) approx (BDF2) = %.4f\n', der_1);

% Controllo BDF1
der_1_bdf1 = (u_full(end) - u_full(end-1)) / h;
fprintf('u''(1) approx (BDF1) = %.4f\n', der_1_bdf1);


%% === TEST - ESERCIZIO 10 (Calore) ===
fprintf('\n=== ESERCIZIO 10 ===\n');
mu = 2;
a = 0; b = 1;
g0 = @(x) 20*sin(pi*x);
T = 0.003; 
h = 0.1; dt = 1e-3;
theta = 0; 

[u_mat, ~, ~] = EqCalore_DiffFin_Theta(mu, @(x,t) 0, a, b, @(t) 0, @(t) 0, g0, T, h, dt, theta);
fprintf('u(0.1, 0.003) = %.4f\n', u_mat(2, 4));


%% === ESERCIZIO GRANDE (Sistema ODE) ===
fprintf('\n=== ESERCIZIO GRANDE (17pt) ===\n');
A_mat = [-2/3, -4, 0, 4/3; 
          1,    0, 0, 0; 
          0,    4, -3, -4; 
          0,    0, 1, 0];
      
f1 = @(t) -(6349/25)*cos(4*t)*exp(-t/5) - 36*sin(3*t)*exp(-t/5) - (112/5)*sin(4*t)*exp(-t/5);
f2 = @(t) -(28*cos(4*t) - (351/5)*cos(3*t) + (1251/25)*sin(3*t)) * exp(-t/5);
g_fun = @(t) [f1(t)/3; 0; f2(t); 0];
y0 = [-7/5; 7; 27; 0];
tf = 50; h = 0.1;

F_sys = @(t,y) A_mat * y + g_fun(t);
[t_h, u_h] = Heun(F_sys, tf, y0, h);

fprintf('u1 (Heun) = \n'); disp(u_h(:, 2));
fprintf('uNh (Heun) = \n'); disp(u_h(:, end));

x1_num = u_h(2,:); x2_num = u_h(4,:);
mask = (abs(x1_num) < 0.5) & (abs(x2_num) < 0.5);
last_false = find(~mask, 1, 'last');
tm = t_h(last_false + 1);
fprintf('tm = %.1f\n', tm);

% Errori e ordine
hs = [1e-1, 5e-2, 2.5e-2, 1.25e-2];
errs = zeros(length(hs), 1);
for i = 1:length(hs)
    hi = hs(i);
    [ti, ui] = Heun(F_sys, tf, y0, hi);
    x1_ex = 7 * cos(4*ti) .* exp(-ti/5);
    errs(i) = max(abs(ui(2,:) - x1_ex));
end
fprintf('Errori Heun: \n'); disp(errs);
fprintf('Ordine stimato p = %.4f\n', log2(errs(end-1)/errs(end)));

% RK
Butcher.A = [0 0 0; 1 0 0; 1/4 1/4 0];
Butcher.b = [1/6 1/6 2/3];
Butcher.c = [0; 1; 1/2];
[~, u_rk] = rk_generic(Butcher, F_sys, [0, tf], y0, struct('h', 0.1));
fprintf('u1 (RK) = \n'); disp(u_rk(:, 2));
fprintf('uNh (RK) = \n'); disp(u_rk(:, end));

% Interpolazione
fprintf('\n=== Punto 7 (Interpolazione) ===\n');
% Stima errore interpolazione lineare a tratti
% e <= H^2/8 * max|f''(x)|
% Qui la funzione è x1(t)
% x1(t) = 7 cos(4t) e^(-t/5)
% Calcoliamo x1''(t) simbolicamente o numericamente fine
% x1'(t) = 7 * (-4 sin(4t) - 1/5 cos(4t)) e^(-t/5)
% x1''(t) = 7 * (-16 cos(4t) + 4/5 sin(4t) + 4/25 cos(4t) - 1/25 (-4 sin(4t) - 1/5 cos(4t))) ? No
% x1''(t) = 7 * e^(-t/5) * [ -16 cos - 1/5(-4 sin) - 4/5 sin + 1/25 cos ]
%         = 7 * e^(-t/5) * [ (-16 + 1/25) cos + (4/5 - 4/5) sin ]
%         = 7 * e^(-t/5) * [ -399/25 cos(4t) ]
% max |x1''| su [0, 50]. 
% Il termine esponenziale decresce. Il coseno oscilla.
% Massimo a t=0? |cos(0)|=1, e^0=1 -> |7 * -399/25| = 7 * 15.96 = 111.72
% Controlliamo se ci sono picchi. L'inviluppo è decrescente.
% Quindi max |x1''| = |x1''(0)| = 111.72
H = 0.1;
max_der2 = abs(7 * (-399/25)); 
err_interp = (H^2 / 8) * max_der2;
fprintf('Errore interpolazione stimato = %.4f\n', err_interp);

fprintf('\n=== ESECUZIONE COMPLETATA ===\n');
