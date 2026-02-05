%% ESAME_15_02_2024_PARTE1_FIX.m
% Soluzione corretta Appello Parte 1 del 15/02/2024
% Correzione parametro gamma=2 (dedotto dalle soluzioni)

clear; clc; close all;
format long e;

%% === TEST - ESERCIZIO 1 (Numeri Floating Point) ===
fprintf('\n=== ESERCIZIO 1 ===\n');
beta = 2; t = 8; L = -20; U = 20;
xmin = beta^(L-1);
xmax = beta^U * (1 - beta^(-t));
fprintf('xmin = %.16e\n', xmin);
fprintf('xmax = %.16e\n', xmax);

%% === TEST - ESERCIZIO 2 (Serie di Taylor) ===
fprintf('\n=== ESERCIZIO 2 ===\n');
x_val = 2; target_err = 1e-3; true_val = exp(x_val);
N = 0; approx = 0;
while true
    approx = approx + x_val^N / factorial(N);
    if abs(true_val - approx) < target_err, break; end
    N = N + 1;
end
fprintf('N minimo = %d\n', N);

%% === TEST - ESERCIZIO 3 (Sistemi) ===
fprintf('\n=== ESERCIZIO 3 ===\n');
fprintf('Numero operazioni totale = %d\n', 10 * 90^2);

%% === TEST - ESERCIZIO 4 (Richardson) ===
fprintf('\n=== ESERCIZIO 4 ===\n');
A = [8 -2 0; -2 10 -2; 0 -2 11];
b = [1; 1; 1];
P = [6 0 0; 0 6 0; 0 0 9];
M = P \ A;
eigs_M = eig(M);
alpha_opt = 2 / (min(eigs_M) + max(eigs_M));
fprintf('alpha_opt = %.4f\n', alpha_opt);
x_curr = b;
B_rich = eye(3) - alpha_opt * M;
g_rich = alpha_opt * (P \ b);
for k = 1:6, x_curr = B_rich * x_curr + g_rich; end
fprintf('x(6) = \n'); disp(x_curr);

%% === TEST - ESERCIZIO 5 (QR) ===
fprintf('\n=== ESERCIZIO 5 ===\n');
fprintf('lambda4 ottimale = %.4f\n', sqrt(18));

%% === TEST - ESERCIZIO 6 (QR Iter) ===
fprintf('\n=== ESERCIZIO 6 ===\n');
A6 = [5 0 0; 3 2 0; 4 1 1.9999999];
Ak = A6; iter = 0;
while abs(Ak(3,2)) >= 1e-6 && iter < 10000
    [Q, R] = qr(Ak); Ak = R * Q; iter = iter + 1;
end
fprintf('Iterazioni QR = %d\n', iter);
fprintf('Approssimazione lambda3 = %.16f\n', Ak(3,3));

%% === TEST - ESERCIZIO 7 (Newton) ===
fprintf('\n=== ESERCIZIO 7 ===\n');
fprintf('Molteplicità m = 3, Ordine = 1\n');

%% === TEST - ESERCIZIO 8 (Corde) ===
fprintf('\n=== ESERCIZIO 8 ===\n');
f8 = @(x) (x-4).*log(x-3).*sin(pi*x);
x0 = 3.5; qc = (f8(4.5) - f8(3.5)) / (4.5 - 3.5);
x1 = x0 - f8(x0)/qc;
x2 = x1 - f8(x1)/qc;
x3 = x2 - f8(x2)/qc;
fprintf('x(2) = %.16f\n', x2);
fprintf('x(3) = %.16f\n', x3);

%% === TEST - ESERCIZIO 9 (Punto Fisso) ===
fprintf('\n=== ESERCIZIO 9 ===\n');
fprintf('gamma in (0, 2)\n');

%% === ESERCIZIO GRANDE (Con gamma=2) ===
fprintf('\n=== ESERCIZIO GRANDE (gamma=2) ===\n');
n = 100;
gamma_val = 2; % CORREZIONE FONDAMENTALE
e = ones(n, 1);
A_mat = spdiags([-e, gamma_val*e, -e], -1:1, n, n);

% --- Punto 4 ---
% Stima errore ||e(10000)||_A
% x_esatta = 1. b = A*1. x0 = b.
% e0 = x_esatta - x0 = 1 - A*1.
x_ex = ones(n, 1);
b_vec = A_mat * x_ex;
x0 = b_vec;
e0 = x_ex - x0;
norm_e0_A = sqrt(e0' * A_mat * e0);

% Condizionamento A
lam = gamma_val - 2*cos((1:n)'*pi/(n+1));
lam_min = min(lam); lam_max = max(lam);
K = lam_max / lam_min;
rho = (K-1)/(K+1);
err_est = rho^10000 * norm_e0_A;
fprintf('Stima errore Punto 4 = %.4f\n', err_est);

% --- Punto 5 ---
% Beta ottimale. Se gamma=2, A ha diagonale 2.
% P ha diagonale beta. Se beta=2 -> P=A -> K=1 -> convergenza immediata.
fprintf('Miglior beta = 2 (rende P=A)\n');

% --- Punto 6 ---
% K(A^-1) = K(A)
fprintf('Stima K(A^-1) = %.4f\n', K);

% --- Punto 7 (Broyden) ---
fprintf('\n--- Punto 7 (Broyden) ---\n');
F_fun = @(x) A_mat * x + 1 - exp(-x/50);
x_curr = 0.1 * ones(n, 1);
B_curr = A_mat; % B0 = A

% Memorizziamo le prime componenti
x1_hist = [];

for k = 1:3
    F_k = F_fun(x_curr);
    
    % Risolvo sistema lineare. 
    % Nota: B_curr diventa densa (rank-1 updates).
    % Per n=100 \ è ok.
    delta = B_curr \ (-F_k);
    
    x_next = x_curr + delta;
    x1_hist(k) = x_next(1);
    
    % Update Broyden
    F_next = F_fun(x_next);
    y_k = F_next - F_k; % In algo 1 z(k+1) = F(x+) - F(x) - Bk d. Ma Bk d = -F(x). Quindi z = F(x+)
    % Verifica Algoritmo 1 nel testo:
    % z(k+1) = F(x_new) - F(x_old) - Bk * delta
    % Ma Bk*delta = -F(x_old)
    % Quindi z(k+1) = F(x_new)
    z_vec = F_next; 
    
    B_next = B_curr + (z_vec * delta') / (delta' * delta);
    
    x_curr = x_next;
    B_curr = B_next;
end

fprintf('x(1)_1 = %.4f\n', x1_hist(1));
fprintf('x(2)_1 = %.4f\n', x1_hist(2));
fprintf('x(3)_1 = %.4f\n', x1_hist(3));

fprintf('\n=== ESECUZIONE COMPLETATA ===\n');
