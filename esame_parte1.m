%% ESAME_15_02_2024_PARTE1.m
% Soluzione automatica Appello Parte 1 del 15/02/2024
% Richiede la libreria MatlabLib nel path.

clear; clc; close all;
format long e;

%% === TEST - ESERCIZIO 1 (Numeri Floating Point) ===
fprintf('\n=== ESERCIZIO 1 ===\n');
% F(beta, t, L, U) = F(2, 8, -20, 20)
beta = 2; t = 8; L = -20; U = 20;

% x_min = beta^(L-1)
xmin = beta^(L-1);
fprintf('xmin = %.16e\n', xmin);

% x_max = beta^U * (1 - beta^(-t))
xmax = beta^U * (1 - beta^(-t));
fprintf('xmax = %.16e\n', xmax);


%% === TEST - ESERCIZIO 2 (Serie di Taylor) ===
fprintf('\n=== ESERCIZIO 2 ===\n');
x_val = 2;
target_err = 1e-3;
true_val = exp(x_val);

% Cerco N minimo
N = 0;
approx = 0;
while true
    approx = approx + x_val^N / factorial(N);
    err = abs(true_val - approx);
    if err < target_err
        break;
    end
    N = N + 1;
end
fprintf('N minimo = %d\n', N);


%% === TEST - ESERCIZIO 3 (Sistemi Lineari Triangolari) ===
fprintf('\n=== ESERCIZIO 3 ===\n');
% 10 sistemi lineari Ax = b_j, A 90x90 triangolare inferiore
n = 90;
num_systems = 10;

% Costo sostituzione in avanti per un sistema: n^2
costo_uno = n^2; 
costo_totale = num_systems * costo_uno;
fprintf('Numero operazioni totale = %d\n', costo_totale);


%% === TEST - ESERCIZIO 4 (Richardson Precondizionato) ===
fprintf('\n=== ESERCIZIO 4 ===\n');
A = [8 -2 0; -2 10 -2; 0 -2 11];
b = [1; 1; 1];
P = [6 0 0; 0 6 0; 0 0 9];

% Calcolo autovalori di P^-1 A
P_inv = inv(P);
M = P_inv * A;
eigs_M = eig(M);
lambda_min = min(eigs_M);
lambda_max = max(eigs_M);

% Parametro ottimale alpha_opt = 2 / (lambda_min + lambda_max)
% Attenzione: se gli autovalori sono reali e positivi.
% Verifica: A e P sono simmetriche definite positive?
% A: 8>0, det(A_2)=76>0, det(A)=836-44-32 = 760 > 0 (OK)
% P: diagonale positiva (OK)
% Autovalori M reali positivi
alpha_opt = 2 / (lambda_min + lambda_max);
fprintf('alpha_opt = %.4f\n', alpha_opt);

% Iterazione Richardson
x_curr = b; % x(0) = b
B_rich = eye(3) - alpha_opt * M;
g_rich = alpha_opt * P_inv * b;

for k = 1:6
    x_curr = B_rich * x_curr + g_rich;
end
fprintf('x(6) = \n'); disp(x_curr);


%% === TEST - ESERCIZIO 5 (Convergenza QR) ===
fprintf('\n=== ESERCIZIO 5 ===\n');
% Autovalori: 10, 8, 6, lambda4, 3
% lambda4 in (3, 6)
% Convergenza QR dipende dal rapporto |lambda_i / lambda_{i-1}|
% Vogliamo massimizzare la velocità, quindi minimizzare il massimo rapporto.
% I rapporti critici che coinvolgono lambda4 sono:
% r1 = lambda4 / 6  (essendo lambda4 < 6)
% r2 = 3 / lambda4  (essendo 3 < lambda4)
% Dobbiamo minimizzare max(r1, r2)
% r1 = r2 => lambda4 / 6 = 3 / lambda4 => lambda4^2 = 18 => lambda4 = sqrt(18) = 3*sqrt(2)
lam4_opt = sqrt(18);
fprintf('lambda4 ottimale = %.4f\n', lam4_opt);


%% === TEST - ESERCIZIO 6 (QR Iterativo) ===
fprintf('\n=== ESERCIZIO 6 ===\n');
A6 = [5 0 0; 3 2 0; 4 1 1.9999999];
tol = 1e-6;
% Metodo QR
Ak = A6;
iter = 0;
while true
    iter = iter + 1;
    [Q, R] = qr(Ak);
    Ak = R * Q;
    
    % Criterio d'arresto: elemento sottodiagonale a(n, n-1)
    % Qui n=3. Elemento (3,2).
    % Nota: la matrice tende a triangolare superiore.
    % Gli autovalori sono sulla diagonale.
    % La convergenza è guidata dal rapporto lambda3/lambda2.
    % lambda1=5, lambda2=2, lambda3=1.9999999
    % lambda3 e lambda2 sono molto vicini! Convergenza lentissima.
    % Monitoriamo l'elemento (3,2)
    if abs(Ak(3,2)) < tol
        break;
    end
    
    if iter > 10000 % Safety break
        fprintf('Raggiunto limite iterazioni\n');
        break;
    end
end
fprintf('Iterazioni QR = %d\n', iter);
fprintf('Approssimazione lambda3 (A(3,3)) = %.16f\n', Ak(3,3));


%% === TEST - ESERCIZIO 7 (Newton Molteplicità) ===
fprintf('\n=== ESERCIZIO 7 ===\n');
% f(x) = (x-4) log(x-3) sin(pi*x)
% Zero alpha = 4
% Analisi ordine:
% f(x) = (x-4) * log(x-3) * sin(pi*x)
% A x=4:
% (x-4) -> ordine 1
% log(x-3) -> log(1)=0. Derivata 1/(x-3)->1. Ordine 1.
% sin(pi*x) -> sin(4pi)=0. Derivata pi*cos(4pi)=pi. Ordine 1.
% Totale molteplicità = 1 + 1 + 1 = 3.
% Verifica numerica o simbolica.
% f'(x) = log(x-3)sin(pi*x) + (x-4)/(x-3)sin(pi*x) + (x-4)log(x-3)pi*cos(pi*x)
% f'(4) = 0 + 0 + 0 = 0
% f''(x) ... sarà 0
% f'''(x) ... diverso da 0.
m = 3;
fprintf('Molteplicità m = %d\n', m);
% Newton classico su radice multipla -> ordine 1 (lineare)
fprintf('Ordine convergenza Newton classico = 1\n');


%% === TEST - ESERCIZIO 8 (Metodo Corde) ===
fprintf('\n=== ESERCIZIO 8 ===\n');
f8 = @(x) (x-4).*log(x-3).*sin(pi*x);
a = 3.5; b = 4.5;
x0 = 3.5;
% Metodo delle corde: x(k+1) = x(k) - f(x(k)) / q_c
% q_c = (f(b) - f(a)) / (b - a)
qc = (f8(b) - f8(a)) / (b - a);

x1 = x0 - f8(x0) / qc; % x(1)
x2 = x1 - f8(x1) / qc; % x(2)
x3 = x2 - f8(x2) / qc; % x(3)

fprintf('x(2) = %.16f\n', x2);
fprintf('x(3) = %.16f\n', x3);


%% === TEST - ESERCIZIO 9 (Punto Fisso Parametrico) ===
fprintf('\n=== ESERCIZIO 9 ===\n');
% phi(x) = 1/x * exp(gamma*(x-1))
% alpha = 1 è punto fisso: phi(1) = 1*1 = 1. OK.
% Convergenza locale: |phi'(1)| < 1
% phi'(x) = -1/x^2 * exp(...) + 1/x * exp(...) * gamma
% phi'(1) = -1 * 1 + 1 * 1 * gamma = gamma - 1
% Condizione: |gamma - 1| < 1
% -1 < gamma - 1 < 1
% 0 < gamma < 2
fprintf('Intervallo convergenza: gamma in (0, 2)\n');


%% === ESERCIZIO GRANDE (Punti 4-7) ===
fprintf('\n=== ESERCIZIO GRANDE (17pt) ===\n');
n = 100;
gamma_val = 10;

% Costruzione matrice A tridiagonale
% Diagonale principale: gamma
% Diagonali sup/inf: -1
e = ones(n, 1);
A_mat = spdiags([-e, gamma_val*e, -e], -1:1, n, n);
% Nota: spdiags mette -1 su tutta la colonna, anche dove non serve, ma l'intersezione con la matrice quadrata è corretta.
% Verifichiamo A(1,2)=-1, A(2,1)=-1.
% A(1,1)=gamma.

b_vec = ones(n, 1); % x_esatta = 1 -> b = A*1.
% A*ones: riga interna (2...n-1): -1 + gamma -1 = gamma - 2
% riga 1: gamma - 1
% riga n: -1 + gamma
% Il testo dice "La soluzione esatta ... è x = 1".
% Quindi b deve essere calcolato come A*x_esatta.
% Ma il testo al Punto 4 dice "considerando x(0) = b".
% Se b non è dato esplicitamente ma solo "sistema lineare Ax=b... soluzione esatta x=1", allora b = A*ones(n,1).
b_vec = A_mat * ones(n, 1);

% --- Punto 4: Gradiente ---
% Errore in norma A dopo k=10000 iterazioni.
% x(0) = b_vec.
% Metodo gradiente:
% r(k) = b - A x(k)
% alpha_k = (r' r) / (r' A r)
% x(k+1) = x(k) + alpha_k r(k)
% Stima errore senza eseguire?
% ||e(k)||_A <= ((K-1)/(K+1))^k ||e(0)||_A
% K = cond(A) spettrale.
% Autovalori di A (tridiagonale Toeplitz-like):
% lambda_j = gamma - 2*cos(j*pi/(n+1))
j = 1:n;
lam = gamma_val - 2*cos(j*pi/(n+1));
lam_min = min(lam);
lam_max = max(lam);
K = lam_max / lam_min;
fprintf('K(A) = %.4f\n', K);

factor = (K-1)/(K+1);
k_iter = 10000;
reduction = factor^k_iter;

% Errore iniziale ||e(0)||_A
x0 = b_vec;
x_ex = ones(n, 1);
e0 = x_ex - x0;
norm_e0_A = sqrt(e0' * A_mat * e0);

err_est = reduction * norm_e0_A;
fprintf('Stima errore ||e(10000)||_A = %.4f\n', err_est);


% --- Punto 5: Gradiente Precondizionato ---
% P = tridiag(-1, beta, -1)
% beta in {2, 4, 6, 8, 10}
% Convergenza più rapida -> K(P^-1 A) minimo.
% A ha diagonale 10. P ha diagonale beta.
% P è "simile" ad A se beta è vicino a 10.
% Se beta = 10, P = A => convergenza in 1 iterazione (K=1).
% Quindi beta = 10 dovrebbe essere la scelta migliore.
% Verifichiamo numericamente.
betas = [2, 4, 6, 8, 10];
min_cond = inf;
best_beta = 0;

for beta_i = betas
    P_mat = spdiags([-e, beta_i*e, -e], -1:1, n, n);
    % Calcolo cond(P^-1 A) generalizzato (autovalori di P^-1 A)
    % Poiché sono simmetriche definite positive, sono autovalori di P^-1 A.
    % Problema generalizzato: A v = lambda P v
    lams_gen = eig(full(A_mat), full(P_mat));
    cond_gen = max(lams_gen) / min(lams_gen);
    fprintf('beta=%d, cond=%.4f\n', beta_i, cond_gen);
    
    if cond_gen < min_cond
        min_cond = cond_gen;
        best_beta = beta_i;
    end
end
fprintf('Miglior beta = %d\n', best_beta);


% --- Punto 6: Stima K(A^-1) ---
% K(A^-1) = K(A) (per matrice simm def pos).
% Algoritmo efficiente: metodo delle potenze per lambda_max e potenze inverse per lambda_min.
% Oppure potenze su A per lam_max e potenze su A^-1 (inverse) per lam_min.
% Qui si chiede K(A^-1). Gli autovalori di A^-1 sono 1/lam_i(A).
% Max autovalore A^-1 = 1 / min(lam(A))
% Min autovalore A^-1 = 1 / max(lam(A))
% K(A^-1) = (1/lam_min) / (1/lam_max) = lam_max / lam_min = K(A).
% L'algoritmo deve stimare questo.
% Algoritmo proposto: Metodo delle potenze su A (per lam_max) e Potenze Inverse su A (per lam_min).
% Implementiamo Potenze e Potenze Inverse.

% Potenze per lam_max(A)
z = ones(n, 1);
for iter=1:50
    w_vec = A_mat * z;
    z = w_vec / norm(w_vec);
    lam_max_est = z' * A_mat * z;
end

% Potenze inverse per lam_min(A)
z = ones(n, 1);
% Fattorizzazione per risolvere efficientemente
[L_chol, U_chol] = lu(A_mat); % O chol se simmetrica
for iter=1:50
    % w = A \ z
    w_vec = U_chol \ (L_chol \ z);
    z = w_vec / norm(w_vec);
    % lam_min_est = 1 / (z' * A^-1 z) ? No, Rayleigh su A
    lam_min_est = z' * A_mat * z;
end

K_est = lam_max_est / lam_min_est;
fprintf('Stima K(A^-1) = %.4f\n', K_est);


% --- Punto 7: Metodo Quasi-Newton (Broyden) ---
% Algoritmo 1 è il metodo di Broyden (rango 1).
% F(x) = Ax + 1 - exp(-x/50)
% x0 = 0.1 * ones
% B0 = A
fprintf('\n--- Punto 7 (Broyden) ---\n');

% Funzione F
F_fun = @(x) A_mat * x + 1 - exp(-x/50);

x_curr = 0.1 * ones(n, 1);
B_curr = A_mat;

% Iterazioni
for k = 1:3
    F_k = F_fun(x_curr);
    
    % Risolvo Bk * delta = -Fk
    delta = B_curr \ (-F_k);
    
    x_next = x_curr + delta;
    
    % Aggiornamento Broyden
    % z = F(x_next) - F(x_curr) - B_curr * delta 
    % Ma B_curr * delta = -F_k
    % Quindi z = F(x_next) - F_k - (-F_k) = F(x_next)
    % La formula nel testo è z = F(x_next) - F(x_curr) - B_curr * delta
    % Se facciamo il calcolo esplicito:
    F_next = F_fun(x_next);
    z_vec = F_next - F_k - B_curr * delta; % = F_next
    
    % B_next = B_curr + (z * delta') / (delta' * delta)
    % Attenzione: questa è una correzione di rango 1 densa. B diventa piena.
    % Per n=100 va bene.
    B_next = B_curr + (z_vec * delta') / (delta' * delta);
    
    fprintf('Iter %d: (x(%d))_1 = %.16f\n', k, k, x_next(1));
    
    x_curr = x_next;
    B_curr = B_next;
end

fprintf('\n=== ESECUZIONE COMPLETATA ===\n');
