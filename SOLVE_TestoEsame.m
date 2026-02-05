%% SOLVE_TestoEsame.m
% Risolutore automatico per gli esercizi in TestoEsame.pdf.
% Esegue calcoli numerici e stampa i risultati.
%
% Requisiti: aggiungere MatlabLib al path (si assume esecuzione dalla root del repo).
clear; clc;

fprintf('=== SOLVE_TestoEsame ===\n');

%% Helpers
doPlot = false;  % metti true per disegnare i cerchi di Gershgorin via gershcircles.m

%% ES 1) Floating point F(2,6,-10,10)
fprintf('\n[ES1] Floating point F(2,6,-10,10)\n');
beta = 2; t = 6; L = -10; U = 10;
x = pi^2;

if exist('floating_point_representation','file')
    [flx, abs_err, rel_err, u] = floating_point_representation(x, beta, t, L, U);
else
    error('floating_point_representation.m non trovato');
end

% secondo numero reale più grande rappresentabile esattamente:
% mantissa max = (beta^t - 1)/beta^t, second max = (beta^t - 2)/beta^t, exponent = U
m2 = (beta^t - 2) / beta^t;
second_largest = m2 * beta^U;

fprintf('fl(pi^2) = %.16g\n', flx);
fprintf('Errore assoluto = %.16g\n', abs_err);
fprintf('Errore relativo = %.16g\n', rel_err);
fprintf('u (unità arrotondamento teorica) = %.16g\n', u);
fprintf('Secondo più grande rappresentabile = %.16g\n', second_largest);

%% ES 2) Minori principali + LU pivoting per righe
fprintf('\n[ES2] Minori principali + LU con pivoting\n');
A = [1 1 1 1;
     1 1 1 2;
     1 1 0 1;
     0 1 1 1];

detA1 = det(A(1:1,1:1));
detA2 = det(A(1:2,1:2));
detA3 = det(A(1:3,1:3));
fprintf('det(A1)=%.16g, det(A2)=%.16g, det(A3)=%.16g\n', detA1, detA2, detA3);

% LU con pivoting per righe (parziale): P*A = L*U
[Llu, Ulu, P] = lu(A);
fprintf('P =\n'); disp(P);

% ricavo permutazione righe
perm = zeros(size(A,1),1);
for i=1:size(A,1)
    perm(i) = find(P(i,:)==1);
end
fprintf('Permutazione righe (P*A): [%s]\n', num2str(perm.'));
if all(perm == (1:size(A,1))')
    fprintf('Nessuno scambio di righe (pivoting non effettuato).\n');
else
    fprintf('Righe scambiate dove perm differisce dall''identità.\n');
end

%% ES 3) Vero/Falso (A SPD)
fprintf('\n[ES3] Affermazioni vere (A simmetrica definita positiva)\n');
% Risposte teoriche:
% 1) ||A^{-1}||_2 = 1/lambda_min(A) -> quindi FALSA
% 2) K2(A^{-1}) = K2(A) = lambda_max/lambda_min -> VERA
% 3) per SPD: singolari = autovalori (positivi) -> VERA
% 4) Gauss-Seidel converge per SPD -> VERA
% 5) Richardson stazionario con alpha_opt converge per SPD -> VERA
% 6) Gradiente coniugato converge per SPD (in aritmetica esatta in <=n) -> VERA
ans3 = [false true true true true true];
fprintf('V/F = [%s]\n', string(ans3));

%% ES 4) Metodi autovalori (spettri dati)
fprintf('\n[ES4] Affermazioni su metodi autovalori (spettri)\n');
% 1) Inverse power per minimo modulo di A: sì (min |lambda| = 0.1 unico) -> VERA
% 2) QR iterations per tutto spettro di B: sì -> VERA
% 3) Inverse power per minimo modulo di A e B: per B min |lambda| = 1 con tie (±1) => non garantito -> FALSA
% 4) Power method per massimo modulo di B e C: sì (4 unico) -> VERA
ans4 = [true true false true];
fprintf('V/F = [%s]\n', string(ans4));

%% ES 5) Gershgorin per colonne
fprintf('\n[ES5] Cerchi di Gershgorin per colonne\n');
A = [ 1  2  3  1;
      1 -2  1 -1;
     -1 -3 -2  1;
      2  2  2 10];

center = diag(A);
radiic = sum(abs(A - diag(center))); % raggi per colonne
rmax = max(radiic);

fprintf('Raggi (colonne) = [%s]\n', num2str(radiic));
fprintf('Raggio massimo (colonne) = %.16g\n', rmax);

% cerchio colonna con centro 10 (a44)
j = find(center == 10, 1);
rj = radiic(j);
fprintf('Cerchio colonna con centro 10 ha raggio = %.16g\n', rj);

% Numero di autovalori nel cerchio: se disgiunto dagli altri => 1 autovalore.
% Check disgiunzione su asse reale (matrice reale, dischi centrati su Re=diag):
intervals = [center - radiic; center + radiic]'; % [left,right] per ciascun disco
disc4 = intervals(j,:);
others = intervals; others(j,:) = [];
isDisjoint = all(disc4(1) > others(:,2)) || all(disc4(2) < others(:,1));
fprintf('Disco centro 10 disgiunto dagli altri? %d\n', isDisjoint);
if isDisjoint
    fprintf('=> Contiene esattamente 1 autovalore (criterio di Gershgorin).\n');
else
    fprintf('=> Non si può garantire 1; serve analisi di unione dischi.\n');
end

if doPlot && exist('gershcircles','file')
    gershcircles(A);
end

%% ES 6) Punto fisso (phi è Newton per f)
fprintf('\n[ES6] Punto fisso: phi(x) = x - f/f'' (Newton)\n');
phi = @(x) x - (x.^3 + 4*x - 8)./(3*x.^2 + 4);
x0 = 0.5; tol = 1e-8; nmax = 1000;

[x_seq, k] = fixed_point_driver(x0, phi, nmax, tol);
alpha = x_seq(end);
fprintf('alpha ~= %.16g\n', alpha);
fprintf('iterazioni = %d\n', k);
fprintf('ordine teorico atteso = 2 (Newton su radice semplice)\n');

%% ES 7) Punto fisso dipendente da gamma
fprintf('\n[ES7] Punto fisso: phi(x)=x + gamma/(2pi)*cos(pi x), alpha=1/2\n');
% Risposte teoriche:
% 1 vero: |phi'(alpha)|<1 <=> 0<gamma<4
% 2 vero: monotona locale se 0<phi'(alpha)<1 <=> 0<gamma<2
% 3 falso: gamma=2 => phi'(alpha)=0 e primo derivato non nullo è il terzo => ordine 3, non 1
ans7 = [true true false NaN];

% 4: calcolo numerico x^(10)
gamma = 1;
phi_g = @(x) x + gamma/(2*pi)*cos(pi*x);
x = 1/4;
for it=1:10
    x = phi_g(x);
end
ans7(4) = (abs(x - 0.5) < 1e-3);
fprintf('V/F = [%s]\n', string(ans7));
fprintf('Per gamma=1, x0=1/4: x10=%.16g, |x10-1/2|=%.3e\n', x, abs(x-0.5));

%% ES 8) Newton per minimo di Phi(y1,y2)=sin(y1)+cos(y2)
fprintf('\n[ES8] Newton per minimizzazione 2D\n');
gradPhi = @(y) [cos(y(1)); -sin(y(2))];
hessPhi = @(y) [-sin(y(1)) 0; 0 -cos(y(2))];

xk = [4;2];
fprintf('x0 = [%g,%g]^T\n', xk(1), xk(2));
for it=1:2
    g = gradPhi(xk);
    H = hessPhi(xk);
    step = H \ g;
    xk = xk - step;
end
fprintf('x2 = [%0.16g, %0.16g]^T\n', xk(1), xk(2));
fprintf('gradPhi(y) = [cos(y1); -sin(y2)]\n');
fprintf('HessPhi(y) = diag(-sin(y1), -cos(y2))\n');

%% PARTE: A pentadiagonale n=100
fprintf('\n[PARTE PENTADIAG] n=100, A=pentadiag(1,-4,9,-4,1)\n');
n = 100;
A = pentadiag_sp(1,-4,9,-4,1,n);
x_true = ones(n,1);
b = A*x_true;

% Punto 1: Cholesky + errore
fprintf('\n[Punto 1] Cholesky solve + errore\n');
R = chol(A);                 % A = R'*R
x_hat = R \ (R' \ b);         % solve
err2 = norm(x_hat - x_true);
relerr2 = err2 / norm(x_true);
fprintf('||x_hat-x||_2 = %.3e, rel = %.3e\n', err2, relerr2);

% (flop count qui lasciato come commento / stima: fattorizzazione bandata ~ O(n))
fprintf('Nota: matrice pentadiagonale => Cholesky bandato ha costo O(n).\n');

% Punto 2: Gauss-Seidel + 100 iter da x0=b + stima ordine
fprintf('\n[Punto 2] Gauss-Seidel 100 iter, x0=b\n');
x0 = b;
if exist('gs','file')
    [x_gs, k_gs, ~] = gs(A, b, x0, 0, 100); % tol=0 => forza 100 iter se implementato così
else
    [x_gs, k_gs] = gs_basic(A,b,x0,100);
end
fprintf('Iterazioni effettuate (richieste 100): %d\n', k_gs);

% Sequenza errori per stimare ordine (qui calcolo errori ad ogni iterazione con implementazione interna)
[x_seq_gs, err_seq] = gs_history(A,b,x0,100,x_true);
if exist('data_analysis_toolbox','file')
    [p_est, C_est, ~] = data_analysis_toolbox('roots', err_seq);
    fprintf('Stima ordine (asintotico) p ~= %.3f, C ~= %.3g\n', p_est, C_est);
else
    fprintf('data_analysis_toolbox.m non trovato: salto stima p.\n');
end

% Punto 3: precondizionatore P(gamma) e richardson (gradiente precondizionato)
fprintf('\n[Punto 3] Gamma_opt numerico + richardson precondizionato\n');
gammas = linspace(2,4,81);
kappa = zeros(size(gammas));
for i=1:numel(gammas)
    P = precond_P(gammas(i), n);
    % cond del problema generalizzato A v = lambda P v
    % Per n=100, dense eig è abbastanza veloce e più robusto di eigs
    eigvals = eig(full(A), full(P));
    lmax = max(eigvals);
    lmin = min(eigvals);
    kappa(i) = lmax / lmin;
end
[~, idx] = min(kappa);
gamma_opt = gammas(idx);
fprintf('gamma_opt (scan) ~= %.4f, kappa ~= %.4g\n', gamma_opt, kappa(idx));

Popt = precond_P(gamma_opt, n);
x0 = b;
tol = 1e-6; nmax = 5000;
if exist('richardson','file')
    [xk, k] = richardson(A, b, Popt, x0, tol, nmax); % alpha omesso => dinamico
else
    error('richardson.m non trovato');
end
errA = sqrt((xk-x_true)' * A * (xk-x_true));
resN = norm(A*xk - b) / norm(b);
fprintf('k = %d, ||xk-x||_A = %.3e, residuo norm = %.3e\n', k, errA, resN);

% Punto 4: Rayleigh quotient gradient per lambda_min
fprintf('\n[Punto 4] Rayleigh-grad per lambda_min (50 iter)\n');
x0 = (1:n)';  iters = 50;
[lambda50, lambda_hist] = rayleigh_grad_min(A, x0, iters);
fprintf('lambda^(50) ~= %.16g\n', lambda50);

% stima ordine su |lambda_k - lambda_50| (proxy)
errL = abs(lambda_hist - lambda50);
errL = errL(errL>0);
if numel(errL) >= 6 && exist('data_analysis_toolbox','file')
    [p_estL, ~, ~] = data_analysis_toolbox('roots', errL);
    fprintf('Stima ordine su lambda (proxy) p ~= %.3f\n', p_estL);
end

%% Interpolazione / integrazione / LSQ
fprintf('\n[Parte dati]\n');

% Interpolazione grado 4
xD = [0 1 2 2.5 3]; yD = [1.5 -0.5 1.5 -0.5 2];
p4 = polyfit(xD, yD, 4);
val = polyval(p4, 1.5);
p4_dd = polyder(polyder(p4));
val2 = polyval(p4_dd, 1.5);
fprintf('Pi4(1.5)=%.16g, d2Pi4/dx2(1.5)=%.16g\n', val, val2);

% Integrale interpolante lineare a tratti = trapezi
xH = [0 0.25 0.5 0.75 1]; yH = [2 0.5 -1.5 0.5 2];
IH = trapz(xH, yH);
fprintf('Integral_0^1 Pi1^H(x) dx = %.16g\n', IH);

% LSQ grado 3 su 7 punti
xLS = [0 0.5 1 1.5 2 2.5 3];
yLS = [4 0.25 0.5 -1.5 0.5 1 -0.5];
if exist('data_analysis_toolbox','file')
    [coeffs, resid_sq] = data_analysis_toolbox('lsq', xLS, yLS, 3);
else
    coeffs = polyfit(xLS, yLS, 3);
    yfit = polyval(coeffs, xLS);
    resid_sq = sum((yLS - yfit).^2);
end
fprintf('LSQ deg3 residuo quadratico = %.16g\n', resid_sq);
fprintf('m per residuo nullo (interpolazione esatta) = 6\n');

% Trapezio composito su e^{-x} con raddoppio N finché |IN - I_{N/2}|<tol
fprintf('\n[Trapezio adattivo] e^{-x}\n');
tol = 1e-6;
f = @(x) exp(-x);
N = 2;
Iprev = trap_comp(0,1,N,f);
while true
    N = 2*N;
    Icur = trap_comp(0,1,N,f);
    if abs(Icur - Iprev) < tol, break; end
    Iprev = Icur;
end
fprintf('N minimo = %d\n', N);

%% Heun stabilità (vero/falso)
fprintf('\n[Heun stabilita] V/F\n');
% 1 falso, 2 falso, 3 vero (stesso intervallo su asse reale negativo di FE), 4 falso
ansHeun = [false false true false];
fprintf('V/F = [%s]\n', string(ansHeun));

%% BVP (0,4) con Neumann in x=4
fprintf('\n[BVP 0-4] -u'''' + 6u'' + 2(4-x)u=0, u(0)=1, u''(4)=0 (Neumann)\n');
h = 0.1;
[u, xgrid] = bvp_fd_0_4(h);
u40 = u(end); % x=4
fprintf('u_40 ~= u(4) ~= %.16g\n', u40);

% u'(0.5) centrata
j = round(0.5/h) + 1; % index MATLAB (x=0 -> 1)
up05 = (u(j+1) - u(j-1)) / (2*h);
fprintf('u''(0.5) (centrata) ~= %.16g\n', up05);

%% BVP (0,1) upwind
fprintf('\n[BVP 0-1 upwind] -2u'''' + 500u'' = 12, u(0)=0, u(1)=2\n');
h = 1/20;
[u2, x2] = bvp_upwind_0_1(h);
u19 = u2(end-1); % x19
fprintf('u_19 ~= u(x19=0.95) ~= %.16g\n', u19);

%% PDE calore: usa EqCalore_DiffFin_Theta (CN = theta=0.5)
fprintf('\n[PDE Calore CN]\n');
mu = 1/2;
fxt = @(x,t) 0*x;
a=0; bdom=2;
us = @(t) 0*t;
ud = @(t) 0*t;
g0 = @(x) 4*x.^2 .* (2 - x);
T = 1; hx = 0.5; dt = 0.05; theta = 0.5;

if exist('EqCalore_DiffFin_Theta','file')
    [U, xg, tg] = EqCalore_DiffFin_Theta(mu, fxt, a, bdom, us, ud, g0, T, hx, dt, theta);
    % u(0.5,1): x index 2, t index 21
    u_05_1 = U(2, 21);
    fprintf('u(0.5,1) ~= %.16g\n', u_05_1);
else
    fprintf('EqCalore_DiffFin_Theta.m non trovato: salto.\n');
end

%% ODE sistema lineare: CN efficiente + Simpson su p(t)
fprintf('\n[ODE lineare] CN efficiente + Simpson su p(t)\n');
A4 = [-4 -17 0 5;
       1   0 0 0;
       0   2 -3 -23;
       0   0 1 0];
g = @(t) [(-5*cos(t) + (285/4)*sin(t))*exp(-t/2);
          0;
          (83*cos(t) - 18*sin(t))*exp(-t/2);
          0];
y0 = [5;0;-2;4];
tf = 10; h = 0.2;

[tv, Ucn] = cn_linear_affine(A4, g, y0, tf, h);
u1 = Ucn(:,2);
uN = Ucn(:,end);
fprintf('u1 = \n'); disp(u1);
fprintf('uNh = \n'); disp(uN);

pvals = Ucn(2,:); % p(t) è la seconda componente
I_p = simpson_samples(tv, pvals);
fprintf('Simpson int_0^{tf} p(t) dt ~= %.16g\n', I_p);

%% ODE Punto 2: errori Eh su h_i e stima ordine
fprintf('\n[ODE Punto2] Errori Eh e stima ordine\n');
hs = (2.^(-(1:4))) / 10; % 0.05, 0.025, ...
Eh = zeros(size(hs));
for ii=1:numel(hs)
    [tvi, Ui] = cn_linear_affine(A4, g, y0, tf, hs(ii));
    % Yi rimosso perche non definito e non usato
    % calcolo max errore
    Ei = 0;
    for k=1:numel(tvi)
        yex = exact_y_vec(tvi(k));
        Ei = max(Ei, norm(Ui(:,k) - yex));
    end
    Eh(ii) = Ei;
end
fprintf('h = [%s]\n', num2str(hs));
fprintf('Eh = [%s]\n', num2str(Eh));
if exist('data_analysis_toolbox','file')
    [p_est, ~, ~] = data_analysis_toolbox('roots', Eh);
    fprintf('Stima ordine CN p ~= %.3f\n', p_est);
end

%% ODE Punto 3: stabilità assoluta (g=0) per FE, CN, Heun
fprintf('\n[ODE Punto3] Stabilita assoluta per g=0\n');
lam = eig(A4);
fprintf('eig(A) = \n'); disp(lam);

% CN A-stable: se Re(lam)<=0 => stabile per ogni h
if all(real(lam) <= 0)
    fprintf('CN: stabile per ogni h>0 (Re(lambda)<=0).\n');
else
    fprintf('CN: non garantita per ogni h (ci sono Re(lambda)>0).\n');
end

% Per FE e Heun stampo h_max numerico via scan
hscan = linspace(0,5,20001);
hmax_FE = max_h_stable(hscan, lam, @(z) 1+z);
hmax_Heun = max_h_stable(hscan, lam, @(z) 1+z+0.5*z.^2);
fprintf('FE: h_max ~ %.4f (scan)\n', hmax_FE);
fprintf('Heun: h_max ~ %.4f (scan)\n', hmax_Heun);

%% ODE Punto 4: metodo multipasso implicito lineare (AM2-like)
fprintf('\n[ODE Punto4] Metodo multipasso implicito lineare\n');
% è implicito perché compare f_{n+1} = A u_{n+1} + g(t_{n+1})
fprintf('Il metodo è IMPLICITO (compare f_{n+1}).\n');

% Calcolo u1 con CN, poi applico multipasso
h = 0.2;
[tv, Ucn] = cn_linear_affine(A4, g, y0, tf, h);
u0 = y0;
u1 = Ucn(:,2);

[tv2, Uam] = am2_like_linear(A4, g, u0, u1, tf, h);
fprintf('u2 =\n'); disp(Uam(:,3));
fprintf('uNh =\n'); disp(Uam(:,end));

fprintf('\n=== FINE ===\n');

%% ====== FUNZIONI LOCALI ======
function [x_seq, k] = fixed_point_driver(x0, phi, nmax, tol)
    if exist('ptofis','file')
        [x_seq, k] = ptofis(x0, phi, nmax, tol);
        if iscell(x_seq), x_seq = x_seq{1}; end
        if size(x_seq,2) > 1, x_seq = x_seq(:); end
    else
        x_seq = zeros(nmax+1,1); x_seq(1)=x0;
        for k=1:nmax
            x_seq(k+1) = phi(x_seq(k));
            if abs(x_seq(k+1)-x_seq(k)) < tol, x_seq = x_seq(1:k+1); return; end
        end
    end
end

function A = pentadiag_sp(a2,a1,a0,b1,b2,n)
    % diagonali: -2,-1,0,+1,+2
    A = spdiags([a2*ones(n,1), a1*ones(n,1), a0*ones(n,1), b1*ones(n,1), b2*ones(n,1)], ...
                -2:2, n, n);
end

function P = precond_P(gamma,n)
    % tridiagonale con 2 ai bordi e gamma interno, -1 off-diag
    d = gamma*ones(n,1); d(1)=2; d(end)=2;
    P = spdiags([-ones(n,1), d, -ones(n,1)], -1:1, n, n);
end

function [x_seq, err_seq] = gs_history(A,b,x0,kmax,x_true)
    x = x0;
    n = length(b);
    err_seq = zeros(kmax,1);
    x_seq = zeros(n,kmax+1);
    x_seq(:,1)=x0;
    for k=1:kmax
        for i=1:n
            s1 = A(i,1:i-1)*x(1:i-1);
            s2 = A(i,i+1:n)*x(i+1:n);
            x(i) = (b(i)-s1-s2)/A(i,i);
        end
        x_seq(:,k+1)=x;
        err_seq(k) = norm(x - x_true);
    end
end

function [x, k] = gs_basic(A,b,x0,kmax)
    x=x0;
    n=length(b);
    for k=1:kmax
        for i=1:n
            s1 = A(i,1:i-1)*x(1:i-1);
            s2 = A(i,i+1:n)*x(i+1:n);
            x(i) = (b(i)-s1-s2)/A(i,i);
        end
    end
end

function I = trap_comp(a,b,N,f)
    x = linspace(a,b,N+1);
    y = f(x);
    h = (b-a)/N;
    I = h*(0.5*y(1) + sum(y(2:end-1)) + 0.5*y(end));
end

function [u,x] = bvp_fd_0_4(h)
    a=0; b=4;
    x = (a:h:b).';
    m = length(x); % include endpoints
    % Unknowns u0..u_{m-1}, but u0 fixed Dirichlet and Neumann at end
    u0 = 1;

    % Build system for unknowns u1..u_{m-1} (last included)
    n = m-1;
    A = sparse(n,n);
    rhs = zeros(n,1);

    % Neumann at x=b: u'(b)=0 => (u_m - u_{m-1})/h = 0 but u_m doesn't exist (b is last)
    % With nodes x_{m} absent; we have last node is x_{m-1}=b:
    % use backward diff: (u_{m-1} - u_{m-2})/h = 0 => u_{m-1} = u_{m-2}
    % We'll enforce this as last equation.
    for j=2:m-1  % interior nodes indices in MATLAB: 2..m-1 correspond to x=h..b-h
        xi = x(j);
        row = j-1; % maps u_j to unknown index (u1->1, ..., u_{m-1}->n)

        % centered second derivative u'' ~ (u_{j-1}-2u_j+u_{j+1})/h^2
        % centered first derivative u' ~ (u_{j+1}-u_{j-1})/(2h)
        % equation: -u'' + 6u' + 2(4-x)u = 0
        cjm1 = -(1/h^2) + 6*(-1/(2*h));
        cj   = -(-2/h^2) + 2*(4-xi);
        cjp1 = -(1/h^2) + 6*( 1/(2*h));

        % u_{j-1}
        if j-1 == 1
            rhs(row) = rhs(row) - cjm1*u0; % u0 known
        else
            A(row, (j-2)) = A(row, (j-2)) + cjm1; % u_{j-1} unknown
        end
        % u_j
        A(row, (j-1)) = A(row, (j-1)) + cj;
        % u_{j+1}
        A(row, (j))   = A(row, (j))   + cjp1;
    end

    % first equation for j=1 (x=h): use same stencil with u0 known and u2 unknown
    j=2; xi=x(j); row=1;
    cjm1 = -(1/h^2) + 6*(-1/(2*h));
    cj   = -(-2/h^2) + 2*(4-xi);
    cjp1 = -(1/h^2) + 6*( 1/(2*h));
    rhs(row) = rhs(row) - cjm1*u0;
    A(row,1) = A(row,1) + cj;
    A(row,2) = A(row,2) + cjp1;

    % Neumann last: u_{m-1} - u_{m-2} = 0
    A(n,n) = 1;
    A(n,n-1) = -1;
    rhs(n) = 0;

    u_unknown = A\rhs; % gives u1..u_{m-1}
    u = [u0; u_unknown];
end

function [u,x] = bvp_upwind_0_1(h)
    a=0; b=1;
    x = (a:h:b).';
    m = length(x);
    u0 = 0; uN = 2;

    % unknowns u1..u_{m-2} internal, size n=m-2
    n = m-2;
    A = sparse(n,n);
    rhs = 12*ones(n,1);

    p = 500; mu = 2; % equation: -2 u'' + 500 u' = 12
    for j=2:m-1 % internal nodes (MATLAB)
        row = j-1; % maps u_j -> unknown index (u1->1)
        % u'' centered: (u_{j-1}-2u_j+u_{j+1})/h^2
        % u' upwind (p>0): (u_j - u_{j-1})/h
        cjm1 = -mu*(1/h^2) + p*(-1/h);
        cj   = -mu*(-2/h^2) + p*(1/h);
        cjp1 = -mu*(1/h^2);

        % u_{j-1}
        if j-1 == 1
            rhs(row) = rhs(row) - cjm1*u0;
        else
            A(row, (j-2)) = A(row, (j-2)) + cjm1;
        end
        % u_j
        A(row, (j-1)) = A(row, (j-1)) + cj;
        % u_{j+1}
        if j+1 == m
            rhs(row) = rhs(row) - cjp1*uN;
        else
            A(row, (j)) = A(row, (j)) + cjp1;
        end
    end

    u_int = A\rhs;
    u = [u0; u_int; uN];
end

function [t, U] = cn_linear_affine(A, g, y0, tf, h)
    N = round(tf/h);
    t = (0:N)*h;
    n = length(y0);
    U = zeros(n, N+1);
    U(:,1) = y0;

    I = speye(n);
    M = (I - (h/2)*A);
    K = (I + (h/2)*A);
    % LU una volta (efficiente)
    [L,Uu,P] = lu(M);

    for k=1:N
        tk = t(k); tk1 = t(k+1);
        rhs = K*U(:,k) + (h/2)*(g(tk) + g(tk1));
        U(:,k+1) = Uu \ (L \ (P*rhs));
    end
end

function I = simpson_samples(t, y)
    % Simpson composito su griglia uniforme
    h = t(2)-t(1);
    N = numel(t)-1;
    if mod(N,2) ~= 0
        error('Simpson richiede N pari (numero sottointervalli).');
    end
    I = h/3*(y(1) + y(end) + 4*sum(y(2:2:end-1)) + 2*sum(y(3:2:end-2)));
end

function y = exact_y_vec(t)
    % y = (p_dot, p, q_dot, q)
    p = 5*exp(-t/2)*sin(t);
    pd = 5*exp(-t/2)*(cos(t) - 0.5*sin(t));
    q = 4*exp(-t/2)*cos(t);
    qd = 4*exp(-t/2)*(-sin(t) - 0.5*cos(t));
    y = [pd; p; qd; q];
end

function hmax = max_h_stable(hscan, lam, R)
    ok = false(size(hscan));
    for i=1:numel(hscan)
        z = hscan(i)*lam;
        ok(i) = all(abs(R(z)) <= 1 + 1e-12);
    end
    idx = find(ok, 1, 'last');
    if isempty(idx), hmax = 0; else, hmax = hscan(idx); end
end

function [lambdaK, lambda_hist] = rayleigh_grad_min(A, x0, K)
    x = x0;
    lambda_hist = zeros(K+1,1);
    y = x / norm(x);
    lambda = y'*(A*y);
    lambda_hist(1) = lambda;

    for k=1:K
        r = (A - lambda*speye(size(A,1))) * y;
        num = (r'*r);
        den = (r'*(A*r));
        alpha = num/den;
        x = y - alpha*r;
        y = x / norm(x);
        lambda = y'*(A*y);
        lambda_hist(k+1) = lambda;
    end
    lambdaK = lambda_hist(end);
end

function [t, U] = am2_like_linear(A, g, u0, u1, tf, h)
    % u_{n+1} = u_n + h/12(-f_{n-1}+8f_n+5f_{n+1})
    % (I - 5h/12 A) u_{n+1} = u_n + h/12(-f_{n-1}+8f_n+5 g_{n+1}) + (5h/12)g_{n+1} already in f_{n+1}
    % con f_n = A u_n + g(t_n)
    N = round(tf/h);
    t = (0:N)*h;
    n = length(u0);
    U = zeros(n, N+1);
    U(:,1)=u0; U(:,2)=u1;

    I = speye(n);
    M = (I - (5*h/12)*A);
    [L,Uu,P] = lu(M);

    f0 = A*U(:,1) + g(t(1));
    f1 = A*U(:,2) + g(t(2));

    for k=2:N
        tkp1 = t(k+1);
        gkp1 = g(tkp1);
        rhs = U(:,k) + (h/12)*(-f0 + 8*f1 + 5*gkp1); % parte nota, manca 5*(A u_{k+1})
        % Risolvo: (I - 5h/12 A)u_{k+1} = rhs
        U(:,k+1) = Uu \ (L \ (P*rhs));

        % aggiorno f_{k-1}, f_k
        f0 = f1;
        f1 = A*U(:,k+1) + gkp1;
    end
end