function [L, U, x] = thomas(A, b)
%
% [L, U, x] = thomas(A, b)
% 
% Applica l'algoritmo di Thomas per risolvere un sistema lineare Ax=b con
% matrice A tridiagonale e determina i fattori L ed U di A applicando
% l'algoritmo (senza pivoting).
% 
% Parametri di ingresso:
% A    Matrice tridiagonale quadrata
% b    Termine noto, vettore colonna
% 
% Parametri di uscita:
% L    Matrice triangolare inferiore della fattorizzazione LU di A
% U    Matrice triangolare superiore della fattorizzazione LU di A
% x    Vettore soluzione del sistema Ax = b 
%
%                                         Politecnico di Milano, 04/04/2024
%

%%%%%%%%%%%%%%%%%
% PREPROCESSING %
%%%%%%%%%%%%%%%%%
% dimensione del sistema
N = size (A, 1);

% crea vettore dei parametri alpha della fattorizzazione di Thomas
alpha = zeros (N, 1);

% crea vettore dei parametri delta della fattorizzazione di Thomas
delta = zeros (N-1, 1);

% estrae sopradiagonale, sottodiagonale e diagonale principale di A
c = diag (A, 1);
e = diag (A, -1);
a = diag (A, 0);


%%%%%%%%%%%%%%%%%%%
% FATTORIZZAZIONE %
%%%%%%%%%%%%%%%%%%%
alpha(1) = a(1);

for i = 2:N
      delta(i-1) = e(i-1) / alpha(i-1);
      alpha(i) = a(i) - delta(i-1) * c(i-1);
end

L = diag (ones (N,1), 0) + diag (delta, -1);
U = diag (alpha, 0) + diag (c, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLUZIONE DEL SISTEMA %
%%%%%%%%%%%%%%%%%%%%%%%%%
y = zeros(N,1);
y(1) = b(1);

for i = 2:N
      y(i) = b(i) - delta(i-1) * y(i-1);
end

x = zeros(N,1);
x(N) = y(N) / alpha(N);

for i = N-1:-1:1
      x(i) = (y(i) - c(i) * x(i+1)) / alpha(i);
end