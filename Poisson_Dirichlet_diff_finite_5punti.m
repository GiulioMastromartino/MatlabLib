function [ X, Y, U ] = Poisson_Dirichlet_diff_finite_5punti( mu, f, g, Omega, hx, hy )
%
% [ X, Y, U ] = Poisson_Dirichlet_diff_finite_5punti( mu, f, g, Omega, hx, hy )
% 
% Schema alle differenze finite a 5 punti per l'approssimazione del 
% problema di Poisson-Dirichlet a coefficienti costanti in un 
% dominio rettangolare  
%
% - mu * laplaciano(u) = f(x,y)   in Omega = ( a, b ) x ( c, d )
% u(x,y) = g(x,y)                 su bordo di Omega
%
% Parametri di ingresso:
% mu      Coefficiente diffusione (scalare positivo)
% f       Function handle @(x,y) per termine forzante
% g       Function handle @(x,y) per dato di Dirichlet al bordo
% Omega   Vettore dei limiti del dominio rettangolare [ a, b, c, d ]
% hx      Passo griglia in direzione x
% hy      Passo griglia in direzione y
% 
% Parametri di uscita:
% X, Y    Matrici (Nx+2)x(Ny+2) delle coordinate dei nodi
% U       Matrice (Nx+2)x(Ny+2) della soluzione approssimata nei nodi
%
%                                         Politecnico di Milano, 12/06/2024
%

% --- 1. Griglia computazionale ---
x_start = Omega(1); x_end = Omega(2);
y_start = Omega(3); y_end = Omega(4);

% Numero intervalli
Nx = round((x_end - x_start) / hx);
Ny = round((y_end - y_start) / hy);

% Nodi interni + bordo (Nx+1 nodi in x, Ny+1 nodi in y)
% La funzione originale usava Nx come nodi interni, qui adattiamo 
% la logica per essere consistenti con Nx come numero di intervalli
xnodes = linspace(x_start, x_end, Nx + 1);
ynodes = linspace(y_start, y_end, Ny + 1);

% Meshgrid per output e calcoli
[X, Y] = meshgrid(xnodes, ynodes);
X = X'; Y = Y'; % Trasponiamo per avere indici (i,j) coerenti con (x,y)
% Ora X e Y sono (Nx+1) x (Ny+1)

% Numero totale nodi
num_nodes = (Nx + 1) * (Ny + 1);

% --- 2. Assemblaggio Sistema Lineare Au = b ---
% Inizializzazione sparse per efficienza
% Usiamo la numerazione lessicografica: nodo (i,j) -> k = i + (j-1)*(Nx+1)
% con i=1..Nx+1, j=1..Ny+1

% Stima non-zeri: 5 diagonali per nodi interni, 1 per nodi bordo
A = spalloc(num_nodes, num_nodes, 5*num_nodes);
b = zeros(num_nodes, 1);

dx2 = hx^2;
dy2 = hy^2;
coeff_center = 2*mu/dx2 + 2*mu/dy2;
coeff_x = -mu/dx2;
coeff_y = -mu/dy2;

for j = 1 : Ny + 1
    for i = 1 : Nx + 1
        k = i + (j - 1) * (Nx + 1); % Indice globale
        
        % Coordinate nodo
        xi = xnodes(i);
        yi = ynodes(j);
        
        % Controllo se nodo è al bordo
        if (i == 1 || i == Nx + 1 || j == 1 || j == Ny + 1)
            % Equazione identità per Dirichlet: u_k = g(xi, yi)
            A(k, k) = 1.0;
            b(k) = g(xi, yi);
        else
            % Nodo interno: schema a 5 punti
            % -mu * ( (u_{i+1,j}-2u+u_{i-1,j})/hx^2 + (u_{i,j+1}-2u+u_{i,j-1})/hy^2 ) = f
            
            % Indici vicini
            k_nord = k + 1;             % (i+1, j) -> errore nella logica originale, corretto sotto
            k_sud  = k - 1;             % (i-1, j)
            k_est  = k + (Nx + 1);      % (i, j+1)
            k_ovest= k - (Nx + 1);      % (i, j-1)
            
            % Nota: la meshgrid trasposta implica X(i,j) = x_i, Y(i,j) = y_j
            % k corrisponde a (i,j).
            % Vicino est (i+1, j) è k+1 solo se ordiniamo per x veloce
            % La numerazione k = i + (j-1)*(Nx+1) scorre prima le x (i) poi le y (j)
            
            k_xp1 = i + 1 + (j - 1) * (Nx + 1); % (i+1, j)
            k_xm1 = i - 1 + (j - 1) * (Nx + 1); % (i-1, j)
            k_yp1 = i + (j) * (Nx + 1);         % (i, j+1)
            k_ym1 = i + (j - 2) * (Nx + 1);     % (i, j-1)
            
            A(k, k)     = coeff_center;
            A(k, k_xp1) = coeff_x;
            A(k, k_xm1) = coeff_x;
            A(k, k_yp1) = coeff_y;
            A(k, k_ym1) = coeff_y;
            
            b(k) = f(xi, yi);
        end
    end
end

% --- 3. Soluzione ---
u_vec = A \ b;

% --- 4. Reshape in matrice U ---
U = zeros(Nx + 1, Ny + 1);
for j = 1 : Ny + 1
    for i = 1 : Nx + 1
        k = i + (j - 1) * (Nx + 1);
        U(i, j) = u_vec(k);
    end
end

end