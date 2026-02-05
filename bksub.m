function x = bksub(A,b)
%
% x = bksub(A,b)
%
% Algoritmo delle sostituzioni all'indietro
%
% Parametri di ingresso:
% A    Matrice quadrata triangolare superiore del sistema lineare Ax=b
% b    Termine noto, vettore colonna
%
% Parametri di uscita:
% x    Soluzione del sistema lineare Ax=b
%
%                                         Politecnico di Milano, 04/04/2024
%

n=length(b);

if((size(A,1)~=n)||(size(A,2)~=n))
  error('ERRORE: dimensioni incompatibili')
end

if(~isequal(A,triu(A)))
  error('ERRORE: matrice non triangolare superiore')
end


if(prod(diag(A))==0)
% almeno un elemento diagonale nullo
  error('ERRORE: matrice singolare')
end

x=zeros(n,1);

x(n) = b(n)/A(n,n);

for i=n-1:-1:1
   x(i)=(b(i)-A(i,i+1:n)*x(i+1:n))/A(i,i);
end







