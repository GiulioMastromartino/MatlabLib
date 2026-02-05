function x = fwsub(A,b)
%
% x = fwsub(A,b)
%
% Algoritmo delle sostituzioni in avanti
%
% Parametri di ingresso:
% A    Matrice quadrata triangolare inferiore del sistema lineare Ax=b
% b    Termine noto, vettore colonna
%
% Parametri di uscita:
% x    Soluzione del sistema lineare Ax=b
%
%                                         Politecnico di Milano, 04/04/2024
%

n=length(b);

if ((size(A,1)~=n)||(size(A,2)~=n))
  error('ERRORE: dimensioni incompatibili')
end

if ~isequal(A,tril(A))
  error('ERRORE: matrice non triangolare inferiore')
end


if (prod(diag(A))==0)
% almeno un elemento diagonale nullo
  error('ERRORE: matrice singolare')
end

x=zeros(n,1);

x(1)=b(1)/A(1,1);

for i=2:n
  x(i)=(b(i)-A(i,1:i-1)*x(1:i-1))/A(i,i);
end



