function [x, matriz_tridiagbloque] = tridiagThomas( a, b, c, f )
%
%[ΔX, Fik] = tridiagThomas(Bik, Aik, Cik, Dk).
%
%Función que recibe vectores o matrices bloque Bik, Aik, Cik, Dk con los 
%cuales construye la matriz tridiagonal Fik del sistema de ecuaciones matricial 
%Fik• ΔX=Dk. Los vectores o matrices Bik, Aik, Cik deben ser fila, horizontales
% e internamente se distribuyen a los elementos correspondientes en la matriz 
%tridiagonal. El sistema es resuelto con  ΔX=Fik\Dk
%
%  Solve the  n x n  tridiagonal system for x:
%
%  [ a(1)  c(1)                                  ] [  x(1)  ]   [  f(1)  ]
%  [ b(2)  a(2)  c(2)                            ] [  x(2)  ]   [  f(2)  ]
%  [       b(3)  a(3)  c(3)                      ] [        ]   [        ]
%  [            ...   ...   ...                  ] [  ...   ] = [  ...   ]
%  [                    ...    ...    ...        ] [        ]   [        ]
%  [                        b(n-1) a(n-1) c(n-1) ] [ x(n-1) ]   [ f(n-1) ]
%  [                                 b(n)  a(n)  ] [  x(n)  ]   [  f(n)  ]
%
%  f must be a vector or block of square matrices (row or column) of length n
%  a must be a vector or block of square matrices of length n 
%  b must be a vector or block of square matrices of length n-1| b(2)...b(n)
%  c must be a vector or block of square matrices of length n-1| c(1)...c(n-1)


tamano = size(a);
if tamano(1) > tamano(2);
    a = a';
end
tamano = size(c);
if tamano(1) > tamano(2);
    c = c';
end
tamano = size(f);
if tamano(1) > tamano(2);
    f = f';
end
tamano = size(a);
if tamano(1) > tamano(2)
    b = b';
end
tamano = size(a);
m = tamano(2);
len_sub_j = tamano(1);
intervalo = len_sub_j - 1;
n = m/tamano(1);
matriz_tridiagbloque = zeros(m,m);
for iter = 0:n-1
    j = len_sub_j*iter + 1;
    matriz_tridiagbloque(j:j+intervalo,j:j+intervalo) = a(1:1+intervalo,j:j+intervalo);
end
iter = 0;
for i = 1:n-1
    j = len_sub_j*(i) + 1;
    k = len_sub_j*(iter) + 1;
    iter = iter + 1;
    matriz_tridiagbloque(j:j+intervalo,k:k+intervalo) = b(1:1+intervalo,k:k+intervalo);
end
iter = 0;
for i = 1:n-1
    j = len_sub_j*(i) + 1;
    k = len_sub_j*(iter) + 1;
    iter = iter + 1;
    matriz_tridiagbloque(k:k+intervalo,j:j+intervalo) = c(1:1+intervalo,k:k+intervalo);
end
if tamano(1) < tamano(2);
    f = f';
end
x = matriz_tridiagbloque\f;

