function [ obj ] = interpolLagrange2( x, Ex, Ey, A, B)
%InterpolLagrange2 Realiza la interpolación Lagrange de donde corta la
%recta que pasa por los puntos A,B: [xA, yA] y [xB, yB] a la curva de 2do
%orden de Lagrande que tiene los puntos Ex y Ey
%
%Función objetivo para realizar mediante un solve() interpolación de 
%Lagrange cuadrática que halle el punto horizontal x donde corta una 
%línea recta que atraviesa los vectores de coordenadas A, B [1x2] sobre 
%una curva construida construyendo pares de coordenadas [x, y] tomando
% valores de los vectores Ex, Ey de mínimo 3 elementos y arbitrariamente 
%largos. Utilizado en la construcción del diagrama de Hengstebeck para la
% ubicación de la recta q en alimentaciones bifásicas.
%
%
%   Luis Jesús Díaz

m = (A(2) - B(2))/(A(1) - B(1));
b = -m*A(1) + A(2);
Eyobj = m*x + b;
max = find(Ex > x, 1, 'first');

try 
    obj = Eyobj - (Ey(max - 3)*(((x - Ex(max - 2))*(x - Ex(max)))/((Ex(max - 3) - Ex(max - 2))*(Ex(max - 3) - Ex(max)))) + Ey(max - 2)*(((x - Ex(max - 3))*(x - Ex(max)))/((Ex(max - 2) - Ex(max - 3))*(Ex(max - 2) - Ex(max)))) + Ey(max)*(((x - Ex(max - 3))*(x - Ex(max - 2)))/((Ex(max) - Ex(max-3))*(Ex(max) - Ex(max - 2))))); 
catch ME 
    min = max;
    obj = Eyobj - (Ey(min)*(((x - Ex(min + 1))*(x - Ex(min + 3)))/((Ex(min) - Ex(min + 1))*(Ex(min) - Ex(min + 3)))) + Ey(min + 1)*(((x - Ex(min))*(x - Ex(min + 3)))/((Ex(min + 1) - Ex(min))*(Ex(min + 1) - Ex(min + 3)))) + Ey(min + 3)*(((x - Ex(min))*(x - Ex(min + 1)))/((Ex(min + 3) - Ex(min))*(Ex(min + 3) - Ex(min + 1))))); 
end
end

