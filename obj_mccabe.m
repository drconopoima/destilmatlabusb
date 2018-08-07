function [ salida ] = obj_mccabe(y, zf, q, P, mezcla, light, heavy, RM, kij)
%
%“obj_mccabe.m”. [salida] = obj_mccabe(y, zf, q, P, mezcla, light, heavy, RM, kij).
%Función objetivo para obtener mediante un fsolve() la fracción molar de clave
% liviano en la mezcla vapor, según el equilibrio de fases el cual predice la curva 
%y según el modelo termodinámico provisto “RM”, para obtener la ordenada y en dónde
% una recta q de alimentación de fracción molar de liviano zf, condición térmica q 
%y a presión P corta el equilibrio. 
%
%
%
%   Luis Jesús Díaz Manzo
comp = Sustancia.empty(0,2);
comp(1) = mezcla.comp(light);
comp(2) = mezcla.comp(heavy);
if nargin < 10
    kij = 0;
end
gas = Mezcla(comp, [y, 1-y], kij);
[~, x, ~, flag] = RM.DewT(P, gas);
x = x(1);
salida = (y(1) - (q.*x)./(q - 1) + (zf./(q-1)));
end

