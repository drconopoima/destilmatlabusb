function [ salida ] = obj_mccabe(y, zf, q, P, mezcla, light, heavy, RM, kij)
%
%�obj_mccabe.m�. [salida] = obj_mccabe(y, zf, q, P, mezcla, light, heavy, RM, kij).
%Funci�n objetivo para obtener mediante un fsolve() la fracci�n molar de clave
% liviano en la mezcla vapor, seg�n el equilibrio de fases el cual predice la curva 
%y seg�n el modelo termodin�mico provisto �RM�, para obtener la ordenada y en d�nde
% una recta q de alimentaci�n de fracci�n molar de liviano zf, condici�n t�rmica q 
%y a presi�n P corta el equilibrio. 
%
%
%
%   Luis Jes�s D�az Manzo
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

