function [f, Z] = fugF(self, T, P, Comp, fase, element)
%[f, Z] = fugF(T, P, Comp, fase) calcula la fugacidad de un fluido a la temperatura
%T, la presión P y fase puede ser 'liq' para líquido y 'vap' para vapor.
%
%Referencias
%
%Figueira, Freddy (2005). "Desarrollo de ecuaciones de estado del tipo Van
%der Waals para fluidos puros polares y no polares". Universidad Simón
%Bolívar. Venezuela
%
%   Luis Jesús Díaz Manzo

Z = self.cmpr(T, P, Comp);
if nargin > 5 && ~isempty(element)
    element = 0;
end
if strcmp(fase, 'liq')
    Z = min(Z);
elseif strcmp(fase, 'vap')
    Z = max(Z);
else
    error(['El valor "' fase '" del parámetro "fase" es incorrecto. Introduzca ''liq'' o ''vap'''])
end
R = 8.3145;
B = self.b_fun(Comp).*P./(R.*T);
A = self.a_fun(T, Comp).*(P./(R.*T)^2);
sigma = self.propied(1);
epsilon = self.propied(2);
k1 = sigma + epsilon;
k2 = epsilon.*sigma;
DELTA = (k1).^2 - 4.*(k2);
if DELTA > 0
    lg = log((2.*Z + B.*(k1 - sqrt(DELTA)))./(2.*Z + B.*(k1 + sqrt(DELTA))));
    Iv = 1./sqrt(DELTA).*lg;
elseif DELTA == 0
    Iv = - (2.*B)./(2.*Z + k1.*B);
else
    Iv = 2./sqrt(-DELTA).* atan((2.*Z + k1.*B)./(B.*sqrt(DELTA)));
end
f = P.*exp(Z-1-log(Z-B)+(A./B).*Iv);

if element ~= 0
    f = f(element);
end
f = abs(f);
end