function [f] = fugF(self, T, P, mezcla, fase, element)
%
%IdealEdE.fugF(T, P, Comp, fase) 
%
%Esta implementación obtiene fugacidad por solución ideal para fase = 'liq'
%igual a Pvap/P para líquidos y retorna f = 1 para vapor (gas ideal)
%
%Parámetros:
%Comp es objeto Sustancia.m que maneja propiedades de sustancias puras 
%también puede ser un objeto de tipo Sustancia.m
%
%fase: es string 'liq' para cálculo de fugacidad de un líquido, o string
%'vap' para cálculo de fugacidad de vapor
%
%    Luis Jesús Díaz Manzo
Comp = mezcla.comp;
num_sust = mezcla.num_sust;
f = ones(1, num_sust);
if nargin > 5 && ~isempty(element)
    element = 0;
end
if strcmp(fase,'vap') == 1
elseif strcmp(fase,'liq') == 1
    for i = 1:num_sust
        Psat = Comp(i).psat{1};
        Psat = Psat(T);
        f(i) = Psat/P;
    end
else
    error(['El valor "' fase '" de parámetro "fase" es incorrecto.']);
end

if element ~= 0
    f = f(element);
end

f = abs(f);

end