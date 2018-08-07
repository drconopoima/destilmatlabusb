function [f] = fug(self, T, P, Mezcla, fase)
%f = IdealEdE.fug(T, P, Mezcla, fase) 
%
%Esta implementación obtiene fugacidad 
%
%Parámetros:
%Mezcla: es objeto Mezcla.m que maneja arreglos de objetos de tipo Sustancia
%
%También puede Mezcla ser un array de clases Sustancia.m para compuestos puros
%
%fase: es string 'liq' para cálculo de fugacidad de un líquido, o string
%'vap' para cálculo de fugacidad de vapor
%
%Referencia: 
%
%Perry, R., & Green, D. (1999). Perry's Chemical Engineering Handbook. 
%McGraw-Hill Companies. Capítulo 2. 'Physical Properties'
%
%   Luis Jesús Díaz
if isa(Mezcla, 'Mezcla')
    if (strcmp(fase,'vap') == 1) || (strcmp(fase,'liq') == 1)
        f = self.fugF(T, P, Mezcla, fase);
    else
        error(['El valor "' fase '" del parámetro ''fase'' es incorrecto'])
    end
end
if isa(Mezcla, 'Sustancia')
    num_sust = length(Mezcla);
    if (strcmp(fase,'vap') == 1) || (strcmp(fase,'liq') == 1)
        for i = 1: num_sust
            f(i) = self.fugF(T, P, Mezcla(i), fase);
        end
    else
        error(['El valor "' fase '" del parámetro ''fase'' es incorrecto'])
    end
end
end
