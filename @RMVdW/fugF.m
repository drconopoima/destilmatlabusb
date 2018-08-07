function [f] = fugF(self, T, P, Mezcla, fase, element)
% [f, Z] = RMVdW.fugF(T, P, Mezcla, fase) 
%
%Esta implementación obtiene fugacidad de un conjunto de sustancias 
%por ecuaciones cúbicas
%
%Parámetros:
%Mezcla: es objeto Mezcla.m que maneja arreglos de objetos de tipo Sustancia
%
%fase: es string 'liq' para cálculo de fugacidad de un líquido, o string
%'vap' para cálculo de fugacidad de vapor
%
%    Luis Jesús Díaz Manzo

if nargin < 6 
    element = 0;
end
if (strcmpi(fase,'vap') == 1) || strcmpi(fase,'liq') == 1
    R = 8.3145;
    Comp = Mezcla.comp;
    num_sust = Mezcla.num_sust;
    x = Mezcla.conc;
    Z = self.cmpr(T, P, Mezcla);
    indice = [];
    for rre3 = 1:length(Z)
        if Z(rre3) == 0
            indice(length(indice)+1) = rre3;
        end
    end
    Z(indice) = [];
    if strcmpi(fase, 'liq')
        Z = min(Z, [], 2);
    elseif strcmpi(fase, 'vap')
        Z = max(Z, [], 2);
    end
    ai = self.ai_fun(Mezcla, T);
    b_mix = self.b_fun(Mezcla);
    bi = zeros(1, num_sust);
    for i = 1: num_sust;
        bi(i) = self.EdE.b_fun(Comp(i));
    end
    aij = self.aij_fun(Mezcla, ai);
    am = 0;
    for i = 1:num_sust
       for j = 1:num_sust
          am = am + x(i)*x(j)*aij(i,j);
       end
    end
    B = b_mix*P./(R.*T);
    A = am.*(P./(R.*T).^2);
    sigma = self.EdE.propied(1);
    epsilon = self.EdE.propied(2);
    k1 = sigma + epsilon;
    k2 = epsilon.*sigma;
    DELTA = (k1).^2 - 4.*(k2);
    if DELTA > 0
        lg = log((2.*Z + B.*(k1 - sqrt(DELTA)))./(2.*Z + B.*(k1 + sqrt(DELTA))));
        Iv = 1./sqrt(DELTA).*lg;
    elseif DELTA == 0
        Iv = - (2.*B)./(2.*Z + k1.*B);
    else
        Iv = 2./sqrt(-DELTA).* atan((2.*Z + k1.*B)./(B.*aqrt(DELTA)));
    end
    
    f = zeros(length(T), num_sust);
    for sust = 1:num_sust
        sumat = 0;
        for i = 1:num_sust
            sumat = sumat + x(i)*aij(i,sust);
        end
        f(:,sust) = P.*exp(bi(sust)./b_mix.*(Z-1) - log(Z-B)+(A./B).*(2.*sumat./am - bi(sust)./b_mix).*Iv);
    end
else
    error(['El valor "' fase '" de par�metro "fase" es incorrecto. Introduzca ''liq'' o ''vap''']);
end

if element ~= 0
    f = f(:,element);
end

end





