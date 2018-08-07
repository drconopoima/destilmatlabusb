classdef EdECubica
    %EdECubica Clase que maneja todas las ecuaciones de estado cúbicas
    %
    %
    %Clase: EdECubica  < handle
    %Esta clase implementa el cálculo de las propiedades termodinámicas de sustancias
    %       puras. La clase provee métodos para calcular todas esas propiedades, que 
    %       son todas heredadas por las clases EdE derivadas de ellas
    %
    %
    %presion(self, T, v, comp, parametros, alfa_func). Calcula la presión a partir
    %    de la temperatura y el volumen específico para una clase sustancia en ‘comp’,
    %    utilizando parámetros de la ecuación de estado cúbica generalizada, 
    %    ecuaciones (3), (4) y (5) y el function handle en ‘alfa_func’ para el cálculo
    %    del parámetro de atracción de la ecuación de estado 
    %
    %
    %
    %
    %Referencia:
    %
    %Smith, Van Ness, Abbott. Introduction to Chemical Engineering Thermodynamics. 
    % 7th edition 
    %
    %Poling, B., Prausnitz, J., O'Connell, J. (2001)"The Properties of 
    %Gases And Liquids". Quinta Edición. McGraw-Hill
    %
    %   Luis Jesús Díaz Manzo
    
    properties
        familia = 'Cubicas';
    end
    
    methods
        function self = EdECubica()
        end
        function pcalc = presion(self, T, v, Comp, parametros, alfa_func) 
            %alfa_func function handle que reemplaza al default si es
            %suministrado
            R = 8.314;
            b = self.b_fun(Comp);
            if nargin > 4              
                sigma = parametros(1);
                epsilon = parametros(2);
            else
                sigma = self.propied(1);
                epsilon = self.propied(2);
            end
            if nargin > 5
                a = self.a_fun(T, Comp, alfa_func);
            else
                a = self.a_fun(T, Comp);
            end
            pcalc = R.*T./(v - b) - (a)./((v + epsilon.*b).*(v + sigma.*b));
        end
        function [vliq, vvap] = vsatT(self, T, Comp)
            
        end
        function [vliq, vvap] = vsatP(self, P, Comp)
            
        end
        function Z = cmpr(self, T, P, Comp) %Calcula el factor de compresibilidad
            %Las letras griegas se refieren a la nomenclatura utilizada por
            %Poling, Prausnitz, O'Connell en la referencia
            R = 8.3145;
            b = self.b_fun(Comp);
            delta = 2.*b;
            epsilon = -b.^2;
            THETA = self.a_fun(T, Comp);
            Bprima = b.*P./(R.*T);
            deltaprima = delta.*(P)./(R.*T);
            THETAprima = THETA.*P./(R.*T).^2;
            epsilonprima = epsilon.*(P./(R.*T)).^2;
            a0 = 1;
            a1 = (deltaprima - Bprima - 1);
            a2 = (THETAprima + epsilonprima - deltaprima.*(Bprima + 1));
            a3 = -(epsilonprima.*(Bprima + 1) + THETAprima.*(Bprima));
            Z = roots([a0, a1, a2, a3]);
            Z = Z(imag(Z)==0); % Solo raices reales
            Z = Z(Z > 0); %Solo raices positivas
        end
        function [Tsat, error, ite] = isofugT(self, P, Comp, varargin)
            %Por igualación de fugacidad de componente puro halla la Psat
            Ts = Comp.tsat{1};
            Tmin = Comp.psat{2};
            Tmax = Comp.psat{3};
            Pmin = Comp.tsat{2};
            Pmax = Comp.tsat{3};
            if nargin < 4 
                tol = 1e-4;
            else
                tol = varargin{1};
                if isempty(tol)
                    tol = 1e-4;
                end
            end
            if nargin < 5
                maxite = 100;
            else
                maxite = varargin{2};
                if isempty(maxite)
                    maxite = 100;
                end
            end
            ite = 0;       
            Tsemilla = (Tmax - Tmin)./(Pmax - Pmin).*(P - Pmin) + Tmin;
            options = optimset('Display', 'off');
            Tsupi = fzero(@(T) Ts(T, P), Tsemilla, options);
            dT = Tsupi.*0.03;         
            if 1e-1 > tol
                error = 1e-1;
            else
                error = abs(2*tol);
            end
            while abs(error)> tol
                ite = ite + 1;
                fugL = self.fugF(Tsupi, P, Comp, 'liq');
                fugV = self.fugF(Tsupi, P, Comp, 'vap');
                error = (fugV - fugL)/P;            
                if abs(error) > tol
                    if ite > 2 
                        if sign(error) ~= sign(errorant)
                            Tsup = (Tsupi - Tsupiant)./(error-errorant).*(-errorant)+ Tsupiant;
                            dT = dT ./ 5;
                        elseif sign(error) ~= sign(errorantant)
                            Tsup = (Tsupi - Tsupiantant)./(error-errorantant).*(-errorantant)+ Tsupiantant;
                            dT = dT ./ 5;
                        else
                            if fugV > fugL
                                Tsup = Tsupi + dT;
                            else                        
                                Tsup = Tsupi - dT;
                            end   
                        end
                        
                    elseif ite > 1                        
                        if sign(error) ~= sign(errorant)
                            Tsup = (Tsupi - Tsupiant)./(error - errorant).*(-errorant)+ Tsupiant;
                            dT = dT ./ 5;
                        else
                            if fugV > fugL
                                Tsup = Tsupi + dT;
                            else                        
                                Tsup = Tsupi - dT;
                            end 
                        end
                    else
                        if fugV > fugL
                            Tsup = Tsupi + dT;
                        else                        
                            Tsup = Tsupi - dT;
                        end
                    end
                    
                    if ite >= 2
                        Tsupiantant = Tsupiant;
                        errorantant = errorant;
                    end
                    errorant = error;
                    Tsupiant = Tsupi;
                    Tsupi = Tsup;
                end
                if ite == maxite
                    warning('No se ha logrado con el numero de iteraciones un resultado')
                    break
                end
            end
            Tsat = Tsupi;
        end
        function [Psat, error, ite] = isofugP(self, T, Comp, varargin)
            %Por igualación de fugacidad de componente puro halla la Psat
            R = 8.314;
            Ps = Comp.psat{1};
            Psupi = Ps(T);
            dP = Psupi.*0.1;
            if nargin < 4 
                tol = 1e-3;
            else
                tol = varargin{1};
                if isempty(tol)
                    tol = 1e-3;
                end
            end
            if nargin < 5
                maxite = 100;
            else
                maxite = varargin{2};
                if isempty(maxite)
                    maxite = 100;
                end
            end
            error = 100000;
            ite = 0;
            while (error)> tol
                ite = ite + 1;
                fugL = self.fugF(T, Psupi, Comp, 'liq');
                fugV = self.fugF(T, Psupi, Comp, 'vap');
                error = abs(fugV - fugL)/Psupi;
                if fugV > fugL
                    Psupi = Psupi - dP;
                else                        
                    Psupi = Psupi + dP;
                end                    
                if dP > 0.5*error * Psupi;
                    dP = dP ./ (3);
                end
                if ite == maxite
                    warning('No se ha logrado con el numero de iteraciones un resultado')
                    break
                end
            end
            Psat = Psupi;
        end
    end
    
end

