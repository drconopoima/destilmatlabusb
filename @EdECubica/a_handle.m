function [ a_factor ] = a_handle(self, T, Comp, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%
% English Explanation
% a_factor Calcula el parámetro de atracción a de la EdE cúbica correspondiente
% opcional el input: handle_alfa (default PREdE.alfa_fun(T, Comp)) para el cálculo del factor alfa de la ecuación de PR.
%
% Referencias:
%
% Smith, Van Ness, Abbott. Introducción a la Ingeniería Química Termodinámica.
% 7ª edición
% Storvick, Sandler. Equilibrios de fase y propiedades de los fluidos en la sustancia química
% Industria: Estimación y Correlación. 1977
%   Luis Jesús Díaz Manzo
%%%%%%%%%%%%%%%%%%%%%%%%
% EXPLICACIÓN EN ESPAÑOL
%
% a_factor Calcula el parámetro de atracción a de la EdE cúbica correspondiente
% opcional el input: handle_alfa (default PREdE.alfa_fun(T, Comp)) para el cálculo del factor alfa de la ecuación de PR.
%
% Referencias:
%
% Smith, Van Ness, Abbott. Introduction to Chemical Engineering Thermodynamics.
% 7th edition
% Storvick, Sandler. Phase Equilibria and Fluid Properties in the Chemical
% Industry: Estimation and Correlation. 1977
%   Luis Jesús Díaz Manzo

PSI = self.propied(4);
R = 8.3145;
unit_R = 'kJ/(kg-mol.K)';
Tc = Comp.tcri;
Pc = Comp.pcri;
omega = Comp.w_acent;
if nargin > 3
    alfa = varargin{1};
    if ~isempty(alfa)
        if isa(alfa, 'function_handle')
            ninput = nargin(alfa);
        else
            ninput = 0;
        end
        if ninput == 1
            alfa = alfa(T);
        elseif ninput == 2
            alfa = alfa(T, Tc);
        elseif ninput == 3
            if nargin > 4
                m = varargin{2};
                if ~isempty(m)
                    if isa(m, 'function_handle')
                        minput = nargin(m);
                    else
                        minput = 0;
                    end
                    if minput > 0
                        m = m(omega);
                    end
                else
                    m = self.m_fun(omega);
                end
                alfa = alfa(T, Tc, m);
            end
        end
    elseif nargin > 4
        m = varargin{2};
        if ~isempty(m)
            if isa(m, 'function_handle')
                minput = nargin(m);
            else
                minput = 0;
            end
            if minput == 0
                m = @(w) m;
            end
            alfa = self.alfa_fun(T, Comp, m);
        end
    end
else
    alfa = self.alfa_fun(T, Comp);
end
a_factor = PSI.*(alfa.*R.^2.*Tc.^2)./(Pc);
end
