function [ alpha ] = alfa_fun(self, T, Comp, handle_m )
%alfa_fun calcula el parámetro alfa de Peng Robinson
%opcional el input: handle_m (default PREdE.m_PREdE(w)) para el cálculo del
%factor m del factor alfa de la ecuación de PR.
%
%Referencia:
%
%Smith, Van Ness, Abbott. Introduction to Chemical Engineering Thermodynamics. 
% 7th edition 
%Storvick, Sandler. Phase Equilibria and Fluid Properties in the Chemical
%Industry: Estimation and Correlation. 1977
%   Luis Jesús Díaz Manzo

Tc = Comp.tcri;
Tr = T/Tc;
w = Comp.w_acent;

if nargin > 3
    m = handle_m(w);
else
    m = self.m_fun(w);
end
t = 1 - sqrt(Tr);
alpha = (1 + m.*(t)).^2; %Tabla 3.1 Smith Van Ness

end