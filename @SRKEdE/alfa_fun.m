function [ alpha ] = alfa_fun(self, T, Comp, handle_m )
%alfa_fun calcula el par�metro alfa de Soave-Redlich-Kwong
%opcional el input: handle_m (default SRKEdE.m_SRKEdE(w)) para el c�lculo del
%factor m del factor alfa de la ecuaci�n de SRK.
%
%Referencia:
%
%Smith, Van Ness, Abbott. Introduction to Chemical Engineering Thermodynamics. 
% 7th edition 
%Storvick, Sandler. Phase Equilibria and Fluid Properties in the Chemical
%Industry: Estimation and Correlation. 1977
%   Luis Jes�s D�az Manzo

Tc = Comp.tcri;
Tr = T/Tc;
w = Comp.w_acent;

if nargin > 3
    m = handle_m(w);
else
    m = self.m_fun(w);
end
t = 1 - sqrt(Tr);
alpha = (1 + m*(t))^2; %Tabla 3.1 Smith Van Ness

end