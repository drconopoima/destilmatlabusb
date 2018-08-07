function [ m_factor ] = m_fun( self, omega )
%m_PREoS Evalua el m en el alfa de la ecuación de estado de Peng-Robinson
%según modificación hecha por Gasem en 2001
%Referencia:
%Gasem, Gao, Pan & Robinson, Fluid Phase Equilibria, 181, 113-125 (2001).
%   Luis Jesús Díaz Manzo

m_factor = 0.134 + 0.508.*omega - 0.0467.*omega.^2;

end

