function [ m_factor ] = m_fun( self, omega )
%m_PREoS Evalua el m en el alfa de la ecuaci�n de estado de Peng-Robinson
%seg�n modificaci�n hecha por Gasem en 2001
%Referencia:
%Gasem, Gao, Pan & Robinson, Fluid Phase Equilibria, 181, 113-125 (2001).
%   Luis Jes�s D�az Manzo

m_factor = 0.134 + 0.508.*omega - 0.0467.*omega.^2;

end

