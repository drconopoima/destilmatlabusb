function [ m_factor ] = m_fun( self, omega )
%m_PREoS Evalua el m en el alfa de la ecuaci�n de estado de Peng-Robinson.
%   Luis Jes�s D�az Manzo

m_factor = 0.37464 + 1.54226.*omega - 0.26992.*omega.^2; 

end