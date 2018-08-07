function [ m_factor ] = m_fun( self, omega )
%m_fun Evalua el m en el alfa de la ecuaci�n de estado de Soave-Redlich-Kwong.
%   Luis Jes�s D�az Manzo

m_factor = 0.48 + 1.574.*omega - 0.176.*omega.^2; 

end