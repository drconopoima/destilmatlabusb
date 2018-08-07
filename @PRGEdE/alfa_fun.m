function [ alpha ] = alfa_fun(self, T, Comp, handle_m )
%alfa_fun calcula el par�metro alfa de Peng Robinson Gasem
%opcional el input: handle_m (default PRGEdE.m_PRGEdE(w)) para el c�lculo del
%factor m del factor alfa de la ecuaci�n de PRG.
%Referencia:
%Gasem, Gao, Pan & Robinson, Fluid Phase Equilibria, 181, 113-125 (2001).
%La expresi�n para alfa en esta EdE:
% t = 1 - (Tr).^m;
% m = 0.134 + 0.508.*w - 0.0467.*w.^2
% alfa = exp[(A + B.*Tr).*t]
% A = 2.00
% B = 0.836
%
%   Luis Jes�s D�az Manzo

Tc = Comp.tcri;
Tr = T/Tc;

w = Comp.w_acent;

if nargin > 3
    m = handle_m(w);
else
    m = self.m_fun(w);
end
t =  1 - (Tr).^m;
A = 2.00;
B = 0.836;
alpha = exp((A + B.*Tr).*t); %Referencia de introducci�n
end