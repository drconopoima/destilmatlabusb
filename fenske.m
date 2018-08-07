function [ nmin, recup_d, recup_b] = fenske(composic, flow,...
    sup_recup_d, sup_recup_b, light_key, heavy_key, volatildest,...
    volatilfondo)
%FENSKE Resuelve a reflujo total una columna de destilaci�n para hallar el
%n�mero m�nimo de etapas y la distribuci�n de no claves a reflujo total
%
%�Fenske.m�. [nmin, recup_d, recup_b] = fenske(composic, flow, sup_recup_d,�
%sup_recup_b, light_key, heavy_key, volatildest, volatilfondo)
%Funci�n para la resoluci�n del algoritmo de fenske para la estimaci�n del 
%m�nimo n�mero de etapas de una torre de destilaci�n en funci�n de las 
%caracter�sticas de alimentaci�n: composici�n, flujo, suposici�n trivial de
% recuperaci�n de componentes clave liviano y pesado y volatilidades 
%relativas de destilado y fondo.
%
%
%   Luis Jes�s D�az.
KD = volatildest;
KB = volatilfondo;
alfaD = KD./KD(heavy_key);
alfaB = KB./KB(heavy_key);

alfa_prom = (alfaD.*alfaB).^(0.5);
distill = sum(flow.*composic.*sup_recup_d);
residue = sum(flow.*composic.*sup_recup_b);
xD = (flow.*composic.*sup_recup_d)./(distill);
xB = (flow.*composic.*sup_recup_b)./(residue);

nmin = log(((xD(light_key))/(xD(heavy_key)))*((xB(heavy_key))/(xB(...
    light_key))))/(log(alfa_prom(light_key)));

gamma = alfa_prom.^nmin.*((xD(heavy_key)*distill)/(xB(heavy_key)*residue));

D1 = (gamma.*composic.*flow)./(1.+gamma);

B1 = composic.*flow - D1;

recup_d = D1./(flow.*composic);
recup_d(light_key) = sup_recup_d(light_key);
recup_b = B1./(flow.*composic);
recup_b(light_key) = 1 - recup_d(light_key);
recup_b(heavy_key) = sup_recup_b(heavy_key);
recup_d(heavy_key) = 1 - recup_b(heavy_key);
end