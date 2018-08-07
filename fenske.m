function [ nmin, recup_d, recup_b] = fenske(composic, flow,...
    sup_recup_d, sup_recup_b, light_key, heavy_key, volatildest,...
    volatilfondo)
%FENSKE Resuelve a reflujo total una columna de destilación para hallar el
%número mínimo de etapas y la distribución de no claves a reflujo total
%
%“Fenske.m”. [nmin, recup_d, recup_b] = fenske(composic, flow, sup_recup_d,…
%sup_recup_b, light_key, heavy_key, volatildest, volatilfondo)
%Función para la resolución del algoritmo de fenske para la estimación del 
%mínimo número de etapas de una torre de destilación en función de las 
%características de alimentación: composición, flujo, suposición trivial de
% recuperación de componentes clave liviano y pesado y volatilidades 
%relativas de destilado y fondo.
%
%
%   Luis Jesús Díaz.
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