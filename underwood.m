function [ rmin, dflow, incog ] = underwood(composic, flow, sup_recup_d, sup_recup_b,...
    KD, light_key, heavy_key, beta)
%underwood Con el método corto shorcut de dimensionamiento obtiene reflujo
%mínimo
%
%“underwood.m”. [rmin, dflow, incog] = underwood(composic, flow, sup_recup_d,…
%                                   sup_recup_b, KD, light_key, heavy_key, beta)
%
%Función que realiza la distribución a reflujo mínimo de los productos de 
%destilado y fondo así como calcula el valor del reflujo mínimo (rmin), del
% flujo de destilado (dflow) y la distribución de compuestos no claves (incog) 
%de longitud arbitraria. Recibe las composiciones, flujos, una suposición de
% recuperación en destilado y residuo de los componentes clave liviano y pesado,
% un vector de coeficientes de distribución líquido vapor y una fracción vaporizada
% de la alimentación.
%
%
%   Luis Jesús Díaz Manzo
num_comp = numel(composic);
alfa = KD./KD(heavy_key);

for i = 1:num_comp
    if sup_recup_d(i) >= 1
        sup_recup_d(i) = 1 - 1e-15;
        sup_recup_b = 1 - sup_recup_d;
    elseif sup_recup_d(i) <= 0
        sup_recup_d(i) = 1e-15;        
        sup_recup_b = 1 - sup_recup_d;
    end
end

D = sum(flow.*composic.*sup_recup_d);
xD = (flow.*composic.*sup_recup_d)./(D);
dif_alfa_j = (alfa-1);
dif_alfa_lk = (alfa(light_key) - 1);
dif_alfa_lk_j = (alfa(light_key) - alfa);
corroborar_claves = (dif_alfa_j./dif_alfa_lk).*((xD(...
    light_key).*D)./(composic(light_key).*flow))+((dif_alfa_lk_j)./(...
    dif_alfa_lk)).*((xD(heavy_key).*D)./(composic(heavy_key).*flow));
claves = [];
for i=1:light_key
    if corroborar_claves(i) > 0.0000001 && corroborar_claves(i) < ... 
            0.9999999
        claves(end+1) = i;
        if corroborar_claves(i) < sup_recup_d(i)
            sup_recup_d(i) = corroborar_claves(i);
            sup_recup_b = 1 - sup_recup_d;
        end
    end
end
for i = light_key+1:heavy_key
    claves(end + 1) = i;
end
for i = heavy_key:num_comp
    if corroborar_claves(i) > 0.0000001 && corroborar_claves(i) < ... 
            0.9999999
        if ~any(claves == i)
            claves(end+1) = i;      
        end
        if corroborar_claves(i) < sup_recup_d(i)
            sup_recup_d(i) = corroborar_claves(i);
            sup_recup_b = 1 - sup_recup_d;
        end
    end
end
if isempty(claves) || length(claves) == 1
    claves = [light_key, heavy_key];
end
D = sum(flow.*composic.*sup_recup_d);
xD = (flow.*composic.*sup_recup_d)./(D);
num_thetas = length(claves)-1;
alfa_semilla = zeros(1, num_thetas);
options = optimset('Display', 'none', 'TolX', 1e-7, 'TolFun', 1e-7);
if num_thetas == 1    
    theta = fsolve(@(x) Under_theta(x, alfa, flow, composic, beta),...
        ((alfa(light_key)+alfa(heavy_key))/2), options);    
    
    Rmin = sum((alfa.*xD)./(alfa-theta))-1;
    dist_flow = sum(sup_recup_d.*flow.*composic);
    incog = [Rmin, dist_flow];
else
    for i=min(claves):1:max(claves)-1
        alfa_semilla(i + 1 - min(claves)) = (alfa(i) + alfa(i + 1))/2;
    end
    theta = fsolve(@(x) Under_theta(x, alfa, flow, composic, beta),...
        alfa_semilla, options);
    indice_lk = find((claves - light_key) == 0);
    
    indice_hk = find((claves - heavy_key) == 0);
    try 
        claves(indice_lk) = [];
    catch
    end
    try
        claves(indice_hk-1) = [];
    catch
    end
    for i= 1:length(theta) + 1
        if i == length(theta)
            semilla_incog(i) = theta(1);
        elseif i == length(theta) + 1
            semilla_incog(i) = D;
        else
            semilla_incog(i) = xD(claves(i));
        end
    end
    incog = fsolve(@(x) Under_Riguroso(x, theta, alfa, sup_recup_d,...
        composic, flow, claves),semilla_incog, options);
    if abs((semilla_incog(end-1) - incog(end-1)))> 0.5
        semilla_incog(end-1) = incog(end-1);
        incog = fsolve(@(x) Under_Riguroso(x, theta, alfa, sup_recup_d,...
            composic, flow, claves),semilla_incog, options);
    end    
end
dflow = incog(end);
    rmin = incog(end-1);
    incog1 = zeros(1, length(incog) + length(claves));
    incog1(1:length(claves)) = claves;
    incog1(length(claves)+1:end) = incog;
    incog = incog1;
end