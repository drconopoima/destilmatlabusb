function salida = get_vapour_pressure(key_compound)
%get_vapour_pressure Obtiene anonymous func. de Psat @ T de key_compound
%en unidades de kPa y T en Kelvin, además obtiene anonymous func. de Tsat
%
%“get_vapour_pressure.m”. Constructor del function handle que correlaciona 
%la presión de saturación de más de 340 sustancias de la base de datos DIPPR
% en función de la temperatura. Si se utilizan múltiples bases de datos, en 
%el archivo “vapourpressure.dat” se puede cambiar para el componente el 
%identificador de la correlación, que por defecto es “A” y agregar nuevos 
%constructores “B” u otros. Es utilizado por la clase IdealEdE para el cálculo
%de la fugacidad de solución ideal.
%
%   Luis Jesús Díaz Manzo

fid = fopen('vapourpressure.dat');
load('indices', 'database')
key_compound(ismember(key_compound,' .:;!')) = []; %retira puntuación del nombre
key_compound = lower(key_compound); %Minusculas
try 

    indice = database(key_compound);
    F = textscan(fid,'%d %s %s %f %f %f %f %f %f %f %s',1,'headerlines',4+indice);
    fclose(fid);

    [~,~,~, C1, C2, C3 , C4, C5, tmin, tmax, ecuacion] = F{:};

if strcmp(ecuacion,'A')
    presion = @(T) (exp(C1 + C2./T + C3.*log(T) + C4.*T.^C5))/1000;
    temperatura = @(T,P) (C4.*T.^(C5 + 1) + C3.*T.*log(T)-(log(P*1000) - C1).*T) + C2;
elseif strcmp(formula, 'B')
    
end
pmin = presion(tmin);
pmax = presion(tmax);
salida = {presion, temperatura, tmin, tmax, pmin, pmax};
catch
   warning('El compuesto no está presente en la base de datos DIPPR. Puede agregarlo con funciones GET_INDICES.m y base de datos propiedades.dat Y al archivo vapourpressures.dat')

end
end