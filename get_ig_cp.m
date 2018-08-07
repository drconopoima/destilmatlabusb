function salida = get_ig_cp(key_compound)
%get_ig_cp(compuesto) Obtiene la correlación calor específico del compuesto
%input compuesto debe ser un string en minúscula de nombre IUPAC sin 
%espaciados
%
%output es un function handle que suministra el cp en unidades kJ/(kg-mol.K)
%From Perry's Chemical Engineers Handbook
%Specific heat capacities of ideal gas:
%(A) Cpig = (C1 + C2.*(((C3)./T)./(sinh((C3)./T))).^2 +C4.*(((C5./T))./(cosh(C5./T)))^2)./1000
%
% Puede agregarse otra correlación (B) que haga uso de otros datos
%
%   Luis Jesús Díaz Manzo

fid = fopen('cpgasideal.dat');
load('indices', 'database')
key_compound(ismember(key_compound,' .:;!')) = []; %retira puntuación del nombre
key_compound = lower(key_compound); %Minusculas
try 

    indice = database(key_compound);
    F = textscan(fid,'%d %s %s %s %f %f %f %f %f %f %f %f %f %f %s',1,'headerlines',2+indice);
    fclose(fid);

    [~,~,~, ~, ~, C1, C2, C3 , C4, C5, tmin, ~, tmax, ~, ecuacion] = F{:};

    if strcmp(ecuacion,'A')
        salida = {@(T) (C1.*1e5 + C2.*1e5.*(((C3.*1e3)./T)./(sinh((C3.*1e3)./T))).^2 +C4.*1e5.*(((C5./T))./(cosh(C5./T))).^2).*1e-3, tmin, tmax};
    elseif strcmp(ecuacion, 'B')

    end

catch
   warning('El compuesto no está presente en la base de datos de cp de gas ideal. Puede agregarlo con funciones GET_INDICES.m y base de datos properties.dat. Y al archivo cpgasideal.dat')

end
end