function salida = get_latent_heat(key_compound)
%get_latent_heat(compuesto) Obtiene la correlación calor latente del compuesto
%input compuesto debe ser un string en minúscula de nombre IUPAC sin 
%espaciados
%
%“get_latent_heat.m”. Constructor del function handle que correlaciona los calores de
% vaporización para distintas sustancias puras mediante una función de la temperatura.
%Si se utilizan múltiples bases de datos, en el archivo “heats_vaporization.dat” se 
%puede cambiar para el componente el identificador de la correlación, que por defecto 
%es “A” y agregar nuevos constructores “B” u otros.
%
%
%output es un function handle que suministra el calor de vaporización en unidades kJ/(kg-mol)
%From Perry's Chemical Engineers Handbook
%Specific heat capacities of ideal gas:
%(A) Cpig = (C1 + C2.*(((C3)./T)./(sinh((C3)./T))).^2 +C4.*(((C5./T))./(cosh(C5./T)))^2)./1000
%
% Puede agregarse otra correlación (B) que haga uso de otros datos
%
%   Luis Jesús Díaz Manzo

fid = fopen('heats_vaporization.dat');
load('indices', 'database')
key_compound(ismember(key_compound,' .:;!')) = []; %retira puntuación del nombre
key_compound = lower(key_compound); %Minusculas
try 

    indice = database(key_compound);
    F = textscan(fid,'%d %s %s %s %f %f %f %f %f %f %f %f %f %s',1,'headerlines',1+indice);
    fclose(fid);

    [~,~,~, ~, ~, C1, C2, C3 , C4, tmin, ~, Tcri, ~, ecuacion] = F{:};

    if strcmp(ecuacion,'A')
        if strcmpi(key_compound, 'CarbonMonoxide')
            salida = {@(T) C1.*1e7.*(1-T./132.92).^(C2+C3.*(T./132.92)+C4.*(T./132.92).^2).*1e-3, tmin, Tcri};
        elseif strcmpi(key_compound, 'Terephthalicacid')
            salida = {@(T) C1.*1e7.*(1-T./1113).^(C2+C3.*(T./1113)+C4.*(T./1113).^2).*1e-3, tmin, Tcri};
        else
            salida = {@(T) C1.*1e7.*(1-T./Tcri).^(C2+C3.*(T./Tcri)+C4.*(T./Tcri).^2).*1e-3, tmin, Tcri};
        end
    elseif strcmp(ecuacion, 'B')
        
    end

catch
   warning('El compuesto no está presente en la base de datos de cp de gas ideal. Puede agregarlo con funciones GET_INDICES.m y base de datos properties.dat. Y al archivo cpgasideal.dat')

end
end