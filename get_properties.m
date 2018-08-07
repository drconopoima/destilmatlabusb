function salida = get_properties(key_compound)
%%%%%%%%%%%%%%%%%%%%%%%
% ENGLISH
%
% get_properties(compound) Get the critical properties of compoundinput compound must be a lowercase string named IUPAC without
% spaced
%
% get_properties.m. It obtains the critical properties of the DIPPR database for pure substances from the substance whose string is provided.
%
%
% Molar mass properties, Tcri(K), Pcri(kPa), Vcri (m3/kg-mol), Zc and the concentric factor (w)
%   Luis Jesús Díaz Manzo
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPANISH
%
% get_properties(compuesto) Obtiene las propiedades críticas de compuesto
% input compuesto debe ser un string en minúscula de nombre IUPAC sin
% espaciados
%
% @get_properties.m. Obtiene las propiedades críticas de la base de datos
% DIPPR para sustancias puras, de la sustancia cuyo string le sea proveido.
%
%
% Propiedades masa molar, Tcri(K), Pcri(kPa), Vcri (m3/kg-mol), Zc y el
% factor acentrico (w)
%   Luis Jesús Díaz Manzo

fid = fopen('properties.dat');
load('indices', 'database')
key_compound(ismember(key_compound,' .:;!')) = []; %retira puntuaci�n del nombre
key_compound = lower(key_compound);
try
    indice = database(key_compound);
    F = textscan(fid,'%d %s %s %f %f %f %f %f %f',1,'headerlines',2+indice);
    fclose(fid);
    [indice,~,formul,mw,tcri,pcri,~,~,w_acent] = F{:};
catch ME
    warning('El compuesto no est� presente en la base de datos. \nPuede agregarlo con funciones GET_INDICES.m y base de datos properties.dat \n')
end
salida = [formul,mw(1), tcri(1), pcri(1)*1000, w_acent(1), indice];
end
