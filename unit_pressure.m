function [ pres_convert ] = unit_pressure( value, unit1, unit2 )
%unit_pressure convierte la presion desde la unidad unit1 hasta la unidad
%unit2, y retorna el valor de la presión en esas unidades
%   Luis Jesús Díaz Manzo
%   Se han implementado cambios de unidades entre kPa, atm, psi, psig,
%   mmH2O, 

units = {'kPa','atm', 'psi', 'mmHg', 'Pa'};
kpa = [1, 101.325, 101.325/14.696, 101.325/760, 1/1000];
inicial = 0;
final = 0;
if strcmpi(unit1, unit2) == 1
    pres_convert = value;
else
    for i = 1 : length(units)
        if strcmpi(unit1, units{i}) == 1
            inicial = i;
        elseif strcmpi(unit2, units{i}) == 1 
            final = i;
        elseif inicial && final
            break
        end 
    end
    if final == 1
        pres_convert = kpa(inicial)*value;
    end
end
end

