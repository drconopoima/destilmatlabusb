function Hdep = entalpia(self, T, P, mezcla, fase)
% H = entalpia(T, P, mezcla, fase) Calcula la desviaci�n de entalp�a de la 
% fase requerida con respecto al gas ideal. Para la fase vapor H = 0. Para
% la fase l�quida, el resultado es el inverso del calor latente de
% vaporizaci�n
%
%Luis Jes�s D�az Manzo

if strcmpi(fase, 'vap')
    Hdep = 0;
elseif strcmpi(fase, 'liq')
    num_sust = mezcla.num_sust;
    comp = mezcla.comp;
    conc = mezcla.conc;
    Hdep = zeros(1, num_sust);
    for i = 1: num_sust
        latenth = comp(i).Hvap;
        if T < latenth{3}
            Hvap = latenth{1};
            Hdep(i) = Hvap(T).*conc(i);
        end
    end
    Hdep = -sum(Hdep);
else
    error('El valor del par�metro fase es incorrecto');
end

end