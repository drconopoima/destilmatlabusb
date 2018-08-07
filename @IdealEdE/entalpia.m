function Hdep = entalpia(self, T, P, mezcla, fase)
% H = entalpia(T, P, mezcla, fase) Calcula la desviación de entalpía de la 
% fase requerida con respecto al gas ideal. Para la fase vapor H = 0. Para
% la fase líquida, el resultado es el inverso del calor latente de
% vaporización
%
%Luis Jesús Díaz Manzo

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
    error('El valor del parámetro fase es incorrecto');
end

end