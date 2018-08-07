function [ vi ] = volumesp( self, T, P, mezcla, fase)
%volumesp Calcula el volumen específico de cada componente en mezcla
%   Luis Jesús Díaz Manzo

num_sust = mezcla.num_sust;
comp = mezcla.comp;
R = 8.3145;
Zi = zeros(1, num_sust);
vi = zeros(1, num_sust);
for i=1: num_sust
    if strcmpi(fase, 'vap')
        Zi(i) = max(self.EdE.cmpr(T, P, comp(i)));
    elseif strcmpi(fase, 'liq')
        Zi(i) = min(self.EdE.cmpr(T, P, comp(i)));
    end
    vi(i) = R.*T.*Zi(i)./(P);
end

end

