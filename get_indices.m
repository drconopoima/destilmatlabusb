function database = get_indices()
%
%“get_indices.m”. En caso de modificación o expansión de las bases 
%de datos de componentes se fabrican unas nuevas bases de datos “indices.mat” e 
%“índices1.mat” al correr los algoritmos get_indices() y get_indices1(), respectivamente.
%
%“índices.mat”, “índices1.mat”: Bases de datos que relacionan el string de texto o 
%caracteres de nombre de un compuesto con el número de índice en el cual este está
% ubicado en las bases de datos, ordenados alfabéticamente. El primero recibe un 
%string y retorna un índice entre 1 y 345, el segundo recibe un índice entre 1 y 345
% y retorna un string de nombre. Se utiliza el primero para ubicar la línea de los datos 
%de los archivos contenedores de los compuestos, se utiliza el segundo para recuperar el
% string con el cual definir una clase Sustancia y mediante su método de definición 
%interno defina todas las propiedades.
%
%   Luis Jesús Díaz Manzo

fid = fopen('properties.dat');
F = textscan(fid,'%d %s %s %f %f %f %f %f %f','headerlines',3);
fclose(fid);
[id,compound,~,~,~,~,~,~,~] = F{:};
database = containers.Map();
for i=1:length(id)
    compoundi = lower(compound{i});
    database(compoundi) = id(i);
end

save('indices', 'database');

end