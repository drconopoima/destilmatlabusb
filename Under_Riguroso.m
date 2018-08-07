function [ exit ] = Under_Riguroso( soluc, theta, alfa, Recup_D,...
    Z, Flow, claves)
%“Under_Riguroso”. Under_Riguroso(soluc, theta, alfa, Recup_D, Z, Flow, claves).
%
%Función objetivo para la resolución con fsolve() del sistema de ecuaciones 
%de una columna de destilación simple en condición de reflujo mínimo por el 
%algoritmo de Underwood. 
%
%   Luis Jesús Díaz Manzo

exit = zeros(1, length(theta) + 1);
dist_flow = 0;
for i = 1: length(alfa)
    if ~(any(claves == i))
        dist_flow = dist_flow + Recup_D(i).*Flow.*Z(i);
    end
end
for l = 1 : length(exit)-2
    dist_flow = dist_flow + soluc(l)*soluc(end);
end
for i = 1:length(theta)
    indice_soluc = 1;
    aux = sum((alfa.*Recup_D.*Z.*Flow)./(alfa - theta(i)));
    for j = 1:length(alfa)
        if (any(claves == j))
            aux = aux - ((alfa(j).*Recup_D(j).*Z(j).*Flow)./(alfa(...
                j) - (theta(i)))) + ((alfa(j).*(soluc(indice_soluc)...
                )*soluc(end))/((alfa(j) - (theta(i)))));
            indice_soluc = indice_soluc+1;
        end
    end
    exit(i) = sum(aux) - soluc(end)*(1 + soluc(end - 1));
end
exit(end) = dist_flow - soluc(end);
end

