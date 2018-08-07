function [ exit ] = Under_theta( theta, alfa, Flow, Z, beta)
%UNDER_HANDLE Es un function handle que sirve de argumento al FSOLVE para
%hallar raíces theta que satisfagan el sistema de ecuaciones de Underwood
%   Luis Jesús Díaz Manzo
exit = zeros(1, length(theta));
for i = 1:length(theta)
    aux = zeros(1, length(Z));
    for j=1 : length(alfa)
        aux(j) = ((alfa(j)*Z(j)*Flow)/(alfa(j)-theta(i)));
    end
    
    exit(i) = sum(aux) - beta*Flow;
end
end