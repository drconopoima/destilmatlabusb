function sumat = RachfordRice(beta,z,K)
%“RachfordRice.m”. RachfordRice(beta, z, K).
%Función objetivo para obtención de la fracción vaporizada por el método
% de Rachford Rice para su uso en un FSOLVE().


sumat = (z.*(1-K))./(1+beta.*(K-1));
sumat = sum(sumat);