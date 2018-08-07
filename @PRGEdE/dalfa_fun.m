function dalpha = dalfa_fun(self, T, Comp, handle_m)
%
%
%   Luis Jesús Díaz Manzo

Tc = Comp.tcri;
Tr = T/Tc;

w = Comp.w_acent;

if nargin > 3
    m = handle_m(w);
else
    m = self.m_fun(w);
end
dalpha = ((-14.77811.*m./T-6.177251.*(m+1)./Tc).*(2.30712).^(T/Tc).*(T./Tc).^m + ( 6.177251.*2.30712.^(T./Tc))./Tc).*exp(-0.836.*T./Tc - 2).*(T./Tc).^m;

end