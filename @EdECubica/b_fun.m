function [ b_factor ] = b_fun( self, Comp )
%b_fun Calcula el covolumen de Van Der Waals. 
%   Luis Jesús Díaz
R = 8.3145;
unit_R = 'kJ/(kg-mol.K)';
OMEGA = self.propied(3);
Tc = Comp.tcri;
Pc = Comp.pcri;

b_factor = OMEGA.*(R.*Tc)./(Pc);

end

