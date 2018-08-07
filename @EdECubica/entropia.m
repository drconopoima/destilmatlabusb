function Sdep = entropia(self, T, P, Comp, fase)
% Sdep = EdECubica.entalpia(T, P, Comp, fase)
%
%Calcula la entropia residual a temperatura T, presión P, del compuesto Comp
%en la fase 'liq' o 'vap'
%
%Referencia:
%
%Aungier, R.H. "A Fast, Accurate Real Gas Equation of State for Fluid Dynamic Analysis Applications". Journal of Fluids Engineering. 117. 277–281. 1995.
%
%Figueira, Freddy (2005). "Desarrollo de ecuaciones de estado del tipo Van
%der Waals para fluidos puros polares y no polares". Universidad Simón
%Bolívar. Venezuela
%
%	Luis Jesús Díaz Manzo

Z = self.cmpr(T, P, Comp);
if strcmp(fase, 'liq')
    Z = min(Z);
elseif strcmp(fase, 'vap')
    Z = max(Z);
else
    error(['El valor "' fase '" del parámetro "fase" es incorrecto. Introduzca ''liq'' o ''vap'''])
end
R = 8.3145;
unit_R = 'kJ/(kg-mol.K)';
Tc = Comp.tcri;
Pc = Comp.pcri;
omega = Comp.w_acent;
B = self.b_fun(Comp).*P./(R.*T);
b = self.b_fun(Comp);
sigma = self.propied(1);
epsilon = self.propied(2);
k1 = sigma + epsilon;
k2 = epsilon.*sigma;
DELTA = (k1).^2 - 4.*(k2);
if DELTA > 0
    lg = log((2.*Z + B.*(k1 - sqrt(DELTA)))./(2.*Z + B.*(k1 + sqrt(DELTA))));
    Iv = 1./sqrt(DELTA).*lg;
elseif DELTA == 0
    Iv = - (2.*B)./(2.*Z + k1.*B);
else
    Iv = 2./sqrt(-DELTA).* atan((2.*Z + k1.*B)./(B.*sqrt(DELTA)));
end
if isa(self.derivada_alfa, 'function_handle')
    dalfa_dt = self.derivada_alfa;
	ninput = nargin(dalfa_dt);
	if ninput == 1
		dalfa_dt = self.derivada_alfa(T);
	elseif ninput == 2
		dalfa_dt = dalfa_dt(T, Tc);
	elseif ninput == 3
		if isa(self.m_inp, 'function_handle')             
			m = self.m_inp;
			if ~isempty(m)
				if isa(m, 'function_handle')
					minput = nargin(m);
				else
					minput = 0;
				end
				if minput > 0
					m = m(omega);
				end
			else 
				m = self.m_fun(omega);
			end
			dalfa_dt = dalfa_dt(T, Tc, m);
		end
	end
else
	dalfa_dt = self.dalfa_fun(T, Comp);
end

PSI = self.propied(4);
da_dt = PSI.*(dalfa_dt.*R.^2.*Tc.^2)./(Pc);

Sdep = R.*log(Z-B) - ((da_dt)./b).*Iv;

end