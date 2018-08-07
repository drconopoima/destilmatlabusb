function [dalpha] = dalfa_fun(self, T, Comp, handle_m)
		
	Tc = Comp.tcri;
	Tr = T/Tc;
	w = Comp.w_acent;

	if nargin > 3
		m = handle_m(w);
	else
		m = self.m_fun(w);
	end
	dalpha = m.*(m.*sqrt(T/Tc) - m - 1)/(Tc.*sqrt(T./Tc));
	
end