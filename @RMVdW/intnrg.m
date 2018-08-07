function Udep = intnrg(self, T, P, mezcla, fase)
%
%
%Luis Jesús Díaz Manzo

if (strcmpi(fase,'vap') == 1) || strcmpi(fase,'liq') == 1
	R = 8.3145;
    Comp = mezcla.comp;
    num_sust = mezcla.num_sust;
    x = mezcla.conc;
    Z = self.cmpr(T, P, mezcla);
    Tr = zeros(1, num_sust);
	for i = 1:num_sust
		Tr(i) = T/Comp(i).tcri;
	end
    if strcmpi(fase, 'liq')
        Z = min(Z);
    elseif strcmpi(fase, 'vap')
        Z = max(Z);
    end
    ai = self.ai_fun(mezcla, T);
    b_mix = self.b_fun(mezcla);
    bi = zeros(1, num_sust);
    for i = 1: num_sust;
        bi(i) = self.EdE.b_fun(Comp(i));
    end
    aij = self.aij_fun(mezcla, ai);
    am = 0;
    for i = 1:num_sust
       for j = 1:num_sust
          am = am + x(i)*x(j)*aij(i,j);
       end
    end
	alpha = zeros(1, num_sust);
	for i = 1: num_sust
		if isa(self.EdE.alfa_inp, 'function_handle')
			ninput = nargin(self.EdE.alfa_inp);
			if ninput == 1 
				alpha(i) = self.EdE.alfa_inp(T);
			elseif ninput == 2
				alpha(i) = self.EdE.alfa_inp(T, Comp(i).tcri);
			else
				if isa(self.EdE.m_inp, 'function_handle')
					m = self.EdE.m_inp;
					alpha(i) = self.EdE.alfa_inp(T, Comp(i).tcri, m(Comp(i).w_acent));
				else
					alpha(i) = self.EdE.alfa_inp(T, Comp(i));
				end
			end
		else 
			if isa(self.EdE.m_inp, 'function_handle')
				m = self.EdE.m_inp;
				alpha(i) = self.EdE.alfa_fun(T, Comp(i), m);
			else
				alpha(i) = self.EdE.alfa_fun(T, Comp(i));
			end
		end
	end
    B = b_mix*P./(R.*T);
    sigma = self.EdE.propied(1);
    epsilon = self.EdE.propied(2);
    k1 = sigma + epsilon;             
    k2 = epsilon.*sigma;                
    DELTA = (k1).^2 - 4.*(k2);
	k = mezcla.kmatrix;
	w = zeros(1, num_sust);
	for i = 1:num_sust
		w(i) = Comp(i).w_acent;
	end
	m = zeros(1, num_sust);
    for i = 1:num_sust
		if isa(self.EdE.m_inp, 'function_handle')
			m(i) = self.EdE.m_inp(w(i));
		else
			m(i) = self.EdE.m_fun(w(i));
		end
    end
    D = zeros(1, num_sust);
    dalfai = zeros(1, num_sust);
    for i = 1:num_sust
        if ~isa(self.EdE.derivada_alfa, 'function_handle')
            dalfai(i) =  self.EdE.dalfa_fun(T, Comp(i));
        elseif ~isa(self.EdE.derivada_alfa, 'function_handle') && isa(self.EdE.m_inp, 'function_handle')
            dalfai(i) = self.EdE.dalfa_fun(T, Comp(i), self.EdE.m_inp);
        elseif isa(self.EdE.derivada_alfa, 'function_handle')
            dalfai(i) = self.EdE.derivada_alfa(T, Tc(i), m(i));
        end
    end    
	for i = 1:num_sust
        for j = 1:num_sust
            D(i) = D(i) + x(i)*x(j)*(1-k(i,j))*(sqrt(ai(i))*(sqrt(ai(j))/(2.*alpha(j)))*((-dalfai(j))*T));
        end
	end
    D = 2.*D;
    D = sum(D);
    if DELTA > 0
        lg = log((2.*Z + B.*(k1 - sqrt(DELTA)))./(2.*Z + B.*(k1 + sqrt(DELTA))));
        Iv = 1./sqrt(DELTA).*lg;
    elseif DELTA == 0
        Iv = - (2.*B)./(2.*Z + k1.*B);
    else
        Iv = 2./sqrt(-DELTA).* atan((2.*Z + k1.*B)./(B.*sqrt(DELTA)));
    end
	Udep = -(((am)./b_mix).*(1+D/am).*Iv);
else
    error(['El valor "' fase '" de par?metro "fase" es incorrecto. Introduzca ''liq'' o ''vap''']);
end

end