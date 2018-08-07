classdef RMVdW < IdealEdE
    %MezcladoCubicas Clase maneja reglas de mezclado de EdE cúbicas
    %mediante las reglas de mezclado de Van der Waals.
    %
    %Referencia:
    %
    %Kwak, T., Mansoori, G. (1985). "Van der Waals Mixing Rules for Cubic Equations
    %of State. Applications for Supercritical Fluid Extraction Modelling". Ch. Eng. 
    %Sci, 41, no. 5, 1303-1309.
    %
    %Smith, Van Ness, Abbott. Introduction to Chemical Engineering Thermodynamics. 
    % 7th edition 
    %
    %Poling, B., Prausnitz, J., O'Connell, J. (2001)"The Properties of 
    %Gases And Liquids". Quinta Edición. McGraw-Hill
    %
    %   Luis Jesús Díaz Manzo
    properties
        familia = 'Cubicas';
        EdE = false;
    end    

    methods
        
        function self = RMVdW(EcdE)
            if nargin > 0
                self.EdE = EcdE;
            end
        end
        function Z = cmpr(self, T, P, mezcla)
            %Las letras griegas se refieren a la nomenclatura utilizada
            % por Poling, Prausnitz, O'Connell en la referencia
            %
            %Las reglas de mezclado de Van der Waals tomadas de Kwak & Mansoori
            if isa(mezcla, 'Mezcla')
                R = 8.3145;
                num_sust = mezcla.num_sust;
                Comp = mezcla.comp;
                x = mezcla.conc;
                ai = self.ai_fun(mezcla, T);
                bi = zeros(1, num_sust);
                for i = 1: num_sust;
                    bi(i) = self.EdE.b_fun(Comp(i));
                end
                [aij] = self.aij_fun(mezcla, ai);
                am = 0;
                bm = self.b_fun(mezcla);
                for i = 1:num_sust
                   for j = 1:num_sust
                      am = am + x(i)*x(j)*aij(i,j);
                   end
                end
                delta = 2.*bm;
                epsilon = -bm.^2;
                THETA = am;
                Bprima = bm.*P./(R.*T);
                deltaprima = delta.*(P)./(R.*T);
                THETAprima = THETA.*P./(R.*T).^2;
                epsilonprima = epsilon.*(P./(R.*T)).^2;
                a1 = (deltaprima - Bprima - 1);
                a0 = ones(length(a1), 1);
                a2 = (THETAprima + epsilonprima - deltaprima.*(Bprima + 1));
                a3 = -(epsilonprima.*(Bprima + 1) + THETAprima.*(Bprima));
                Z = zeros(length(a0), 3);
                lengt = size(Z);
                lengt = lengt(1); 
                for i = 1:lengt
                    Z(i,:) = roots([a0(i), a1(i), a2(i), a3(i)]);
                end
                Z(imag(Z)~=0)  = 0; % Solo raices reales
                Z(Z < 0) = 0; %Solo raices positivas
            else
                error('Argumento ''Mezcla'' no es de clase ''Mezcla''')
            end
        end
        function b_factor = b_fun(self, Mezcla)
            %b_fun Calcula el covolumen de Van Der Waals. 
            %   Luis Jesús Díaz
            unit_R = 'kJ/(kg-mol.K)';            
            comp = Mezcla.comp;
            x = Mezcla.conc;
            num_sust = Mezcla.num_sust;
            b = zeros(1, num_sust);
            for i = 1: num_sust;
                b(i) = self.EdE.b_fun(comp(i));
            end
            b_factor = sum(b.*x);
        end
        function a_factor = a_fun(self, Mezcla, T)
            comp = Mezcla.comp;
            x = Mezcla.conc;
            num_sust = Mezcla.num_sust;
            kij = Mezcla.kij;
            THETA = self.ai_fun(Mezcla, T);
            THETA_mix = sum(x.^2.*THETA); 
            tamano = size(kij);
            o = 0;
            while o < tamano(2)
                for q = o: num_sust
                    for k = q+1 : num_sust       
                        o = o+1;
                        THETA_mix = THETA_mix + 2.*x(q).*x(k).*sqrt(THETA(q).*THETA(k)).*(1 - kij(o));
                    end
                end
            end
            a_factor = THETA_mix;
        end
        function a_ipuro = ai_fun(self, Mezcla, T)
            num_sust = Mezcla.num_sust;
            comp = Mezcla.comp;
            a_ipuro = zeros(length(T), num_sust);
            for j = 1: num_sust
                a_ipuro(:, j) = self.EdE.a_fun(T, comp(j));
            end
        end
		function dalfa_ipuro = dalfai_fun(self, Mezcla, T)
            num_sust = Mezcla.num_sust;
            comp = Mezcla.comp;
            dalfa_ipuro = zeros(1, num_sust);
			for j = 1: num_sust
				if isa(self.EdE.derivada_alfa, 'function_handle')
					dalfa_puro = self.EdE.derivada_alfa;
					ninput = nargin(dalfa_dt);
					if ninput == 1
						dalfa_ipuro(j) = dalfa_puro(T);
					elseif ninput == 2
						dalfa_ipuro(j) = dalfa_puro(T, Tc);
					elseif ninput == 3
						if isa(self.EdE.m_inp, 'function_handle')             
							m = self.EdE.m_inp;
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
								m = self.EdE.m_fun(omega);
							end
							dalfa_ipuro = dalfa_puro(T, Tc, m);
						end
					end
				else
					dalfa_puro = self.EdE.dalfa_fun(T, comp(j));
					dalfa_ipuro(j) = dalfa_puro;
				end
            end
        end
        function [aij] = aij_fun(self, Mezcla, ai_puro)
            num_sust = Mezcla.num_sust;
            k = Mezcla.kmatrix;
            aij = zeros(num_sust, num_sust);
            for i = 1:num_sust
                for j = 1:num_sust
                    aij(i,j) = (1-k(i,j))*sqrt(ai_puro(i)*ai_puro(j));
               end
            end
        end
    end
end