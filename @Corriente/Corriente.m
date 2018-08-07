classdef Corriente<handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENGLISH EXPLANATION
%
% CORRIENTE Object handling different Comp, its phase or mixture,
%     its characteristics, its composition
%       Luis Jesús Díaz Manzo
%
% Current (Mixture, arg1, strarg1, arg2, strarg2, flow, strux, RMEdE, name, density)
% 		The properties you get are:
%           self.id: Name of the current (input arg. "name")
%
% 	Arguments:
% 		Input:
% 		Mixing: Mixing class.m containing the properties of the mixture as
% composition, x classes Substance.m defining the properties of pure substance
% 	Intensive properties of the current: T[K], P[kPa]
% 		arg1: Intensive property of the current (TEMPERATURE, PRESSURE OR VAPORIZED FRACTION). Specified which is strarg1.
% 		strarg1: String. It can take values of `T', `P', `x', to indicate the intensive variable provided was the temperature, pressure or vaporized fraction respectively.
% 		arg2: Intensive property of the current (TEMPERATURE, PRESSURE OR VAPORIZED FRACTION). Specified which is strarg1.
% 		strarg2: String. It can take values of `T', `P', `x', to indicate the
% 	intensive variable provided was the temperature, pressure or vaporized fraction respectively
%
% 	Extensive current properties: Flow rate[kg-mol/K]
% 		flow: Flow of the current. It can be molar (recommended), mass, or volumetric.
%
% 	Methods can also be used to define intensive and extensive variables.
%
% 	Once the degrees of freedom have been completed, the current is resolved to give the rest of the properties of the current (molar weight, phase, liquid and vapour phase compositions, the missing intensive variable, density, mass, molar and volumetric flows, etc.).
%
% 		One of the methods is to obtain the degrees of freedom of the current current object of current class (Current.m) which accepts multiple composite class objects as input and counts live the degrees of freedom needed to calculate the rest of the properties that compose it. At the moment it only estimates mixtures of liquid and vapour and not compromised and overheated, but it can be extended as soon as
%     has another, more advanced thermodynamic model, and is not limited by any combination of well-defined degrees of freedom
%     For example, by creating the set of composite class objects such as:
%       mixture = [ Substance('propane'), Substance('2-methylpropane'), ...
%           Substance('butane'), Substance('2-methylbutane'), ...
%           Substance ('pentane') ];
%
% 	A'Mix' of a given composition is created
%
%     	Mix = Mix(Mix,[0.05, 0.15, 0.25, 0.2, 0.35])
%
% 	A Peng-Robinson thermodynamic model is created with Van der Waals mixing rules
%
% 		Peng = PREdE();
% 		RM = RMVdW(Peng)
%
% 	to create a current'S001' with the mixture is quoted:
%       Current = Current(Mix, [], [], [], [], [], [], [], RM, 'S001')
%     By default, none of the degrees of freedom are defined, which specifies as gdlibertad = 8
%     If you made a mistake in assigning a variable, the same method used with the value "false" resets the current state
%       micorriente.conc(false);
%     If the molar compositions are not standardized, the program assumes that they have been introduced into molar flows and properly normalizes the composition and assigns the molar molF, reducing the degrees of freedom.
%     Now defining the remaining degrees of freedom, at any given time.
%     combination, the rest of the properties are calculated
%       mycurrent.temp( 79.938+273.15);
%       mycurrent.vaporized( 0.01);
%     The degrees of freedom were reduced to 0, so that the
%     liquid and vapor compositions, and pressure.
%     Once all variables have been calculated, the object may have
%     problems by subsequently resetting the values of the variables a
%     to one, for that there is a method that resets to the original state the
%     current:
%       mycurrent.reset()
%     It was also possible to determine the complete state of the current in its delineatio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXPLICACIÓN EN ESPAÑOL
%%%%%%%%%%%%%%%%%%%%%%
%
% Corriente Objeto que maneja diferentes Comp, su fase o mezcla,
%     sus características, su composición
%       Luis Jesús Díaz Manzo
%
% Corriente(Mezcla, arg1, strarg1, arg2, strarg2, flujo, strflujo, RMEdE, nombre, densidad)
% 		Las propiedades que obtiene son:
%           self.id: Nombre de la corriente (input arg. "nombre")
%
% 	Argumentos:
% 		Entrada:
% 		Mezcla: Clase Mezcla.m que contiene las propiedades de la mezcla como
% composición, x clases Sustancia.m que definen las propiedades de sustancia pura
% 	Propiedades intensivas de la corriente: T [K], P [kPa]
% 		arg1: Propiedad intensiva de la corriente (TEMPERATURA, PRESION O FRACCIÓN VAPORIZADA). Se especifica cual en strarg1.
% 		strarg1: String. Puede tomar valores 'T', 'P', 'x', para indicar la variable intensiva provista fue la temperatura, presión o fracción vaporizada respectivamente
% 		arg2: Propiedad intensiva de la corriente (TEMPERATURA, PRESION O FRACCIÓN VAPORIZADA). Se especifica cual en strarg1.
% 		strarg2: String. Puede tomar valores 'T', 'P', 'x', para indicar la
% 	variable intensiva provista fue la temperatura, presión o fracción vaporizada respectivamente
%
% 	Propiedades extensivas de la corriente: Flujo [kg-mol/K]
% 		flujo: Flujo de la corriente. Puede ser molar (recomendado), másico, o  volumétrico.
%
% 	También puede hacerse uso de los métodos para definir las variables intensivas y extensivas.
%
% 	Una vez se han completado los grados de libertad, la corriente se resuelve para dar el resto de las propiedades de la corriente (peso molar, fase, composiciones de las fases líquida y vapor, la variable intensiva faltante, la densidad, los flujos másicos, molares y volumétricos
%
% 		Uno de los métodos es la obtención de los grados de libertad de la corriente objeto de clase corriente (Corriente.m) que acepta como entrada múltiples objetos de clase compuesto y cuenta en vivo los grados de libertad que se necesitan para calcular el resto de las propiedades que la componen. Por el momento solo estima mezclas de líquido y vapor y no compromido y sobrecalentado, pero se puede extender en cuanto
%     tenga otro modelo termodinámico más avanzado, y no se ve limitado por  ninguna combinación de grados de libertad que se encuentren bien determinados
%     Por ejemplo,  creando el conjunto de objetos de clase compuesto como:
%       mezcla = [ Sustancia('propane'), Sustancia('2-methylpropane'), ...
%           Sustancia('butane'), Sustancia('2-methylbutane'), ...
%           Sustancia('pentane') ];
%
% 	Se crea una mezcla 'Mix' de cierta composición dada
%
%     	Mix = Mezcla(mezcla, [0.05, 0.15, 0.25, 0.2, 0.35])
%
% 	Se crea un modelo termodinámico de Peng-Robinson con reglas de mezclado de Van der Waals
%
% 		Peng = PREdE();
% 		RM = RMVdW(Peng)
%
% 	para crear una corriente 'S001' con la mezcla se cita:
%       micorriente = Corriente(Mix, [], [], [], [], [], [], RM, 'S001')
%     Por defecto no tiene definidas ninguno de los grados de libertad, que especifica como gdlibertad = 8
%     Si se equivocó en la asignación de una variable, el mismo método usado con el valor "false" resetea el estado de la corriente
%       micorriente.conc(false);
%     Si las composiciones molares no se encuentran normalizadas, el programa supone que se ha introducido en flujos molares y apropiadamente normaliza la composición y asigna el molF molar, reduciendo los grados de libertad
%     Ahora definiendo los grados de libertad restantes, en cualquier
%     combinación, se calculan el resto de las propiedades
%       micorriente.temp( 79.938+273.15 );
%       micorriente.vaporizado( 0.01 );
%     Se redujo los grados de libertad a 0, por lo que se calcularon las
%     composiciones de líquido y de vapor, y la presión.
%     Una vez que se calcularon todas las variables el objeto puede tener
%     problemas reseteando posteriormente los valores de las variables una
%     a una, para eso hay un método que resetea al estado original la
%     corriente:
%       micorriente.reset()
%     También se pudo determinar el estado completo de la corriente en su definición

    micorriente = Corriente(Mix, 79.938+273.15, 't', 0.01, 'x', 1, 'm', RM, 'S001')

    properties
        id = '';
        fase = false; %puede tomar valores mezcla, l�q o vap
        mezcla = Mezcla.empty(0,0);
        conc = false;  %Composicion global, z, de los Comp
        molF = false;  %molF en unidades molares
        unidflujo = 'kg-mol / h';
        mw = false;
        unidmasamolar = 'kg / kg-mol';
        masaF = false;
        unidflujomasa = 'kg / h';
        volF = false;
        unidflujovolumen = 'm^3 / h';
        densidad = false;
        uniddensidad = 'kg-mol / m^3';
        T = false;
        unidtemperatura = 'K';
        P = false;
        unidpresion = 'kPa';
        beta = false; %beta fraccion vaporizada de la corriente
		K = false;
        conc_liq = false;
        conc_vap = false;
        gdlibertad = 3;
        comp = false;
        num_sust = false;
        MEdE = false;
        q = false;
        H = false;
    end

    methods
        function self = Corriente(mezc, T_o_P_o_beta, tpb1, P_o_T_o_beta, ptb2, flujo, m_o_w_o_v , RMEdE, nombre, varargin)
            if nargin > 0
                if ~isempty(mezc) && isa(mezc, 'Mezcla')
                    self.mezcla = mezc;
                    self.comp = mezc.comp;
                    self.num_sust = mezc.num_sust;
                    if ~isempty(mezc.comp);
                        self.conc = mezc.conc;
                        self.gdlibertad = 3;
                        pesomolar = 0;
                        for i = 1: self.num_sust
                            pesomolar = pesomolar + self.conc(i).*self.comp(i).molmass;
                        end
                        self.mw = pesomolar;
                    else
                        self.gdlibertad = 3 + self.num_sust - 1;%Los grados de libertad
                %son N-1 por composiciones, 2 de (P,T,beta) y molF extensivo
                    end
                end
            end
            if nargin > 2
                if ~isempty(T_o_P_o_beta)
                    if strcmpi(tpb1, 'T')
                        self.T = T_o_P_o_beta;
                        self.gdlibertad = self.gdlibertad - 1;
                    elseif strcmpi(tpb1, 'P')
                        self.P = T_o_P_o_beta;
                        self.gdlibertad = self.gdlibertad - 1;
                    elseif strcmpi(tpb1, 'x')
						if T_o_P_o_beta < 0
							T_o_P_o_beta = 1e-11;
						elseif T_o_P_o_beta > 1
							T_o_P_o_beta = 1 - 1e-11;
						end
                        if T_o_P_o_beta == 1
							T_o_P_o_beta = 1-1e-11;
						elseif T_o_P_o_beta == 0
							T_o_P_o_beta = 1e-11;
                        end
						self.beta = T_o_P_o_beta;
                        self.gdlibertad = self.gdlibertad - 1;
                    else
                        error('Matlab:Corriente', ['La variable ''T_o_P_o_beta'': "' T_o_P_o_beta '" debe ser ''T'', ''P'', o ''x''. El valor "' tpb1 '" es incorrecto'])
                    end
                end
            end
            if nargin > 4
                if ~isempty(P_o_T_o_beta)
                    if strcmpi(ptb2, 'T')
                        self.T = P_o_T_o_beta;
                        self.gdlibertad = self.gdlibertad - 1;
                    elseif strcmpi(ptb2, 'P')
                        self.P = P_o_T_o_beta;
                        self.gdlibertad = self.gdlibertad - 1;
                    elseif strcmpi(ptb2, 'x')
                        if P_o_T_o_beta == 1
							P_o_T_o_beta = 1-1e-11;
						elseif P_o_T_o_beta == 0
							P_o_T_o_beta = 1e-11;
                        end
						self.beta = P_o_T_o_beta;
                        self.gdlibertad = self.gdlibertad - 1;
                    else
                        error('Matlab:Corriente',['La variable ''P_o_T_o_beta'': "' P_o_T_o_beta '" debe ser ''T'', ''P'', o ''x''. El valor "' ptb2 '" es incorrecto'])
                    end
                end
            end
            if nargin > 6
                if ~isempty(flujo)
                    if strcmpi(m_o_w_o_v, 'm')
                        self.molF = flujo;
                        self.gdlibertad = self.gdlibertad - 1;
                    elseif strcmpi(m_o_w_o_v, 'w')
                        self.masaF = flujo;
                        self.gdlibertad = self.gdlibertad - 1;
                    elseif strcmpi(m_o_w_o_v, 'v')
                        self.volF = flujo;
                        self.gdlibertad = self.gdlibertad - 1;
                    else
                        error('Matlab:Corriente',['La variable ''flujo'': "' flujo '" debe ser ''m'', ''w'', o ''v''. El valor "' m_o_w_o_v '" es incorrecto'])
                    end
                end
            end
            if nargin > 7
                if ~isempty(RMEdE) && isa(RMEdE, 'IdealEdE')
                    self.MEdE = RMEdE;
                end
            end
            if nargin > 8
              self.id = nombre;
            end
            if nargin > 9
                self.densidad = varargin{1};
            end
            self.determina_corriente();
        end
        function set_composic = x(self, z) %asignar concentraci�n molar
            %o composicion molar, si no est�n normalizados, los normaliza y
            %asigna adicionalmente un molF molar asociado.
            if ~islogical(z)
                if ~any(self.conc)
                    self.gdlibertad = self.gdlibertad - self.num_sust + 1;
                end
                tamano = size(z);
                if tamano(1) > tamano(2)
                    z = abs(z'); %vector fila de concentraciones positivas
                else
                    z = abs(z);
                end
                tamano = size(z);
                suma = sum(z);
                if tamano(2) == self.num_sust - 1 %Si se provee las primeras N-1 concentraciones independientes calcula la dependiente
                    if suma > 1 && suma < 100 %Si no se proveen en base unitaria se supone base 100
                        warning('Concentraciones no son fracciones molares en base unitaria. Se supondr� base 100');
                        self.conc = z./100;
                        self.conc(end + 1) = (100 - suma)./100;
                    elseif suma > 100
                        warning('Concentraciones no son fracciones molares en base unitaria. Se supone un(varios) componente(s) no presente(s)');
                        self.conc = z./suma;
                        self.conc(end + 1) = 0;
                    else
                        self.conc = z;
                        self.conc(end + 1) = 1 - suma;
                    end
                elseif tamano(2) == self.num_sust %Si se proveen las N concentraciones
                    if suma ~= 1
                        warning('Las concentraciones deben proveerse en fracciones molares en base unitaria. Se procede a su conversi�n');
                        self.conc = z./suma;
                    else
                        self.conc = z;
                    end
                else
                    error('Mala especificaci�n de las concentraciones molares. Sobran o faltan componentes');
                end
                pesomolar = 0;
                for i = 1: self.num_sust
                    pesomolar = pesomolar + self.conc(i).*self.comp(i).pesomolar;
                end
                self.mw = pesomolar;
            else
                self.conc = false;
                self.mw = false;
                self.gdlibertad = self.gdlibertad + self.num_sust - 1;
            end
            %si se completaron los grados de libertad con la asignaci�n,
            %entonces resuelve la corriente
            self.determina_corriente();
            set_composic = self;%salida es una copia de la Corriente actual
        end
        function set_temperature = temp(self, tmp) %asignar T a la
            %corriente
            if ~islogical(tmp)
                if ~self.T
                    self.gdlibertad = self.gdlibertad - 1;
                end
                self.T = tmp;
            else
                self.T = false;
                self.gdlibertad = self.gdlibertad + 1;
            end
            self.determina_corriente();
            set_temperature = self;
        end
        function set_pressure = pres(self, presion)
            if ~islogical(presion)
                if ~self.P
                    self.gdlibertad = self.gdlibertad - 1;
                end
                self.P = presion;
            else
                self.P = false;
                self.gdlibertad = self.gdlibertad + 1;
            end
            self.determina_corriente();
            set_pressure = self;
        end
        function set_fraccionvapor = vaporizado(self, x)
            if ~islogical(x)
                if ~self.beta
                    self.gdlibertad = self.gdlibertad - 1;
                end
				if x < 0
					x = 1e-11;
				elseif x > 1
					x = 1 - 1e-11;
				end
                if x == 1
                    x = 1-1e-11;
                elseif x == 0
                    x = 1e-11;
                end
                self.beta = x;
                self.fase = 'mezcla';
            else
                self.beta = false;
                self.gdlibertad = self.gdlibertad + 1;
            end
            self.determina_corriente();
            set_fraccionvapor = self;
        end
        function set_flowrate = flujomolar(self, F)
            if ~islogical(F)
                if ~self.molF
                    self.gdlibertad = self.gdlibertad - 1;
                end
                self.molF = F;
            else
                self.molF = false;
                self.gdlibertad = self.gdlibertad + 1;
            end
            self.determina_corriente();
            set_flowrate = self;
        end
        function set_flowrate = flujovolumetrico(self, F)
            if ~islogical(F)
                if ~self.molF
                    self.gdlibertad = self.gdlibertad - 1;
                end
                self.volF = F;
                if ~islogical(self.T)
                    self.molF = F*self.densidad;
                end
            else
                self.molF = false;
                self.gdlibertad = self.gdlibertad + 1;
            end
            self.determina_corriente();
            set_flowrate = self;
        end
        function set_flowrate = flujomasico(self, F)
            if ~islogical(F)
                if ~self.molF
                    self.gdlibertad = self.gdlibertad - 1;
                end
                self.masaF = F;
                self.molF = F/self.mw;
            else
                self.molF = false;
                self.gdlibertad = self.gdlibertad + 1;
            end
            self.determina_corriente();
            set_flowrate = self;
        end
        function determina_corriente(self)
            if self.gdlibertad == 0
                if any(self.conc) && self.P && (self.molF || isa(self.molF, 'double'))...
                        && self.T
                    %%%% ---Calculos de los puntos de Burbuja y de Rocio
                    [Tb, sup_y, K1, flag1] = ...
                        self.MEdE.BubbleT(self.P, self.mezcla);
                    [Tb, y1, K1, flag1] = ...
                        self.MEdE.BubbleT(self.P, self.mezcla, sup_y, [], [], Tb);
                    [Tr, sup_x, K2, flag2] = ...
                        self.MEdE.DewT(self.P, self.mezcla);
                    [Tr, x1, K2, flag2] = ...
                        self.MEdE.DewT(self.P, self.mezcla, sup_x, [],[], Tr);
					[Pb, sup_y, K3, flag3] = ...
						self.MEdE.BubbleP(self.T, self.mezcla);
                    if Pb == 0 || any(isnan(sup_y))
                        Pb = self.P/2;
                    else
                        [Pb, y2, K3, flag3] = ...
                            self.MEdE.BubbleP(self.T, self.mezcla, sup_y);
                    end
					[Pr, sup_x, K4, flag4] = ...
						self.MEdE.DewP(self.T, self.mezcla);
                    if Pr == 0 || any(isnan(sup_y))
                        Pr = self.P*2;
                    else
                        [Pr, x2, K4, flag4] = ...
                            self.MEdE.DewP(self.T, self.mezcla, sup_x);
                    end
                    %%%%% -------    %%%%%
                    flag5 = false; %Si se hizo flash cambia de valor
                    if self.T > Tb && self.T < Tr && self.P < Pb && self.P > Pr
						[self.conc_liq, self.conc_vap, ...
                            self.beta, self.K, flag5, ~] = self.MEdE.flash(...
							self.T, self.P, self.mezcla, x1, y1, 0.5);
						self.fase = 'mezcla';
                        self.q = 1 - self.beta;
                        self.H = self.entalpia();
                    elseif self.T < Tb
                        self.conc_liq = self.conc;
                        self.conc_vap = zeros(1, self.num_sust);
                        self.beta = 0;
						self.K = K1;
                        self.fase = 'liq';
                    elseif self.T > Tr
                        self.conc_vap = self.conc;
                        self.conc_liq = zeros(1, self.num_sust);
                        self.beta = 1;
                        self.fase = 'vap';
						self.K = K2;
                    elseif self.P < Pr
						self.conc_vap = self.conc;
						self.conc_liq = zeros(1, self.num_sust);
						self.beta = 1;
						self.fase = 'vap';
						self.K = K4;
                    elseif self.P > Pb
						self.conc_liq = self.conc;
						self.conc_vap = zeros(1, self.num_sust);
						self.beta = 0;
						self.fase = 'liq';
						self.K = K3;
                    end
                    if islogical(flag5)
                        Hdep_bub = self.MEdE.entalpia(Tb, self.P, self.mezcla, 'liq');
                        Hdep_dew = self.MEdE.entalpia(Tr, self.P, self.mezcla, 'vap');
                        Hvap = Hdep_dew - Hdep_bub;
                        if any(self.conc_liq)
%                                 deltaH = integral(@(t) cp(t), self.T, Tb);
                                Hliqsat = Hdep_bub;
                                Hliq = self.MEdE.entalpia(self.T, self.P, self.mezcla, 'liq');
                                deltaH = -Hliq + Hliqsat;
                                Hlideal = deltaH;
                            self.q = (Hvap + Hlideal)/Hvap;
                        elseif any(self.conc_vap)
%                                 deltaH = integral(@(t) cp(t), self.T, Tr);
                                Hvapsat = Hdep_dew;
                                Hvapor = self.MEdE.entalpia(self.T, self.P, self.mezcla, 'vap');
                                deltaH = Hvapsat - Hvapor;
                                Hgideal = deltaH;
                            self.q = (Hgideal)/Hvap;
                        end
%                         Hgi = zeros(1, self.num_sust);
%                         Href = zeros(1, self.num_sust);
%                         for i = 1:self.num_sust
%                             Href(i) = self.comp(i).href;
%                             try
%                                 cp = self.comp(i).cp_gi{1};
%                             catch ME
%                                 error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
%                             end
%                             deltaH = integral(@(t) cp(t), 273.15, self.T);
%                             Hgi(i) = deltaH*self.conc(i);
%                         end
%                         Href = self.MEdE.entalpia(273.15, self.P, self.mezcla, 'liq');
                        self.H = self.entalpia();
                    end
                elseif any(self.conc) && self.P && (self.molF || isa(self.molF, 'double') )...
                        && ((1 - 1.1e-11) < self.beta)
                    self.conc_vap = self.conc;
					[Tr, sup_x, ~, flag] = ...
                        self.MEdE.DewT(self.P, self.mezcla);
                    [self.T, self.conc_liq, self.K, flag] = ...
                        self.MEdE.DewT(self.P, self.mezcla, sup_x, [], [], Tr);
					self.fase = 'mezcla';
                    self.q = 0;
                    self.H = self.entalpia();
                elseif any(self.conc) && self.P && (self.molF || isa(self.molF, 'double'))...
                        && (1.1e-11 > self.beta)
                    self.conc_liq = self.conc;
                    [Tb, sup_y, ~, flag] = ...
                        self.MEdE.BubbleT(self.P, self.mezcla);
                    [self.T, self.conc_vap, self.K, flag] = ...
                        self.MEdE.BubbleT(self.P, self.mezcla, sup_y, [], [], Tb);
					self.fase = 'mezcla';
                    self.q = 1;
                    self.H = self.entalpia();
				elseif any(self.conc) && self.T && (self.molF || isa(self.molF, 'double') )...
                        && (1.1e-11 > self.beta)
                    self.conc_liq = self.conc;
                    self.fase = 'mezcla';
					[Pb, sup_y, ~, flag] = ...
						self.MEdE.BubbleP(self.T, self.mezcla);
                    [self.P, self.conc_vap, self.K, flag] = ...
						self.MEdE.BubbleP(self.T, self.mezcla, sup_y, [],[], Pb);
                    self.q = 1;
                    self.H = self.entalpia();
                elseif any(self.conc) && self.T &&(self.molF || isa(self.molF, 'double'))...
                        && ((1 - 1.1e-11) < self.beta)
                    self.conc_vap = self.conc;
                    self.fase = 'mezcla';
                    [Pr, sup_x, ~, flag] = ...
						self.MEdE.DewP(self.T, self.mezcla);
                    [self.P, self.conc_liq, self.K, flag] = ...
						self.MEdE.DewP(self.T, self.mezcla, sup_x, [], [], Pr);
                    self.q = 0;
                    self.H = self.entalpia();
                elseif any(self.conc) && self.P && (self.molF || isa(self.molF, 'double')) ...
                        && self.beta
                    max_iter = 100;
					tol = 1e-6;
                    [Tb1, sup_y, ~, flag1] = self.MEdE.BubbleT(self.P,...
						self.mezcla);
                    [Tb1, sup_y, K1, flag2] = self.MEdE.BubbleT(self.P,...
						self.mezcla, sup_y, [], [], Tb1);
                    [Tr1, sup_x, ~, flag3] = self.MEdE.DewT(self.P,...
						self.mezcla);
                    [Tr1, sup_x, K2, flag4] = self.MEdE.DewT(self.P,...
						self.mezcla, sup_x, [], [], Tr1);
                    Tsup = Tb1 + self.beta*(...
                        Tb1 - Tr1);
                    dT = 0.03.*Tsup;
                    iter = 0;
                    valI = 1;
                    while valI > tol && (iter < max_iter)
                        iter = iter + 1;
						y_ant = sup_y;
						x_ant = sup_x;
						[sup_x, sup_y, b, K1] = self.MEdE.flash(Tsup,...
                            self.P, self.mezcla, sup_x, sup_y, self.beta);
                        suma = self.beta - b;
						valI = 0;
						for i=1:self.num_sust
							valI = valI + abs(sup_x(i) - x_ant(i)) + abs(sup_y(i) - y_ant(i));
						end
                        if iter > 1
							if sign(suma)~= sign(suma_ant)
								dT = dT/2;
							end
                        end
						Tsup = Tsup + sign(suma)*dT;
                        suma_ant = suma;
                        self.q = 1 - self.beta;
                    end
					self.T = Tsup;
					self.conc_vap = sup_y;
					self.conc_liq = sup_x;
					self.K = K1;
					self.fase = 'mezcla';
                    self.H = self.entalpia();
                elseif any(self.conc) && self.P && (self.volF || isa(self.volF, 'double')) ...
                        && (self.beta < 1.1e-11 && self.beta > 0)
                    self.conc_liq = self.conc;
                    [Tb, sup_y, ~, flag] = ...
                        self.MEdE.BubbleT(self.P, self.mezcla);
                    [self.T, self.conc_vap, self.K, flag] = ...
                        self.MEdE.BubbleT(self.P, self.mezcla, sup_y, [], [], Tb);
					self.fase = 'mezcla';
                    volumen = self.MEdE.volumesp(self.T, self.P, self.mezcla, 'liq').*self.conc_liq;
                    volumen = sum(volumen);
                    density = (volumen)^-1;
                    self.densidad = density;
                    self.molF = self.volF.*self.densidad;
                    self.q = 1;
                    self.H = self.entalpia();
                elseif any(self.conc) && self.P && (self.volF || isa(self.volF, 'double')) ...
                        && (self.beta > 1-1.1e-11 && self.beta < 1)
                    self.conc_vap = self.conc;
                    [Tr, sup_x, ~, flag] = ...
                        self.MEdE.DewT(self.P, self.mezcla);
                    [self.T, self.conc_liq, self.K, flag] = ...
                        self.MEdE.DewT(self.P, self.mezcla, sup_x, [], [], Tr);
					self.fase = 'mezcla';
                    volumen = self.MEdE.volumesp(self.T, self.P, self.mezcla, 'vap').*self.conc_vap;
                    volumen = sum(volumen);
                    density = (volumen)^-1;
                    self.densidad = density;
                    self.molF = self.volF.*self.densidad;
                    self.q = 0;
                    self.H = self.entalpia();
                elseif any(self.conc) && self.P && (self.volF || isa(self.volF, 'double')) ...
                        && (self.beta)
                    max_iter = 100;
					tol = 1e-6;
                    [Tr1, sup_x, self.K, flag] = ...
                        self.MEdE.DewT(self.P, self.mezcla);
                    [Tr1, sup_x, self.K, flag] = ...
                        self.MEdE.DewT(self.P, self.mezcla, sup_x, [], [], Tr1);
                    [Tb1, sup_y, self.K, flag] = ...
                        self.MEdE.BubbleT(self.P, self.mezcla);
                    [Tb1, sup_y, self.K, flag] = ...
                        self.MEdE.BubbleT(self.P, self.mezcla, sup_y, [], [], Tb1);
					Tsup = Tb1 + self.beta*(...
                        Tb1 - Tr1);
                    dT = 0.03.*Tsup;
                    iter = 0;
                    valI = 1;
                    while valI > tol && (iter < max_iter)
                        iter = iter + 1;
						y_ant = sup_y;
						x_ant = sup_x;
						[sup_x, sup_y, b, K1] = self.MEdE.flash(Tsup,...
                            self.P, self.mezcla, sup_x, sup_y, self.beta);
                        suma = self.beta - b;
						valI = 0;
						for i=1:self.num_sust
							valI = valI + abs(sup_x(i) - x_ant(i)) + abs(sup_y(i) - y_ant(i));
						end
                        if iter > 1
							if sign(suma)~= sign(suma_ant)
								dT = dT/2;
							end
                        end
						Tsup = Tsup + sign(suma)*dT;
                        suma_ant = suma;
                    end
					self.T = Tsup;
					self.conc_vap = sup_y;
					self.conc_liq = sup_x;
					self.K = K1;
                    self.fase = 'mezcla';
                    volumen = self.MEdE.volumesp(self.T, self.P, self.mezcla, 'vap').*self.conc_vap.*(self.beta) + self.MEdE.volumesp(self.T, self.P, self.mezcla, 'liq').*self.conc_liq.*(1-self.beta);
                    volumen = sum(volumen);
                    density = (volumen)^-1;
                    self.densidad = density;
                    self.molF = self.volF.*self.densidad;
                    self.q = 1 - self.beta;
                    self.H = self.entalpia();
                elseif any(self.conc) && self.T && self.molF ...
                        && self.beta
                    max_iter = 100;
					tol = 1e-6;
                    [Pb1, sup_y, ~] = self.MEdE.BubbleP(self.T,...
						self.mezcla);
                    [Pb1, sup_y, ~] = self.MEdE.BubbleP(self.T,...
						self.mezcla, sup_y, [], [], Pb1);
                    [Pr1, sup_x, ~] = self.MEdE.DewP(self.T,...
						self.mezcla);
                    [Pr1, sup_x, ~] = self.MEdE.DewP(self.T,...
						self.mezcla, sup_x, [], [], Pr1);
                    Psup = Pb1 + self.beta*(Pb1 - Pr1);
                    dP = 0.03.*Psup;
                    iter = 0;
                    valI = 1;
                    while (valI > tol) && (iter < max_iter)
                        iter = iter + 1;
						y_ant = sup_y;
						x_ant = sup_x;
						[sup_x, sup_y, b, K1] = self.MEdE.flash(self.T,...
                            Psup, self.mezcla, sup_x, sup_y, self.beta);
                        suma = self.beta - b;
						valI = 0;
						for i=1:self.num_sust
							valI = valI + abs(sup_x(i) - x_ant(i)) + abs(sup_y(i) - y_ant(i));
						end
                        if iter > 1
							if sign(suma)~= sign(suma_ant)
								dP = dP/2;
							end
                        end
						Psup = Psup - sign(suma)*dP;
                        suma_ant = suma;
                        self.q = 1 - self.beta;
                    end
					self.P = Psup;
					self.conc_vap = sup_y;
					self.conc_liq = sup_x;
					self.K = K1;
					self.fase = 'mezcla';
                    self.q = 1 - self.beta;
                    self.H = self.entalpia();
                end
                self.ordenaK();
                pesomolar = 0;
                for i = 1: self.num_sust
                    pesomolar = pesomolar + self.conc(i).*self.comp(i).molmass;
                end
                self.mw = pesomolar;
                self.masaF = self.molF*self.mw;
                if islogical(self.densidad) && isa(self.MEdE, 'RMVdW')
                    volumen = self.MEdE.volumesp(self.T, self.P, self.mezcla, 'liq').*self.conc_liq.*(1-self.beta) + self.MEdE.volumesp(self.T, self.P, self.mezcla, 'vap').*self.conc_vap.*(self.beta);
                    volumen = sum(volumen);
                    density = (volumen)^-1;
                    if islogical(self.densidad)
                        self.densidad = density;
                    end
                end
                if islogical(self.volF) && isa(self.MEdE, 'RMVdW')
                    self.volF = self.molF./self.densidad;
                end
%             density = 0;
%             for i = 1: self.num_sust
%                 density = density + self.conc(i).*self.Comp(i).densidad(...
%                     % self.T);
            % end
            % self.densidad = density;
            % self.volF = self.molF ./ self.densidad;
            % self.masaF = self.molF .* self.mw;
            else
                if any(self.conc_liq) || any(self.conc_vap)
                   self.gdlibertad = self.gdlibertad - 1;
                   self.conc_liq = false;
                   self.conc_vap = false;
                   self.fase = '';
                end
            end
        end
        function ordenaK(self)
            if self.gdlibertad == 0
                [self.K, indice] = sort(self.K, 'descend');
                if self.K ~= 0
                    Compo = Sustancia.empty(0, length(indice));
                    concent = zeros(1, length(indice));
                    concentliq = zeros(1, length(indice));
                    concentvap = zeros(1, length(indice));
                    for i = 1:length(indice)
                        Compo(i) = self.comp(indice(i));
                        concent(i) = self.conc(indice(i));
                        concentliq(i) = self.conc_liq(indice(i));
                        concentvap(i) = self.conc_vap(indice(i));
                    end
                    kij = self.mezcla.kij;
                    self.mezcla.comp = Compo;
                    self.mezcla.conc = concent;
                    self.mezcla.kij = kij;
                    self.mezcla.matrixkij();
                    self.comp = Compo;
                    self.conc = concent;
                    self.conc_liq = concentliq;
                    self.conc_vap = concentvap;
                end
            end
        end
        function H = entalpia(self)
            if self.gdlibertad == 0
                if strcmp(self.fase,'mezcla')
                    Hgi = zeros(1, self.num_sust);
                    Href = zeros(1, self.num_sust);
                    for i = 1:self.num_sust
                        Href(i) = self.comp(i).href;
                        try
                            cp = self.comp(i).cp_gi{1};
                        catch ME
                            error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                        end
                        deltaH = integral(@(t) cp(t), 273.15, self.T);
                        Hgi(i) = deltaH*self.conc(i) + Href(i)*self.conc(i);
                    end
                    H = sum(Hgi);
%                     Liquido = Mezcla(self.mezcla.comp, self.conc_liq, self.mezcla.kij);
%                     Gas = Mezcla(self.mezcla.comp, self.conc_vap, self.mezcla.kij);
                    Hdep_bub = self.MEdE.entalpia(self.T, self.P, self.mezcla, 'liq');
                    Hdep_dew = self.MEdE.entalpia(self.T, self.P, self.mezcla, 'vap');

                    Hdep_ref = self.MEdE.entalpia(273.15, 101.325, self.mezcla, 'liq');
                    Hvap = -(Hdep_dew - Hdep_bub);
                    H = H + self.beta*(Hvap) - Hdep_bub + Hdep_ref;
                elseif strcmp(self.fase,'liq') || strcmp(self.fase,'vap')
                    Hdep = self.MEdE.entalpia(self.T, self.P, self.mezcla, self.fase);

                    Hdep_ref = self.MEdE.entalpia(273.15, 101.325, self.mezcla, 'liq');
                    Hgi = zeros(1, self.num_sust);
                    Href = zeros(1, self.num_sust);
                    for i = 1:self.num_sust
                        Href(i) = self.comp(i).href;
                        try
                            cp = self.comp(i).cp_gi{1};
                        catch ME
                            error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                        end
                        deltaH = integral(@(t) cp(t), 273.15, self.T);
                        Hgi(i) = deltaH*self.conc(i) + Href(i)*self.conc(i);
                    end
                    H = sum(Hgi);
                    H = H - Hdep + Hdep_ref;
                end
            else
                H = false;
                warning('La corriente a�n no puede ser determinada porque no se suministraron suficientes datos');
            end
            self.H = H;
        end
        function reset(self)
            self.id = '';
            self.fase = false; %puede tomar valores mezcla, l�q o vap
            self.mezcla = false;
            self.conc = false;  %Composicion global, z, de los Comp
            self.molF = false;  %molF en unidades molares
            self.unidflujo = 'kg-mol / h';
            self.mw = false;
            self.unidmasamolar = 'kg / kg-mol';
            self.masaF = false;
            self.unidflujomasa = 'kg / h';
            self.volF = false;
            self.unidflujovolumen = 'm^3 / h';
            self.densidad = false;
            self.uniddensidad = 'kg-mol / m^3';
            self.T = false;
            self.unidtemperatura = 'K';
            self.P = false;
            self.unidpresion = 'kPa';
            self.beta = false; %beta fraccion vaporizada de la corriente
            self.K = false;
            self.conc_liq = false;
            self.conc_vap = false;
            self.gdlibertad = 3;
            self.comp = false;
            self.num_sust = false;
            self.MEdE = false;
            self.q = false;
            self.H = false;
        end
    end

end
