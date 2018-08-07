classdef Shortcut < handle
    %Shorcut Objeto que hereda de columna y que maneja los métodos cortos
    %de Underwood, Fenske, Gilliland y Kirkbride para el dimensionamiento 
    %de columnas de destilación.
    %   Luis Jesús Díaz Manzo
    %   
    
    properties
        id = '';
        nreal = false;
        feedtray = false;
        num_sust = false;
        flujo = 1;
        units_flujo = 'kg-mol / h';
        Dflow = false;
        units_Dflow = 'kg-mol / h';
        Dconc = false;
        Bflow = false;
        units_Bflow = 'kg-mol / h';
        Bconc = false;
        reco_lk = false;
        reco_hk = false;
        reco_top = false;
        reco_bot = false;
        cfeed = false;
        cdest = false;
        cbot = false;
        pres_drop_top = 0;
        units_press_drop_top = 'kPa';
        pres_drop_bot = 0;
        units_press_drop_bot = 'kPa';
        lk = false;
        hk = false;
        nmin = false;
        rmin = false;
        reflujo = false;
        nteo = false;
        kirkbride = false;
        gdlibertad = 3;
        volat = false;
        volat_top = false;
        volat_bot = false;
        concFeed = false;
        concliqtope = false;
        concvapfondo = false;
        beta_feed = false;
        q = false;
        MEdE = false;
        WesterbergNmin = false;
        WesterbergRmin = false;
        SmokerR = false;
        SmokerS = false;
        Smoker = false;
        reco_und = false;
    end
    
    methods
        function self = Shortcut(alim, lightk, recup_top, heavyk, recup_bot, dest, bot, nombre)
            self.gdlibertad = 3;
            if nargin > 0;
                if ~isempty(alim) && isa(alim, 'Corriente')
                    self.cfeed = alim;
                    self.concFeed = alim.conc;
                    self.beta_feed = 1 - alim.q;
                    self.q = alim.q;
                    self.gdlibertad = self.gdlibertad - 1;
                    self.flujo = alim.molF;
                    self.MEdE = alim.MEdE;
                    self.num_sust = alim.num_sust;
                end
            end
            if nargin > 2   
                if ~isempty(recup_top) && ~isempty(lightk)
                    self.set_recover_top(lightk, recup_top);
                end
            end
            if nargin >  4
                if ~isempty(recup_bot) && ~isempty(heavyk)
                    self.set_recover_bottom(heavyk, recup_bot);
                end
			end
			if nargin > 5
				if ~isempty(dest) && isa(dest, 'Corriente')
					self.cdest = dest;
				end
			end
			if nargin > 6
				if ~isempty(bot) && isa(bot, 'Corriente')
					self.cbot = bot;
                end
			end
            if nargin > 7
                if ~isempty(nombre) && isa(nombre, 'char')
                    self.id = nombre;
                end
            end   
            %recuperacion tope clave liviano, recuperación fondo pesado, 
            %y definición de una corriente de alimentación o un conjunto de
            %volatilidades de tope y/o fondo o constantes para la torre.
            %Ademas de las composiciones arbitrarias para los componentes
            self = self.checkdistribute();
        end
        function salida = set_pressure_drop_top(self, pressure_drop, units)
            %Establece la caída de presión total de la región de platos de 
            %rectificación
            if ~islogical(pressure_drop)
                if nargin == 2
                    units = self.units_press_drop_top;
                else
                    pressure_drop = units_pressure(pressure_drop, units, self.units_press_drop_top);
                end
                self.pres_drop_top = pressure_drop; 
                salida = self;
            end
        end
        function salida = set_pressure_drop_bot(self, pressure_drop, units)
            %Establece la caída de presión total de la región de platos de 
            %despojo
            
            if ~islogical(pressure_drop)
                if nargin == 2
                    units = self.units_press_drop_bot;
                else
                    pressure_drop = units_pressure(pressure_drop, units, self.units_press_drop_bot);
                end
                self.pres_drop_bot = pressure_drop;
                salida = self;
            end
        end
        function salida = set_corriente_feed(self, alim)
            i= length(self.cfeed);            
            %Pueden haber varias alimentaciones, cfeed almacena
            %en un array todas las alimentaciones dadas.
            self.cfeed(i+1) = currentfeed;
            if currentfeed.gdlibertad == 0 
                self.cfeed = alim;
                self.concFeed = alim.conc;
                self.volat = alim.K;
                self.beta_feed = 1 - alim.q;
                self.q = alim.q;
                self.gdlibertad = self.gdlibertad - 1;
                self.flujo = alim.molF;
                self.MEdE = alim.MEdE;
            end
            salida = self;
        end
        function salida = set_corriente_distil(self, currentdist)
            i= length(self.cdest);            
            %Pudiera tener hasta 2 corrientes de destilado, una vaporizada
            %y otra a condensado líquido.
            self.cdest(i+1) = currentdist;
            salida = self;
        end
        function salida = set_corriente_residuo(self, current_bottom)
            self.cbot = current_bottom; 
            salida = self;
        end   
        function salida = set_recover_top(self, clave_liviano, recup_lk)
            if ischar(clave_liviano) || length(clave_liviano)==1
                if ~islogical(recup_lk) || ~islogical(clave_liviano)
                    if islogical(self.lk)
                        self.gdlibertad = self.gdlibertad - 1;
                    end
                    %Si usuario introduce valor numérico
                    if ischar(clave_liviano)
                        for i = 1:self.num_sust
                            if strcmpi(clave_liviano, self.cfeed.comp(...
                                    i).id)
                                self.lk = i; %se identifica ubicación del clave
                            end
                        end
                    else
                        self.lk = clave_liviano;
                    end
                    self.reco_lk = recup_lk; %Se identifica la recuperación del clave
                else
                    %Si el usuario introduce valor lógico y no ha sido asignado
                    %aún un valor para recuperación del liviano, es decir, el
                    %valor default que es "false" es interpretado como 0
                    self.reco_lk = false;
                    self.reco_top = false;
                    self.lk = false;
                    self.gdlibertad = self.gdlibertad + 1;
                end
            else 
                if ~islogical(recup_lk) || ~islogical(clave_liviano)
                    if islogical(self.lk)
                        self.gdlibertad = self.gdlibertad - 1;
                    end
                    %Si usuario introduce valor numérico
                    self.lk = zeros(1, length(clave_liviano));
                    for j=1:length(clave_liviano)
                        if ischar(clave_liviano{j})      
                            for i = 1:self.num_sust
                                if strcmpi(clave_liviano{j}, self.cfeed.comp(...
                                        i).id)
                                    self.lk(j) = i; %se identifica ubicación del clave
                                end
                            end
                        else
                            self.lk = clave_liviano;
                            break
                        end
                    end                          
                    self.reco_lk = recup_lk;
                else
                    %Si el usuario introduce valor lógico y no ha sido asignado
                    %aún un valor para recuperación del liviano, es decir, el
                    %valor default que es "false" es interpretado como 0
                    self.reco_lk = false;
                    self.reco_top = false;
                    self.lk = false;
                    self.gdlibertad = self.gdlibertad + 1;
                end
            end
            salida = self.checkdistribute();
        end
        function salida = set_recover_bottom( self, clave_pesado, recup_hk)    
            if ischar(clave_pesado) || length(clave_pesado)==1
                if ~islogical(recup_hk) || ~islogical(clave_pesado)
                    if islogical(self.hk)
                        
                        self.gdlibertad = self.gdlibertad - 1;
                    end
                    if ischar(clave_pesado)
                        for i = 1:self.num_sust
                            if strcmp(clave_pesado, self.cfeed.comp(...
                                    i).id)
                                self.hk = i;
                            end    
                        end
                    else
                        self.hk = clave_pesado;
                    end
                    self.reco_hk = recup_hk;
                elseif ~islogical(self.reco_hk)
                    %Si el usuario introduce valor lógico y no ha sido asignado
                    %aún un valor para recuperación del pesado, es decir, el
                    %valor default que es "false" es interpretado como 0
                    self.reco_hk = false;
                    self.hk = false;
                    self.gdlibertad = self.gdlibertad + 1;
                end
            else 
                if ~islogical(recup_hk) || ~islogical(clave_pesado)
                    if islogical(self.hk)
                        self.gdlibertad = self.gdlibertad - 1;
                    end
                    %Si usuario introduce valor numérico                    
                    self.hk = zeros(1, length(clave_pesado));
                    for j=1:length(clave_pesado)
                        if ischar(clave_pesado{j})
                            for i = 1:self.num_sust
                                if strcmpi(clave_pesado{j}, self.cfeed.comp(...
                                        i).id)
                                    self.hk(j) = i; %se identifica ubicación del clave
                                end
                            end
                        else
                            self.hk = clave_pesado;
                            break
                        end
                    end
                    self.reco_hk = recup_hk;
                else
                    %Si el usuario introduce valor lógico y no ha sido asignado
                    %aún un valor para recuperación del liviano, es decir, el
                    %valor default que es "false" es interpretado como 0
                    self.reco_hk = false;
                    self.reco_bot = false;
                    self.hk = false;
                    self.gdlibertad = self.gdlibertad + 1;
                end
            end
            salida = self.checkdistribute();
        end
        function salida = checkdistribute(self)
            if ~isempty(self.cfeed) && ~islogical(self.cfeed) && isa(self.cfeed, 'Corriente' )
                if self.gdlibertad == 0 && self.cfeed.gdlibertad == 0
                    self.concFeed = self.concFeed;
                    if length(self.hk) == 1 && length(self.lk) == 1
                        exit = self.distribute(...
                            self.reco_lk, self.reco_hk, max(self.lk), ...
                            min(self.hk));
                        self.nmin= exit{1};
                        self.reco_top = exit{2};
                        self.reco_bot = exit{3};            
                        self.Dflow = sum(self.concFeed.* ... 
                            self.reco_top.*self.flujo);
                        self.Bflow = sum(self.concFeed.* ...
                            self.reco_bot.*self.flujo);
                        self.Dconc = self.flujo.*self.reco_top.*self.concFeed...
                            ./ self.Dflow;                
                        self.Bconc = self.flujo.*self.reco_bot.*self.concFeed...
                            ./ self.Bflow;
                        self.rmin = self.get_underwood(self.reco_top, self.reco_bot);
                    else
                        self.tantearecuperaciones()
                    end
                    
                end
                
            elseif ~islogical(self.concFeed) && ~islogical(self.reco_lk) &&...
                ~islogical(self.reco_hk) && self.gdlibertad == 0 ...
                && ~islogical(self.volat_top) && ~islogical(...
                self.volat_bot)
                exit = distribute(self, self.reco_lk, self.reco_hk,...
                    self.lk, self.hk, self.volat_top,...
                    self.volat_bot);
                self.nmin= exit{1};
                self.reco_top = exit{2};
                self.reco_bot = exit{3};            
                self.Dflow = sum(self.concFeed.* ... 
                    self.reco_top.*self.flujo);
                self.Bflow = sum(self.concFeed.* ...
                    self.reco_bot.*self.flujo);
                self.Dconc = self.flujo.*self.reco_top.*self.concFeed...
                    ./ self.Dflow;                
                self.Bconc = self.flujo.*self.reco_bot.*self.concFeed...
                    ./ self.Bflow;
            end
            salida = self;
        end
        function salida = distribute(self, recuplk,...
                recuphk, light, heavy, volatil_top, volatil_bottom)
            if nargin ==5
                coefKhk = self.cfeed.K(heavy);
                %coefKlk y coefKhk verifican que haya suministrado el liviano
                %y el pesado en su orden adecuado de volat
                coefKlk = self.cfeed.K(light);
                if coefKhk > coefKlk
                    %Si se eligieron mal los livianos y pesados, se corrige
                    nuevo_lk = heavy;
                    heavy = light;
                    light = nuevo_lk;
                    recup_hk = recuplk;
                    recup_lk = recuphk;
                else                
                    recup_hk = recuphk;
                    recup_lk = recuplk;
                end
                recuplk = recup_lk;
                recuphk = recup_hk;
                if ~isa(self.volat, 'double')
                    if isa(self.volat_top, 'double') && isa(self.volat_bot, 'double')
                        self.volat = (self.volat_top.*self.volat_bot).^0.5;
                    elseif isa(self.volat_top, 'double') && ~isa(self.volat_bot, 'double')
                        self.volat = self.volat_top;
                    elseif isa(self.volat_bot, 'double') && ~isa(self.volat_top, 'double')
                        self.volat = self.volat_bot;
                    else
                        self.volat = self.cfeed.conc_vap ./ self.cfeed.conc_liq;
                    end
                end
                    alfa = self.volat ./ self.volat (heavy);
            elseif nargin == 6
                self.volat = volatil_top;
                self.volat_top = volatil_top;
                self.volat_bot = volatil_top;
                alfa = volatil_top ./ volatil_top(heavy);
            elseif nargin < 5 
                try
                    alfa = self.volat ./ self.volat(self.hk);
                catch ME
                    alfa = self.cfeed.K ./ self.cfeed.K(self.hk);
                end
                recuplk = self.reco_lk;
                recuphk = self.reco_hk;
            else 
                self.volat = (volatil_top.* volatil_bottom).^0.5;
                self.volat_top = volatil_top;
                self.volat_bot = volatil_bottom;
                alfa_top = volatil_top ./ volatil_top(heavy);
                alfa_bottom = volatil_bottom ./ volatil_bottom(heavy);
                alfa = (alfa_top.*alfa_bottom).^0.5;
            end
            %Ahora se va a calcular un estimado inicial de la distribución
            recup_top = zeros(1, self.num_sust);
            recup_bottom = zeros(1, self.num_sust);
            for i = 1 : self.num_sust
                if alfa(i) > alfa(light)
                    %Distribución de más ligeros que LK en base a
                    %volat relativa
                    recup_top(i) = (recuplk+(1 - recuplk)*(1 -...
                        1/(alfa(i))^(7)));
                    recup_bottom(i) = (1 - recup_top(i));
                elseif i == light
                    %Distribución del clave liviano
                    recup_bottom(i) = 1 - recuplk;
                    recup_top(i) = 1 - recup_bottom(i);
                elseif i == heavy 
                    %Distribución del clave pesado
                    recup_top(i)= 1 - recuphk;
                    recup_bottom(i) = 1 - recup_top(i);
                elseif alfa(i) < 1
                    %Distribución de más pesador que HK en base a
                    %volat relativa
                    recup_top(i) = (1 - recuphk)*(alfa(i)^(7));
                    recup_bottom(i) = 1 - recup_top(i);
                else
                    %Si hay compuestos de volat intermedia, se estima
                    %su distribución como 0.5 y se recalculará con Fenske y
                    %Underwood
                    recup_top(i) = 0.5;
                    recup_bottom(i) = 0.5;
                end    
            end
            flag = 1;
            iter = 0;
            while flag ~= 0
                iter = iter + 1;
                %Obtuve los estimados de distribución recupdist y
                %recupresidue los cuales corrigen los estimados utilizados
                if ~isempty(self.cfeed) && ~islogical(self.cfeed) && isa(self.cfeed, 'Corriente' )
                    if self.cfeed.gdlibertad == 0
                        self.Dflow = sum(self.flujo.*...
                            self.concFeed.*recup_top);
                        self.Bflow = sum(self.flujo.*...
                            self.concFeed.*recup_bottom);
                        xD = (self.flujo.*self.concFeed...
                            .*recup_top)./(self.Dflow);
                        xB = (self.flujo.*self.concFeed...
                            .*recup_bottom)./(self.Bflow);
                        [~, xtope, K1, flag1] = self.MEdE.DewT(self.cfeed.P - self.pres_drop_top, self.cfeed.mezcla, xD);
                        [~, yfondo, K2, flag2] = self.MEdE.BubbleT(self.cfeed.P + self.pres_drop_top, self.cfeed.mezcla, xB); 
                        KD = K1;
                        KB = K2;            
                        self.volat_top = KD ./ KD(heavy);
                        self.volat_bot = KB ./ KB(heavy);
                        self.Dconc = xtope;
                        self.Bconc = yfondo;
                        self.concliqtope = xtope;
                        self.concvapfondo = yfondo;
                    end
                end
                [mintrays, recupdist, recupresidue] = ...
                    self.getfenske(recup_top, recup_bottom, light, heavy,...
                    self.volat_top, self.volat_bot);
                for i=1:self.num_sust
                    if abs(recupdist(i) - recup_top(i)) > 5e-6 || (recup_top(i) > recupdist(i))
                        recup_top = recupdist;
                        recup_bottom = 1 - recup_top;
                        break
                    end
                end                
                if flag == 1
                    break
                end    
                flag = 1;
            end
            salida = {mintrays, recup_top, recup_bottom};
        end
        function [num_min_trays, recupdist, recupresidue] = getfenske(...
                self, recup_dist, recup_resid, light,...
                heavy, volatil_top, volatil_bottom)
            [num_min_trays, recupdist, recupresidue] = fenske(...
                self.concFeed, self.flujo, recup_dist,...
                recup_resid, light, heavy, volatil_top, volatil_bottom);            
        end
        function [minreflux, Dflow] = get_underwood(self, sup_recup_d, sup_recup_b,...
                volatil_top, volatil_bottom, light, heavy, beta)
            if nargin < 4
                if any(sup_recup_d == 0) || any(sup_recup_b == 1)
                    for i =1 : self.num_sust
                        if sup_recup_d(i) == 0
                            sup_recup_d(i) = 1e-15;
                        elseif sup_recup_d(i) == 1
                            sup_recup_d(i) = 1-1e-15;
                        end
                    end
                end
                if any(sup_recup_b == 0) || any(sup_recup_b == 1)
                    for i =1 : self.num_sust
                        if sup_recup_b(i) == 0
                            sup_recup_b(i) = 1e-15;
                        elseif sup_recup_d(i) == 1
                            sup_recup_d(i) = 1-1e-15;
                        end
                    end
                end
                
                if ~islogical(self.volat_top) && islogical(self.volat_bot)
                    alfa = self.volat_top ./ self.volat_top(self.hk);
                elseif ~islogical(self.volat_bot) && islogical(self.volat_top)
                    alfa = self.volat_bot ./ self.volat_bot(self.hk);
                elseif ~islogical(self.volat_bot) && ~islogical(self.volat_top)
                    alfa = (self.volat_bot.*self.volat_top).^0.5 ./ (self.volat_bot(self.hk).*self.volat_top(self.hk)).^0.5;
                else
                    D = sum(self.flujo.*self.concFeed.*sup_recup_d);
                    xD = (self.flujo.*self.concFeed.*sup_recup_d)./(D);
                    self.volat_top = xD./self.concliqtope;
                    alfa = self.volat_top ./ self.volat_top(self.hk);
                end
                beta = self.beta_feed;
                light = self.lk;
                heavy = self.hk;
            elseif nargin < 5                
                light = self.lk;
                heavy = self.hk;
                alfa = volatil_top ./ volatil_top(heavy);
                beta = self.beta_feed;
            elseif nargin < 6                
                light = self.lk;
                heavy = self.hk;
                beta = self.beta_feed;
                self.volat = (volatil_top);                
                self.volat_top = volatil_top;                
                self.volat_bot = volatil_bottom;                
                alfa = ((volatil_top .* volatil_bottom) ./ (volatil_top(self.hk)*volatil_bottom(self.hk))).^0.5;
            else
                self.volat = (volatil_top.* volatil_bottom).^0.5;
                self.volat_top = volatil_top;
                self.volat_bot = volatil_bottom;
                alfa = self.volat./self.volat(self.hk);
                self.beta_feed = beta;
                self.q = 1 - beta;
            end
                [minreflux, Dflow, incog] = underwood(self.concFeed, self.flujo,...
                    sup_recup_d, sup_recup_b, alfa, light, heavy, self.beta_feed);
            self.reco_und = incog(1:end-2);
        end
        function salida = get_gilliland(self, reflux, eficiencia)
            self.rmin = self.get_underwood(self.reco_top,...
                self.reco_bot,self.beta_feed, self.volat_top, ...
                self.volat_bot);
            if reflux > self.rmin 
                if nargin == 2
                    eficiencia = 0.5;
                end
                self.reflujo = reflux;
                abscisa = (reflux - self.rmin)/(reflux + 1);
                ordenada = 1 - exp(((1+54.4*abscisa)/(11+117.2*abscisa))*(...
                    (abscisa - 1)/(abscisa).^0.5));
                self.nteo = -(ordenada + self.nmin)/(ordenada - 1);
                self.nreal = ceil(self.nteo / eficiencia); 
                self.feedtray = round(2*self.get_kirkbride());
                salida = self;
            end
        end
        function salida = get_kirkbride(self)
            relacion_nr_ns = ((self.concFeed(self.hk)...
                /self.concFeed(self.lk))*(...
                self.Bconc(self.lk)/self.concFeed(...
                self.hk))^2*((self.Dflow)/(self.Bflow...
                )))^0.206;
            self.kirkbride = ((relacion_nr_ns)*(self.nteo))/((relacion_nr_ns)+1);
            salida = self.kirkbride;
        end
        function salida = set_composic_feed(self, composic)
            if islogical(self.cfeed) || islogical(self.cfeed) || ~isa(self.cfeed,'Corriente')
                self.concFeed = composic;
                self.num_sust = length(composic);
            else
                self.cfeed.conc = composic;
                self.concFeed = composic;
            end
            salida = self;
        end
        function salida = set_flow(self, flow)
            self.flujo = flow;
            salida = self;
        end
        function salida = set_constant_volatility(self, volatil)
            if islogical(self.volat)                
                self.gdlibertad = self.gdlibertad - 1;
            end
            self.volat = volatil;
            self.volat_top = volatil;
            self.volat_bot = volatil;
            if ~islogical(self.concFeed) && ~islogical(self.reco_lk) &&...
                    ~islogical(self.reco_hk)
                exit = distribute(self, self.reco_lk, self.reco_hk,...
                    self.lk, self.hk, self.volat_top, ...
                    self.volat_bot);
                self.nmin= exit{1};
                self.reco_top = exit{2};
                self.reco_bot = exit{3};      
                self.Dconc = self.concFeed.*... 
                    self.reco_top;
                self.Dconc = self.Dconc ./ sum(self.Dconc);
                self.Bconc = self.concFeed.*... 
                    self.reco_bot;                
                self.Bconc = self.Bconc ./ sum(self.Bconc);
                self.Dflow = sum(self.concFeed.*... 
                    self.reco_top.*self.flujo);
                self.Bflow = sum(self.concFeed.* ...
                    self.reco_bot.*self.flujo);
            end
            salida = self;
        end
        function salida = set_volatilidad_tope(self, volatil_top)
            if ~islogical(self.hk)
                self.volat_top = volatil_top./volatil_top(self.hk);
            else
                self.volat_top = volatil_top;
            end
            if length(self.volat_bot) == length(volatil_top)
                if islogical(self.volat)                
                    self.gdlibertad = self.gdlibertad - 1;
                end
                self.volat = (volatil_top.*self.volat_bot).^0.5;
                if self.gdlibertad == 0
                    if ~islogical(self.concFeed) && ~islogical(self.reco_lk) &&...
                        ~islogical(self.reco_hk)
                        if isa(self.volat_bot, 'double')
                            exit = distribute(self, self.reco_lk, self.reco_hk,...
                                self.lk, self.hk, self.volat_top, ...
                                self.volat_bot);
                        else 
                            exit = distribute(self, self.reco_lk, self.reco_hk,...
                                self.lk, self.hk, self.volat_top);
                        end
                        self.nmin= exit{1};
                        self.reco_top = exit{2};
                        self.reco_bot = exit{3};                    
                        self.Dconc = self.concFeed.*... 
                            self.reco_top;                    
                        self.Dconc = self.Dconc ./ sum(self.Dconc);
                        self.Bconc = self.concFeed.*... 
                            self.reco_bot;                    
                        self.Bconc = self.Bconc ./ sum(self.Bconc);
                        self.Dflow = sum(self.concFeed.*... 
                            self.reco_top.*self.flujo);
                        self.Bflow = sum(self.concFeed.* ...
                            self.reco_bot.*self.flujo);
                    end
                end
            else 
                if self.gdlibertad == 0
                    if ~islogical(self.concFeed) && ~islogical(self.reco_lk) &&...
                        ~islogical(self.reco_hk)
                        exit = distribute(self, self.reco_lk, self.reco_hk,...
                            self.lk, self.hk, self.volat_top);
                        self.nmin= exit{1};
                        self.reco_top = exit{2};
                        self.reco_bot = exit{3};                    
                        self.Dconc = self.concFeed.*... 
                            self.reco_top;                    
                        self.Dconc = self.Dconc ./ sum(self.Dconc);
                        self.Bconc = self.concFeed.*... 
                            self.reco_bot;                    
                        self.Bconc = self.Bconc ./ sum(self.Bconc);
                        self.Dflow = sum(self.concFeed.*... 
                            self.reco_top.*self.flujo);
                        self.Bflow = sum(self.concFeed.* ...
                            self.reco_bot.*self.flujo);
                    end
                end
            end
            salida = self;
        end
        function salida = set_volatilidad_fondo(self, volatil_bot)
            if ~islogical(self.hk)
                self.volat_bot = volatil_bot./volatil_bot(self.hk);
            else
                self.volat_bot = volatil_bot;
            end
            if length(self.volat_top) == length(volatil_bot)
                if islogical(self.volat)                
                    self.gdlibertad = self.gdlibertad - 1;
                end
                self.volat = (volatil_bot.*self.volat_top).^0.5;
                if self.gdlibertad == 0;
                    if ~islogical(self.concFeed) && ~islogical(self.reco_lk) &&...
                        ~islogical(self.reco_hk)
                        if isa(self.volat_top, 'double')
                            exit = distribute(self, self.reco_lk, self.reco_hk,...
                                self.lk, self.hk, self.volat_top, ...
                                self.volat_bot);
                        else 
                            exit = distribute(self, self.reco_lk, self.reco_hk,...
                                self.lk, self.hk, self.volat_bot);
                        end
                            self.nmin= exit{1};
                        self.reco_top = exit{2};
                        self.reco_bot = exit{3};
                        self.Dconc = self.concFeed.*... 
                            self.reco_top;
                        self.Dconc = self.Dconc ./ sum(self.Dconc);
                        self.Bconc = self.concFeed.*... 
                            self.reco_bot;
                        self.Bconc = self.Bconc ./ sum(self.Bconc);
                        self.Dflow = sum(self.concFeed.*... 
                            self.reco_top.*self.flujo);
                        self.Bflow = sum(self.concFeed.* ...
                            self.reco_bot.*self.flujo);
                    end
                end
            else
                if self.gdlibertad == 0
                    if ~islogical(self.concFeed) && ~islogical(self.reco_lk) &&...
                        ~islogical(self.reco_hk)
                        exit = distribute(self, self.reco_lk, self.reco_hk,...
                            self.lk, self.hk, self.volat_bot);
                        self.nmin= exit{1};
                        self.reco_top = exit{2};
                        self.reco_bot = exit{3};
                        self.Dconc = self.concFeed.*... 
                            self.reco_top;
                        self.Dconc = self.Dconc ./ sum(self.Dconc);
                        self.Bconc = self.concFeed.*... 
                            self.reco_bot;
                        self.Bconc = self.Bconc ./ sum(self.Bconc);
                        self.Dflow = sum(self.concFeed.*... 
                            self.reco_top.*self.flujo);
                        self.Bflow = sum(self.concFeed.* ...
                            self.reco_bot.*self.flujo);
                    end
                end
            end
            salida = self;
        end
        function salida = set_fraccion_vapor(self, fraccionvapor)
            self.beta_feed = fraccionvapor;
            salida = self;
        end
        function tantearecuperaciones(self)
            if length(self.lk) > 1
                recup_lk = self.reco_lk;               
                self.reco_lk = recup_lk - ((1 - recup_lk) / length(self.lk)); %Se identifica la recuperación del clave más pesado
                exit = self.distribute(...
                self.reco_lk, self.reco_hk, max(self.lk), ...
                    min(self.hk));
                self.nmin= exit{1};
                self.reco_top = exit{2};
                self.reco_bot = exit{3};            
                self.Dflow = sum(self.concFeed.* ... 
                self.reco_top.*self.flujo);
                self.Bflow = sum(self.concFeed.* ...
                self.reco_bot.*self.flujo);  
                livianos_feed = self.flujo.*self.concFeed(min(self.lk):max(self.lk));
                dest_calc_livianos = livianos_feed.*self.reco_top(min(self.lk):max(self.lk));
                error = sum(dest_calc_livianos)./sum((livianos_feed)) - recup_lk;
                contador = 0;              
                error1 = error;
                delta = abs(3.5*error);
                recup1 = self.reco_lk;
                self.reco_lk = self.reco_lk - sign(error).*delta;
                while abs(error)> 5e-5 || error  > 0
                    contador = contador + 1;                   
                    exit = self.distribute(...
                        self.reco_lk, self.reco_hk, max(self.lk), ...
                        min(self.hk));
                    self.nmin= exit{1};
                    self.reco_top = exit{2};
                    self.reco_bot = exit{3};            
                    self.Dflow = sum(self.concFeed.* ... 
                        self.reco_top.*self.flujo);
                    self.Bflow = sum(self.concFeed.* ...
                        self.reco_bot.*self.flujo);  
                    livianos_feed = self.flujo.*self.concFeed(min(self.lk):max(self.lk));
                    dest_calc_livianos = livianos_feed.*self.reco_top(min(self.lk):max(self.lk));
                    error = sum(dest_calc_livianos)./sum(livianos_feed) - recup_lk;                    
                    delta = abs(3.5*error);
                    contador = contador + 1;
                    if contador>0 && sign(error)~= sign(error1)                                
                        self.reco_lk = (recup1-self.reco_lk)*((-error)/(error1-error))+self.reco_lk;
                        recup1 = self.reco_lk;
                        error1 = error;
                    else 
                        self.reco_lk = self.reco_lk - sign(error).*delta;
                    end
                    self.rmin = self.get_underwood(self.reco_top, self.reco_bot, self.beta_feed, self.volat_top, self.volat_bot);
                end
                self.lk = max(self.lk);
            end
            if length(self.hk) > 1
                recup_hk = self.reco_hk;                
                self.reco_hk = recup_hk - ((1 - recup_hk) / length(self.hk)); %Se identifica la recuperación del clave más liviano
                exit = self.distribute(...
                self.reco_lk, self.reco_hk, max(self.lk), ...
                    min(self.hk));
                self.nmin= exit{1};
                self.reco_top = exit{2};
                self.reco_bot = exit{3};            
                self.Dflow = sum(self.concFeed.* ... 
                self.reco_top.*self.flujo);
                self.Bflow = sum(self.concFeed.* ...
                self.reco_bot.*self.flujo);  
                pesados_feed = self.flujo.*self.concFeed(min(self.hk):max(self.hk));
                resid_calc_pesados = pesados_feed.*self.reco_bot(min(self.hk):max(self.hk));
                error = sum(resid_calc_pesados)./sum((pesados_feed)) - recup_hk;                
                delta = abs(3.5*error);
                contador = 0;              
                error1 = error;
                recup1 = self.reco_hk;
                self.reco_hk = self.reco_hk - sign(error).*delta; 
                while abs(error)> 5e-5 || error  > 0
                    contador = contador + 1;
                    exit = self.distribute(...
                        self.reco_lk, self.reco_hk, max(self.lk), ...
                        min(self.hk));
                    self.nmin= exit{1};
                    self.reco_top = exit{2};
                    self.reco_bot = exit{3};            
                    self.Dflow = sum(self.concFeed.* ... 
                        self.reco_top.*self.flujo);
                    self.Bflow = sum(self.concFeed.* ...
                        self.reco_bot.*self.flujo);  
                    pesados_feed = self.flujo.*self.concFeed(min(self.hk):max(self.hk));
                    resid_calc_pesados = pesados_feed.*self.reco_bot(min(self.hk):max(self.hk));
                    error = sum(resid_calc_pesados)./sum(pesados_feed) - recup_hk;                    
                    delta = abs(3.5*error);
                    contador = contador + 1;
                    if contador>0 && sign(error)~= sign(error1)                                
                        self.reco_hk = (recup1-self.reco_hk)*((-error)/(error1-error))+self.reco_hk;
                        recup1 = self.reco_hk;
                        error1 = error;
                    else 
                        self.reco_hk = self.reco_hk - sign(error).*delta;
                    end
                end
                self.rmin = self.get_underwood(self.reco_top, self.reco_bot, self.beta_feed, self.volat_top, self.volat_bot);
                self.hk = min(self.hk);
            end
        end
        function salida = westerberg(self)
            if self.gdlibertad == 0;
                alfa = zeros(1, self.num_sust);
                for i=1:self.num_sust
                    alfa(i) = (self.volat_top(i).*self.volat_bot(i)).^0.5/((self.volat_top(self.hk).*self.volat_bot(self.hk)).^0.5);
                end
                WesterbergN = zeros(1,2);
                WNmin = zeros(1, 2);
                WRmin = zeros(1, 2);
                WNmin(1, 1) = 12.29 / (((alfa(self.lk)-1).^0.645)*(1-self.reco_lk));
                WNmin(1, 2) = 12.29 / (((alfa(self.lk)-1).^0.645)*(1-self.reco_hk));
                WRmin(1, 1) = 1.38 / (((alfa(self.lk)-1).^0.717)*(1-self.reco_lk));
                WRmin(1, 2) = 1.38 / (((alfa(self.lk)-1).^0.717)*(1-self.reco_hk));
                WesterbergN(1) = max(WNmin);
                WesterbergN(2) = min(WNmin);
                WesterbergR(1) = max(WRmin);
                WesterbergR(2) = min(WRmin);
                self.WesterbergNmin = 0.8*WesterbergN(1) + 0.2*WesterbergN(2);
                self.WesterbergRmin = 0.8*WesterbergR(1) + 0.2*WesterbergR(2);
            end
            salida = self;
        end
        function salida = get_smoker(self, reflux)
            alfa = zeros(1, self.num_sust);
            alfa(self.lk:self.hk) = (self.volat(self.lk:self.hk).*self.volat_top(self.lk:self.hk)).^0.5./(self.volat(self.hk).*self.volat_top(self.hk)).^0.5;
            
            if reflux > self.rmin
                self.reflujo = reflux;
            else
                error('El reflujo es menor que el mínimo');
            end
            ratioLV= reflux./(reflux + 1);
            ratioDV = 1./(reflux + 1);
            %Smoker supone una línea recta de operación, que es mx+b,
            %entonces m es el ratioLV en la sección de enriquecimiento
            a2 = ratioLV.*(alfa(self.lk)-1);
            xDlk = self.reco_lk.*self.flujo.*self.concFeed(self.lk)./(self.reco_lk.*self.flujo.*self.concFeed(self.lk) + (1-self.reco_hk).*self.flujo.*self.concFeed(self.hk));
            
            a1 = (ratioLV + ratioDV.*xDlk.*(alfa(self.lk) - 1) - alfa(self.lk));
            a0 = ratioDV.*xDlk;
            k = roots([a2, a1, a0]);
            k = k( k > 0 );
            k = k ( k < 1);
            c = (1 + (alfa(self.lk) - 1).*k);
            m = ratioLV;
            log_denom = log(alfa(self.lk)./(c.^2.*m));
            b = ratioDV.*xDlk;
            xprima0 = xDlk - k;
            if abs(self.q - 1) < 1e-3
                zf = self.concFeed(self.lk)./(self.concFeed(self.lk) + self.concFeed(self.hk));
                xpriman = zf - k;
            else 
                zf =  self.concFeed(self.lk)./(self.concFeed(self.lk) + self.concFeed(self.hk));
                zf = (-b - ((zf)/(self.q - 1)))./(-(self.q)/(self.q-1) + m);
                xpriman = zf - k;
            end
            log_num = log(((xprima0).*(1 - xpriman.*((c.*(alfa(self.lk)-1).*m)./(alfa(self.lk)-c.^2.*m))))./(xpriman.*(1 - xprima0.*((c.*(alfa(self.lk)-1).*m)./(alfa(self.lk)-c.^2.*m)))));
            self.SmokerR = log_num./log_denom;
            alfa(self.lk:self.hk) = (self.volat(self.lk:self.hk).*self.volat_bot(self.lk:self.hk)).^0.5./(self.volat(self.hk).*self.volat_bot(self.hk)).^0.5;
            %Seccion de despojo 
            xBlk = (1 - self.reco_lk).*self.flujo.*self.concFeed(self.lk)./((1 - self.reco_lk).*self.flujo.*self.concFeed(self.lk) + (self.reco_hk).*self.flujo.*self.concFeed(self.hk));
            m = (reflux.*zf + xDlk - (reflux + 1).*xBlk)./((reflux + 1).*(zf - xBlk));
            b = (xBlk.*(zf - xDlk))./(((reflux + 1).*(zf - xBlk)));
            ratioLV= m;
            ratioDV = b;
            xpriman = zf - k;
            a2 = ratioLV.*(alfa(self.lk)-1);
            a1 = (ratioLV + ratioDV.*(alfa(self.lk) - 1) - alfa(self.lk));
            a0 = ratioDV;
            k = roots([a2, a1, a0]);
            k = k( k > 0 );
            k = k ( k < 1);
            c = (1 + (alfa(self.lk) - 1).*k);
            log_denom = log(alfa(self.lk)./(c.^2.*m));
            xprima0 = zf - k;
            xpriman = xBlk - k;
            log_num = log(((xprima0).*(1 - xpriman.*((c.*(alfa(self.lk)-1).*m)./(alfa(self.lk)-c.^2.*m))))./(xpriman.*(1 - xprima0.*((c.*(alfa(self.lk)-1).*m)./(alfa(self.lk)-c.^2.*m)))));
            self.SmokerS = log_num./(log_denom);
            self.Smoker = self.SmokerR + self.SmokerS;
            salida = self;
        end
        function salida = reset(self)

            self.id = '';
            self.nreal = false;
            self.feedtray = false;
            self.num_sust = false;
            self.flujo = 1;
            self.units_flujo = 'kg-mol / h';
            self.Dflow = false;
            self.units_Dflow = 'kg-mol / h';
            self.Dconc = false;
            self.Bflow = false;
            self.units_Bflow = 'kg-mol / h';
            self.Bconc = false;
            self.reco_lk = false;
            self.reco_hk = false;
            self.reco_top = false;
            self.reco_bot = false;
            self.cfeed = false;
            self.cdest = false;
            self.cbot = false;
            self.pres_drop_top = 0;
            self.units_press_drop_top = 'kPa';
            self.pres_drop_bot = 0;
            self.units_press_drop_bot = 'kPa';
            self.lk = false;
            self.hk = false;
            self.nmin = false;
            self.rmin = false;
            self.reflujo = false;
            self.nteo = false;
            self.kirkbride = false;
            self.gdlibertad = 3;
            self.volat = false;
            self.volat_top = false;
            self.volat_bot = false;
            self.concFeed = false;
            self.concliqtope = false;
            self.concvapfondo = false;
            self.beta_feed = false;
            self.q = false;
            self.MEdE = false;
            self.WesterbergNmin = false;
            self.WesterbergRmin = false;
            self.SmokerR = false;
            self.SmokerS = false;
            self.Smoker = false;
            salida = self;
        end
    end
end