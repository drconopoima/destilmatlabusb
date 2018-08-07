classdef LewisMatheson < handle
    %LewisMatheson Calcula una torre de destilación por el método riguroso
    %de Lewis y Matheson cuando el usuario ha aportado suficientes
    %grados de libertad
    %   Luis Jesús Díaz Manzo
    
    properties
        id = '';
        platos_reales = false;
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
        reco_bottom = false;
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
        reflujofondo = false;
        nteo = false;
        kirkbride = false;
        gdlibertad = 1;
        volat = false;
        volat_top = false;
        volat_bot = false;
        concFeed = false;
        concliqtope = false;
        concvapfondo = false;
        beta_feed = false;
        q = false;
        heat_loss = 0;
        ratioLV = false;
        ratioLsVs = false;
        ratioLG = false;
        MEdE = false;
        TPlatoPlato = false;
        concLiqPlato = false;
        PlatosTot = false;
        etapaoptima = false;
        concVapPlato = false;
        Qc = false;
        Qb = false;
        starttop = true;
        startbottom = false;
        lk_suge = false;
        hk_suge = false;
end
    
    methods
        function self = LewisMatheson(alim, lightk, recup_top, heavyk, recup_bot, startthere, dest, bot, nombre)
            self.gdlibertad = 3;
            if nargin > 0;
                if ~isempty(alim) && isa(alim, 'Corriente')
                    self.num_sust = alim(1).num_sust;
                    self.cfeed = alim;
                    self.concFeed = alim.conc;
                    
                    self.beta_feed = 1 - alim.q;
                    self.q = alim.q;
                    self.flujo = alim.molF;
                    self.gdlibertad = self.gdlibertad - 1;
                    self.MEdE = alim(1).MEdE;
                end
            end
            if nargin > 2   
                if ~isempty(recup_top) && ~isempty(lightk)
                    self.set_reco_top(lightk, recup_top);
                end
            end
            if nargin > 5
                if strcmpi(startthere, 'y') 
                    self.starttop = false;
                    self.startbottom = true;
                end
            end
            if nargin >  4
                if ~isempty(recup_bot) && ~isempty(heavyk)
                    self.set_reco_bot(heavyk, recup_bot);
                end
            end
            if nargin > 6
                if ~isempty(dest) && isa(dest, 'Corriente')
                    self.cdest = dest;
                end
            else 
                    self.cdest = Corriente.empty(0,1);
            end
            if nargin > 7
                if ~isempty(bot) && isa(bot, 'Corriente')
					self.cbot = bot;
                end
            else 
                self.cbot = Corriente.empty(0,1);
            end
            if nargin > 8
                if ~isempty(nombre) && isa(nombre, 'char')
                    self.id = nombre;
                end
            end   
            %recuperacion tope clave liviano, recuperación fondo pesado, 
            %y definición de una corriente de alimentación o un conjunto de
            %volates de tope y/o fondo o constantes para la torre.
            %Ademas de las composiciones arbitrarias para los componentes
%             self = self.checkdistribute();
        end
        function salida = set_pressure_drop_top(self, pressure_drop, unit)
            %Establece la caída de presión por plato de la región de
            %rectificación
                if ~islogical(pressure_drop)
                    if nargin == 2
                        unit = self.units_press_drop_top;
                    end
                    pressure_drop = unit_pressure(pressure_drop, unit, self.units_press_drop_top);
                    self.pres_drop_top = pressure_drop; 
                    salida = self;
                end
        end
        function salida = set_pressure_drop_bot(self, pressure_drop, unit)
            %Establece la caída de presión por plato de la región de
            %despojo
            if ~islogical(pressure_drop)
                if nargin == 2
                        unit = self.units_press_drop_bot;
                end
                pressure_drop = unit_pressure(pressure_drop, unit, self.units_press_drop_bot);
                self.pres_drop_bot = pressure_drop;
                salida = self;
            end
        end
        function salida = set_cfeed(self, streamfeed)
            i = length(self.cfeed);
            %Pueden haber varias alimentaciones, cfeed almacena
            %en un array todas las alimentaciones dadas.
            self.cfeed(i+1) = streamfeed;
            if self.cfeed.gdlibertad == 0
                self.num_sust = streamfeed.num_sust;
                self.gdlibertad = self.gdlibertad - 1;
                %Escribo el flash de la alimentación en un output 
                fid = fopen('Output_LM.dat','w+');
                fprintf(fid,'La alimentación esta a condicion: flujo=%f, beta=%f\n   zi       yi       xi\n', self.flujo, streamfeed.fraccionvapor);
                fprintf(1,'La alimentación esta a condicion: flujo=%f, beta=%f\n   zi       yi       xi\n', self.flujo, streamfeed.fraccionvapor);
                fclose('all');          
                matrix_output = [streamfeed.conc; streamfeed.conc_vap; streamfeed.conc_liq];
                self.printmfract(matrix_output);
                self.concFeed = streamfeed.conc;
                self.volat = streamfeed.K;
                self.beta_feed = 1 - streamfeed.q;
                self.flujo = streamfeed.molF;
                self.MEdE = streamfeed(1).MEdE;
            end
            salida = self;
        end
        function salida = set_cdest(self, currentdist)
            i= length(self.cdest);
            %Pudiera tener hasta 2 corrientes de destilado, una vaporizada
            %y otra a condensado líquido.
            self.cdest(i+1) = currentdist;
            salida = self;
        end
        function salida = set_cfondo(self, current_bottom)
            self.cfondo = current_bottom; 
            salida = self;
        end   
        function salida = set_reco_top( self, clave_liviano, recup_lk)
            if ischar(clave_liviano) || length(clave_liviano)==1
                if ~islogical(recup_lk) || ~islogical(clave_liviano)
                    if islogical(self.lk)
                        %se especificó la recuperación del clave liviano y
                        %por tanto se disminuye los grados de libertad.
                        self.gdlibertad = self.gdlibertad - 1;
                    end
                    if ischar(clave_liviano)
                        %si el usuario introdujo el clave como un string
                        %con el id del clave deseado.
                        for i = 1:self.num_sust
                            if strcmpi(clave_liviano, self.cfeed.compuestos(...
                                    i).id)
                                self.lk = i; %se identifica ubicación del clave
                            end
                        end
                    else                        
                        %Si usuario introduce valor numérico que refleja la
                        %ubicación del clave liviano entre los
                        %ides de compuestos.
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
                                if strcmpi(clave_liviano{j}, self.cfeed.compuestos(...
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
            if ~isempty(self.cfeed)
                if self.gdlibertad == 0 && self.cfeed.gdlibertad == 0
                    self.concFeed = self.cfeed.conc;
                    if length(self.hk) == 1 && length(self.lk) == 1
                        exit = self.distribute(...
                            self.reco_lk, self.reco_hk, max(self.lk), ...
                            min(self.hk));
                        self.nmin= exit{1};
                        self.reco_top = exit{2};
                        self.reco_bottom = exit{3};            
                        self.Dflow = sum(self.cfeed.conc.* ... 
                            self.reco_top.*self.cfeed.flujo);
                        self.Bflow = sum(self.cfeed.conc.* ...
                            self.reco_bottom.*self.cfeed.flujo);
                        self.Dconc = self.flujo.*self.reco_top.*self.concFeed...
                            ./ self.Dflow;                
                        self.Bconc = self.flujo.*self.reco_bottom.*self.concFeed...
                            ./ self.Bflow;
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
                self.reco_bottom = exit{3};            
                self.Dflow = sum(self.concFeed.* ... 
                    self.reco_top.*self.flujo);
                self.Bflow = sum(self.concFeed.* ...
                    self.reco_bottom.*self.flujo);
                self.Dconc = self.flujo.*self.reco_top.*self.concFeed...
                    ./ self.Dflow;                
                self.Bconc = self.flujo.*self.reco_bottom.*self.concFeed...
                    ./ self.Bflow;
            end   
        salida = self;
        end
        function salida = set_reco_bot( self, clave_pesado, recup_hk)    
            if ischar(clave_pesado) || length(clave_pesado)==1
                if ~islogical(recup_hk) || ~islogical(clave_pesado)
                    if islogical(self.hk)
                        
                        self.gdlibertad = self.gdlibertad - 1;
                    end
                    if ischar(clave_pesado)
                        for i = 1:self.num_sust
                            if strcmp(clave_pesado, self.cfeed.compuestos(...
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
                                if strcmpi(clave_pesado{j}, self.cfeed.compuestos(...
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
                    self.reco_bottom = false;
                    self.hk = false;
                    self.gdlibertad = self.gdlibertad + 1;
                end
            end
            if ~isempty(self.cfeed)
                if self.gdlibertad == 0 && self.cfeed.gdlibertad == 0
                    self.concFeed = self.cfeed.conc;               
                    if length(self.hk) == 1 && length(self.lk) == 1
                        exit = self.distribute(...
                            self.reco_lk, self.reco_hk, self.lk, ...
                            self.hk);
                        self.nmin= exit{1};
                        self.reco_top = exit{2};
                        self.reco_bottom = exit{3};            
                        self.Dflow = sum(self.concFeed.* ... 
                            self.reco_top.*self.flujo);
                        self.Bflow = sum(self.concFeed.* ...
                            self.reco_bottom.*self.flujo);                        
                        self.Dconc = self.flujo.*self.reco_top.*self.concFeed...
                            ./ self.Dflow;                
                        self.Bconc = self.flujo.*self.reco_bottom.*self.concFeed...
                            ./ self.Bflow;
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
                self.reco_bottom = exit{3};            
                self.Dflow = sum(self.concFeed.* ... 
                    self.reco_top.*self.flujo);
                self.Bflow = sum(self.concFeed.* ...
                    self.reco_bottom.*self.flujo);                
                self.Dconc = self.flujo.*self.reco_top.*self.concFeed...
                    ./ self.Dflow;                
                self.Bconc = self.flujo.*self.reco_bottom.*self.concFeed...
                    ./ self.Bflow;
            end
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
                self.reco_bottom = exit{3};            
                self.Dflow = sum(self.cfeed.conc.* ... 
                self.reco_top.*self.cfeed.flujo);
                self.Bflow = sum(self.cfeed.conc.* ...
                self.reco_bottom.*self.cfeed.flujo);  
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
                    self.reco_bottom = exit{3};            
                    self.Dflow = sum(self.cfeed.conc.* ... 
                       self.reco_top.*self.cfeed.flujo);
                    self.Bflow = sum(self.cfeed.conc.* ...
                       self.reco_bottom.*self.cfeed.flujo);  
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
                self.reco_bottom = exit{3};            
                self.Dflow = sum(self.cfeed.conc.* ... 
                self.reco_top.*self.cfeed.flujo);
                self.Bflow = sum(self.cfeed.conc.* ...
                self.reco_bottom.*self.cfeed.flujo);  
                pesados_feed = self.flujo.*self.concFeed(min(self.hk):max(self.hk));
                resid_calc_pesados = pesados_feed.*self.reco_bottom(min(self.hk):max(self.hk));
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
                    self.reco_bottom = exit{3};            
                    self.Dflow = sum(self.cfeed.conc.* ... 
                        self.reco_top.*self.cfeed.flujo);
                    self.Bflow = sum(self.cfeed.conc.* ...
                        self.reco_bottom.*self.cfeed.flujo);  
                    pesados_feed = self.flujo.*self.concFeed(min(self.hk):max(self.hk));
                    resid_calc_pesados = pesados_feed.*self.reco_bottom(min(self.hk):max(self.hk));
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
                self.hk = min(self.hk);
            end
        end
        function salida = distribute(self, recuplk,...
                recuphk, light, heavy, volatil_top, volatil_bottom)
            if nargin ==5
                coefKhk = self.cfeed(1).K(heavy);
                %coefKlk y coefKhk verifican que haya suministrado el liviano
                %y el pesado en su orden adecuado de volat
                coefKlk = self.cfeed(1).K(light);
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
                    if ~isa(self.volat_bot, 'double') && isa(self.volat_top, 'double');
                        self.volat = self.volat_top;
                    %self.volat = self.cfeed(1).conc_vap ./ self.cfeed(1).conc_liq;
                    elseif ~isa(self.volat_top, 'double') && isa(self.volat_bot, 'double')
                        self.volat = self.volat_bot;
                    elseif isa(self.volat_top, 'double') && isa(self.volat_bot, 'double')
                        self.volat = (self.volat_top.*self.volat_bot).^0.5;
                    else
                        self.volat = self.cfeed(1).conc_vap ./ self.cfeed(1).conc_liq;
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
                    alfa = self.cfeed(1).K ./ self.cfeed(1).K(self.hk);
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
                    %Si hay comp de volat intermedia, se estima
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
                if ~isempty(self.cfeed)
                    if self.cfeed(1).gdlibertad == 0
                        self.Dflow = sum(self.flujo.*...
                            self.concFeed(1,:).*recup_top);
                        self.Bflow = sum(self.flujo.*...
                            self.concFeed(1, :).*recup_bottom);
                        xD = (self.flujo.*self.concFeed(1, :)...
                            .*recup_top)./(self.Dflow);
                        xB = (self.flujo.*self.concFeed(1, :)...
                            .*recup_bottom)./(self.Bflow);
                        [~, xtope, K1, flag1] = self.MEdE.DewT(self.cfeed(1).P, self.cfeed(1).mezcla, xD);
                        [~, yfondo, K2, flag2] = self.MEdE.BubbleT(self.cfeed(1).P, self.cfeed(1).mezcla, xB); 
                        KD = K1;
                        KB = K2;            
                        self.volat_top = KD ./ KD(heavy);
                        self.volat_bot = KB ./ KB(heavy);
                        self.volat = (self.volat_top.*self.volat_bot).^0.5;
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
            
            %Flujo de destilado y de fondo
            self.Dflow=sum(self.flujo.*self.concFeed.*recup_top);
            self.Bflow=sum(self.flujo.*self.concFeed.*recup_bottom);            
            fid = fopen('Output_LM.dat','a+');
            fprintf(fid,'La recuperación supuesta es: D=%f B=%f\n  Recup. D  Recup. B\n',self.Dflow, self.Bflow);     
            fprintf(1,'La recuperación supuesta es: D=%f B=%f\n  Recup. D  Recup. B\n',self.Dflow, self.Bflow);            
            fclose('all');
            matrix_output = [recup_top; recup_bottom];
            self.printmfract(matrix_output);
            %Concentraciones en fracciones molares de destilado y fondo            
            fid = fopen('Output_LM.dat','a+');
            fprintf(fid,'Las fracciones molares son: \n     yD        xB\n');
            fprintf(1, 'Las fracciones molares son: \n     yD        xB\n');
            fclose('all');    
            matrix_output = [self.Dconc; self.Bconc];
            self.printmfract(matrix_output);
            self.reco_top = recup_top;
            self.reco_bottom = recup_bottom;
            [self.rmin]= self.get_underwood(self.reco_top, self.reco_bottom);
            fid = fopen('Output_LM.dat','a+');
            fprintf(fid,'El número mínimo de etapas es Nmin=%f\ny el reflujo mínimo es Rmin=%f\n',mintrays,self.rmin);
            fprintf(1,'El número mínimo de etapas es Nmin=%f\ny el reflujo mínimo es Rmin=%f\n',mintrays,self.rmin);
            fclose('all');
            salida = {mintrays, recup_top, recup_bottom};
            
        end
        function printmfract(self, matrix_horizontal) %imprime vectores en el archivo 'Output_LM.dat' en 1,2,3,4 columnas
            %Útil para imprimer recuperaciones de Fenske o las fracciones molares
            fid = fopen('Output_LM.dat', 'a+');
            area = (size(matrix_horizontal)); %dimensiones de la matrix dada
            height = area(1);
            wid = area(2);
            if height == 2
                for i=1:wid
                    fprintf(fid,'  %f  %f \n',matrix_horizontal(1,i), matrix_horizontal(2,i));
                    fprintf(1,'  %f  %f \n',matrix_horizontal(1,i), matrix_horizontal(2,i));
                end
            elseif height == 1
                for i=1:wid
                    fprintf(fid,'  %f \n',matrix_horizontal(1,i));
                    fprintf(1,'  %f \n',matrix_horizontal(1,i));
                end
            
            elseif height == 3                
                for i=1:wid
                    fprintf(fid,'  %f  %f  %f  \n',matrix_horizontal(1,i), matrix_horizontal(2,i), matrix_horizontal(3,i));
                    fprintf(1,'  %f  %f  %f  \n',matrix_horizontal(1,i), matrix_horizontal(2,i), matrix_horizontal(3,i));
                end
            elseif height == 4
                
                for i=1:wid
                    fprintf(fid,'  %f  %f  %f  %f  \n',matrix_horizontal(1,i), matrix_horizontal(2,i), matrix_horizontal(3,i), matrix_horizontal(4,i));
                    fprintf(1,'  %f  %f  %f  %f  \n',matrix_horizontal(1,i), matrix_horizontal(2,i), matrix_horizontal(3,i), matrix_horizontal(4,i));
                end
            end        
            fclose('all');
        end
        function [num_min_trays, recupdist, recupresidue] = getfenske(...
                self, recup_dist, recup_resid, light,...
                heavy, volatil_top, volatil_bottom)
            %Resuelve Fenske para hallar el número de etapas mínimo y la
            %distribución de los componentes a los productos
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
                    alfa = (self.volat_bot.*self.volat_top).^0.5 ./ (self.volat_bot(self.hk).*self.volat_top(self.hk))^0.5;
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
            if self.lk~=1 && self.hk ~=self.num_sust
                for i = self.lk+1:self.num_sust
                    [minreflux(i), Dflow(i)] = underwood(self.concFeed, self.flujo,...
                        sup_recup_d, sup_recup_b, alfa, light, i, self.beta_feed);
                end
                [minreflux, indice] = sort(minreflux, 'descend');
                if indice(1) ~= heavy
                    fprintf('Se sugiere cambiar el componente clave pesado por el HK = %d \n', indice(1));
                    self.hk_suge = indice(1);
                end
                for i = 1:indice(1) - 1
                    [minreflux(i), Dflow(i)] = underwood(self.concFeed, self.flujo,...
                        sup_recup_d, sup_recup_b, alfa, i, indice(1), self.beta_feed);
                end
                [minreflux, indice] = sort(minreflux, 'descend');
                if indice(1) ~= light
                    fprintf('Se sugiere cambiar el componente clave liviano por el LK = %d \n', indice(1));
                    self.lk_suge = indice(1);
                end
                minreflux = minreflux(1);
                Dflow = Dflow(indice(1));
            else
                [minreflux, Dflow] = underwood(self.concFeed, self.flujo,...
                        sup_recup_d, sup_recup_b, alfa, light, heavy, self.beta_feed);
            end
        end
        function salida = reflujoreal(self, reflujodato)
            self.reflujo = reflujodato;
            fid = fopen('Output_LM.dat','a+');    
            fprintf(fid, 'Con un reflujo de: %f\n',self.reflujo);
            fprintf(1, 'Con un reflujo de: %f\n',self.reflujo);
            fclose('all');          
            self.platosrect();
            salida = self;            
        end
        function salida = reflux2rmin(self, reflujodato)
            self.reflujo = reflujodato*self.rmin;
            fid = fopen('Output_LM.dat','a+');    
            fprintf(fid, 'Con un reflujo de: %f\n',self.reflujo);
            fprintf(1, 'Con un reflujo de: %f\n',self.reflujo);
            fclose('all');
            self.platosrect();
            salida = self;
        end
        function salida = platosrect(self)
            num_sust = self.num_sust;
			P = self.cfeed(1).P;
            q = 1 - self.beta_feed;
            LKey = self.lk;
            HKey = self.hk;
            Flow = self.flujo;
            alfa = (self.volat_top.*self.volat_bot).^0.5;
            alfa = (alfa)./(alfa(HKey));
            recup_top = zeros(1, num_sust);
            recup_bottom = zeros(1, num_sust);
            for i = 1 : num_sust
                if alfa(i) > alfa(self.lk)
                    %Distribución de más ligeros que LK en base a
                    %volat relativa
                    recup_top(i) = (self.reco_lk+(1 - self.reco_lk)*(1 -...
                        1/(alfa(i))^(10)));
                    recup_bottom(i) = (1 - recup_top(i));
                elseif i == self.lk
                    %Distribución del clave liviano
                    recup_bottom(i) = 1 - self.reco_lk;
                    recup_top(i) = 1 - recup_bottom(i);
                elseif i == self.hk
                    %Distribución del clave pesado
                    recup_top(i)= 1 - self.reco_hk;
                    recup_bottom(i) = 1 - recup_top(i);
                elseif alfa(i) < 1
                    %Distribución de más pesador que HK en base a
                    %volat relativa
                    recup_top(i) = (1 - self.reco_hk)*(alfa(i)^(10));
                    recup_bottom(i) = 1 - recup_top(i);
                else
                    %Si hay comp de volat intermedia, se estima
                    %su distribución como 0.5 y se recalculará con Fenske y
                    %Underwood
                    recup_top(i) = self.reco_top(i);
                    recup_bottom(i) = self.reco_bottom(i);
                end    
            end
            if any(recup_top == 0) || any(recup_top == 1)
                for i =1 : num_sust
                    if recup_top(i) == 0
                        recup_top(i) = 1e-15;
                    elseif recup_top(i) == 1
                        recup_top(i) = 1-1e-15;
                    end
                end
            end
            if any(recup_bottom == 0) || any(recup_bottom == 1)
                for i =1 : self.num_sust
                    if recup_bottom(i) == 0
                        recup_bottom(i) = 1e-15;
                    elseif recup_bottom(i) == 1
                        recup_bottom(i) = 1-1e-15;
                    end
                end
            end
            %Temperatura TDvapor la del destilado antes del condensador
            %total
            self.reco_top = recup_top;
            self.reco_bottom = recup_bottom;
            self.Dconc = (self.reco_top.*self.concFeed.*self.flujo)./sum((self.reco_top.*self.concFeed.*self.flujo));
            self.Bconc = (self.reco_bottom.*self.concFeed.*self.flujo)./sum((self.reco_bottom.*self.concFeed.*self.flujo));
            self.Dflow = sum((self.reco_top.*self.concFeed.*self.flujo));
            self.Bflow = sum((self.reco_bottom.*self.concFeed.*self.flujo));
            comp = self.cfeed(1).comp;
            kij =  self.cfeed(1).mezcla.kij;
            D = self.Dflow;
            if isa(self.reflujo, 'double') && self.reflujo > self.rmin
                reflux = self.reflujo;
            elseif self.reflujo < self.rmin
                error('El valor de reflujo es menor al reflujo mínimo')
            elseif  isa(self.reflujofondo, 'double') && self.reflujo > self.rmin
                refluxfondo = self.reflujofondo;
            end
            L = D*(reflux);
            V = L + D;
            if self.starttop
                if self.lk ~= 1 || self.hk ~= self.num_sust
                    ratioLV = reflux/(reflux + 1);
                    ratioDV = 1/(reflux + 1);
                    fid = fopen('Output_LM.dat','a+');    
                    fprintf(fid,'Para el plato de tope, los flujos molares de liquido y vapor son: L=%f V=%f\n', L , V );
                     fprintf(1,'Para el plato de tope, los flujos molares de liquido y vapor son: L=%f V=%f\n', L , V );
                    fprintf(fid, 'Con ello, se calcula las L/V = %f y D/V = %f\n\n', ratioLV, ratioDV);
                    fprintf(1, 'Con ello, se calcula las L/V = %f y D/V = %f\n\n', ratioLV, ratioDV);

                    Y = self.Dconc;
                    xD = self.Dconc;
                    MezclaPlato_i = Mezcla(comp, Y, kij);
                    [Tr, sup_x, ~, flag1]= self.MEdE.DewT(P, MezclaPlato_i);
                    [TDvapor, X1, K1, flag2, valI, iter]= self.MEdE.DewT(P, MezclaPlato_i, sup_x, [], [], Tr);
                    for i = 1:num_sust
                        alfa1(1, i) = K1(i)./K1(HKey);
                    end
                    %Temperatura TL0 es la del líquido que retorna a la columna
                    %y también la del destilado líquido a retirar de la columna
                    T1 = zeros(1, 100);
                    [Tb, sup_y, ~, flag3] = self.MEdE.BubbleT(P, MezclaPlato_i);
                    [Tb, sup_y, KL0, flag4] = self.MEdE.BubbleT(P, MezclaPlato_i, sup_y, [], [], Tb);
                    T1(1) = TDvapor;
                    marcador = 0;
                    beta = self.beta_feed;
                    x = self.concFeed;
                    terminador = (x(LKey) - ((xD(LKey).*(1 - q))/( 1+reflux)))/(x(HKey)- ((xD(HKey).*(1-q))./(reflux + 1)));
                    check_terminador1 = terminador + 1;            
                    fprintf(fid,'En la alimentacion xF(LKey)/xF(HKey)=%f\n', terminador);
                    fprintf(1,'En la alimentacion xF(LKey)/xF(HKey)=%f\n', terminador);
                    fclose('all');
                    num_sust = num_sust;
                    ite = 0;
        %             ratioLV = 0.445;
        %             ratioDV = 0.555;
        %             X1 = [0.00425, 0.0425, 0.2495, 0.6612, 0.0425, 1e-15];
        %             D = 0.38;
        %             Y = [0.03, 0.07, 0.147, 0.13, 0.003, 1e-15];
        %             Y = Y./D;
        %             xD = Y;
        %             MezclaPlato_i = Mezcla(comp, X1, kij);
        %             [Tr, sup_x, ~, flag1]= self.MEdE.DewT(P, MezclaPlato_i);
        %             [TDvapor, fX1, K1, flag2, valI, iter]= self.MEdE.DewT(P, MezclaPlato_i, sup_x, [], [], Tr);
                    while marcador == 0
                        ite=ite+1;
                        for i=1:num_sust
                            Y(ite+1,:)=(ratioLV.*X1(ite,:)+ratioDV.*xD);
                        end
                        MezclaPlato_i(ite+1) = Mezcla(comp, Y(ite+1, :), kij);
                        [Tr, sup_x, ~, flag3]= self.MEdE.DewT(P, MezclaPlato_i(ite+1));
                        [T1(ite+1), X1(ite+1,:), K, flag3]= self.MEdE.DewT(P, MezclaPlato_i(ite+1), sup_x, [], [], Tr);
                        for i = 1:num_sust
                            alfa1(ite+1, i) = K(i)./K(HKey);
                        end
                        check_terminador1 = (X1(ite+1,LKey)/X1(ite+1,HKey));
                        fid = fopen('Output_LM.dat','a+');  
                        fprintf(fid,'En el plato %d desde el condensador T=%f K\n      X       Y\n',ite,T1(ite));  
                        fprintf(1,'En el plato %d desde el condensador T=%f K\n      X       Y\n',ite,T1(ite));  
                        fclose('all');
                        matrix_output = [X1(ite+1,:); Y(ite+1,:)];
                        self.printmfract(matrix_output);
                        if check_terminador1 < terminador
                            PlatosRect = ite; %Cuando llego a alimentacion tengo el num de platos
                            marcador = 1;
                        end
                        if ite>100
                            error('No se consiguió alcanzar el plato de alimentación en 100 iteraciones');
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Condicion termica de la alimentacion
                    Ls=(1-beta)*Flow+L; %Liquido en la zona de agotamiento
                    Vs=Ls-self.Bflow; %Vapor en la zona de agotamiento
                    ratioLsVs = Vs/Ls;
                    ratioLsB = self.Bflow/Ls;
                    T2 = zeros(1, 100);
                    xB = self.Bconc;     
                    LiquidoW = Mezcla(comp, xB, kij);
                    [Tb, sup_y, ~, flag1]= self.MEdE.BubbleT(P, LiquidoW);
                    [T2(1), yfondo(1,:) , K1, flag2]= self.MEdE.BubbleT(P, LiquidoW, sup_y, [], [], Tb);     
                    Xs(1,:) = xB;
                    fid = fopen('Output_LM.dat','a+');    
                    fprintf(fid,'Para el calderín, las Vs/Ls = %f y B/Ls = %f \n \n', ratioLsVs , ratioLsB );
                    fprintf(1,'Para el calderín, las Vs/Ls = %f y B/Ls = %f \n \n', ratioLsVs , ratioLsB );
                    fclose('all');
                    %Resolución desde el calderín hasta la alimentación 
                    marcador = 0;
                    ite = 0;
                    LiquidoPlato_i = Mezcla.empty(0, 100);
                    while marcador == 0
                        ite=ite+1;
                        for i=1:num_sust
                            Xs(ite+1,i)=(ratioLsB.*xB(i) + ratioLsVs.*yfondo(ite,i));
                        end
                        LiquidoPlato_i(ite) = Mezcla(comp, Xs(ite+1, :), kij);
                        check_terminador2 = (Xs(ite+1,LKey)/Xs(ite+1,HKey));
                        [Tb, sup_y,~, flag5]= self.MEdE.BubbleT(P, LiquidoPlato_i(ite));
                        [T2(ite+1), yfondo(ite+1,:), K2, flag5]= self.MEdE.BubbleT(P, LiquidoPlato_i(ite), sup_y, [], [], Tb);
                        for i = 1:num_sust
                            alfa2(ite, i) = K2(i)./K2(HKey);
                        end
                        fid = fopen('Output_LM.dat','a+');  
                        fprintf(fid,'Subiendo al plato %d (Ncalderin = 1) T=%f K\n      X       Y\n',ite,T2(ite+1));
                        fprintf(1,'Subiendo al plato %d (Ncalderin = 1) T=%f K\n      X       Y\n',ite,T2(ite+1));
                        fclose('all');
                        matrix_output = [Xs(ite,:); yfondo(ite,:)];
                        self.printmfract(matrix_output);                
                        fid = fopen('Output_LM.dat','a+');  
                        fprintf(fid, '\n');
                        fclose('all');
                        if check_terminador2 > terminador
                            PlatosStrip = ite; %Cuando llego a alimentacion tengo el num de platos
                            %Para no contar la alimentación dos veces, se considera
                            %que el plato de alimentación se contó en la zona
                            %superior y se resta 1, pero como se empezó a contar "ite"
                            %cuando ya se hizo equilibrio en el calderín, entonces
                            %se suma el calderín al número de platos.
                            %Se descarta el valor hallado para el siguiente plato                    
                            T1 = T1(T1~=0);
                            T2 = T2(T2~=0);
                            Xs = Xs(1:end-2,:);
                            T2 = T2(1:end-2);
                            yfondo = yfondo(1:end-1,:);
                            fid = fopen('Output_LM.dat','a+');  
                            fprintf(fid,'El plato de alimentación es el Ns = %f \n',ite);
                            fprintf(1,'El plato de alimentación es el Ns = %f \n',ite);
                            fprintf(fid,'Se descarta el valor conseguido para el plato %f \n', ite+1);
                            fprintf(1,'Se descarta el valor conseguido para el plato %f \n', ite+1);
                            fclose('all');                    
                            marcador = 1;
                        end

                        if ite>100
                            error('No se pudo alcanzar el plato de alimentación en la sección de despojo.')
                        end
                    end
                    %%%%%%%%%%%%%%%% 
                
                else
                        for  kui = 1:2
                            try 
                                self.Bconc = Xs(end,:);
                                self.Bconc = self.Bconc./(sum(self.Bconc));
                                RELIVIANO = (self.reco_bottom.*self.concFeed.*self.flujo)./sum((self.reco_bottom.*self.concFeed.*self.flujo));
                                RELIVIANO = RELIVIANO(self.lk);
                                if RELIVIANO > self.Bconc(self.lk)
                                    RevisionLiviano = (self.concFeed.*self.flujo - self.Bflow.*self.Bconc)./(sum(self.concFeed.*self.flujo - self.Bflow.*self.Bconc));
                                    RevisionLiviano = RevisionLiviano(self.lk);
                                    self.Dconc(self.lk) = RevisionLiviano;
                                    self.Dconc = self.Dconc./(sum(self.Dconc));
                                end
                            catch ME
                                self.Dconc = (self.reco_top.*self.concFeed.*self.flujo)./sum((self.reco_top.*self.concFeed.*self.flujo));
                                self.Bconc = (self.reco_bottom.*self.concFeed.*self.flujo)./sum((self.reco_bottom.*self.concFeed.*self.flujo));
                            end
                            ratioLV = (self.reflujo)./(self.reflujo + 1);
                            ratioDV = 1./(self.reflujo + 1);
                            Ls = self.cfeed.molF.*self.cfeed.q + L;
                            Vs = Ls - self.Bflow;
                            ratioLsVs = Vs./Ls;
                            ratioVsB = self.Bflow ./ Ls;
                            fid = fopen('Output_LM.dat','a+');    
                            fprintf(fid,'Para el plato de tope, los flujos molares de liquido y vapor son: L=%f V=%f\n', L , V );
                            fprintf(1,'Para el plato de tope, los flujos molares de liquido y vapor son: L=%f V=%f\n', L , V );
                            fprintf(fid, 'Con ello, se calcula las relaciones V/L = %f y D/V = %f\n\n', ratioLV, ratioDV);
                            fprintf(1, 'Con ello, se calcula las relaciones V/L = %f y D/V = %f\n\n', ratioLV, ratioDV);
                            Y = self.Dconc;
                            x = self.concFeed;
                            xD = self.Dconc;
                            xB = self.Bconc;
                            MezclaPlato_i = Mezcla(comp, Y, kij);
                            [Tr, sup_x, ~, flag3]= self.MEdE.DewT(P, MezclaPlato_i);
                            [T1, X1, K, flag3]= self.MEdE.DewT(P, MezclaPlato_i, sup_x, [], [], Tr);
                            X1(1,:) = self.Dconc;
                            for i = 1:num_sust
                                alfa1(1, i) = K(i)./K(HKey);
                            end
                            terminador = (x(LKey) - ((xD(LKey).*(1 - q))/( 1+reflux)))/(x(HKey)- ((xD(HKey).*(1-q))./(reflux + 1)));
                            check_terminador1 = terminador + 1;            
                            fprintf(fid,'En la alimentacion xF(LKey)/xF(HKey)=%f\n', terminador);
                            fprintf(1,'En la alimentacion xF(LKey)/xF(HKey)=%f\n', terminador);
                            fclose('all');
                            iter = 0;
                            while check_terminador1 >= terminador
                                iter = iter + 1;
                                MezclaPlato_i(iter) = Mezcla(comp, Y(iter, :), kij);
                                [Tr, sup_x, ~, flag3]= self.MEdE.DewT(P, MezclaPlato_i(iter));
                                [T1(iter+1), X1(iter+1,:), K, flag3]= self.MEdE.DewT(P, MezclaPlato_i(iter), sup_x, [], [], Tr);
                                for i = 1:num_sust
                                    alfa1(iter+1, i) = K(i)./K(HKey);
                                end
                                for i = 1:num_sust
                                    Y(iter+1, i) = ratioLV.*X1(iter+1,i) + ratioDV.*xD(i) ;
                                end
                                check_terminador1 = X1(iter + 1, self.lk)./(X1(iter + 1, self.hk));
                                fid = fopen('Output_LM.dat', 'a+');
                                fprintf(fid, 'El plato %d subiendo del calderín tiene una composición: T=%f K\n      X       Y\n\n', iter, T1(iter+1));
                                fprintf(1, 'El plato %d subiendo del calderín tiene una composición: T=%f K\n      X       Y\n\n', iter, T1(iter+1));
                                fclose('all');
                                matrix_output = [X1(iter+1,:); Y(iter+1,:)];
                                self.printmfract(matrix_output);    
                                
                                if check_terminador1 < terminador
                                    PlatosRect = iter;
                                    fid = fopen('Output_LM.dat', 'a+');
                                    fclose('all');
                                    break
                                end
                            end
                                %%%%%%%% Se continúa operando en la región de
                                %%%%%%%% enriquecimiento
                                yfondo = Y(end, :);
                                terminador = self.flujo.*(1 - self.reco_lk).*self.concFeed(self.hk);
                                check_terminador2 = terminador + 1;        
                                fid = fopen('Output_LM.dat', 'a+');
                                fprintf(fid,'Pasando el plato de alimentación con la fracción de liviano a pesado de %f \n', terminador);
                                fprintf(1,'Pasando el plato de alimentación con la fracción de liviano a pesado de %f \n', terminador);
                                fclose('all');
                                iter = 0;
                                while check_terminador2 >= terminador
                                    iter = iter + 1;
                                    MezclaPlato_i(iter) = Mezcla(comp, yfondo(iter, :), kij);
                                    [Tr, sup_x, ~, flag3]= self.MEdE.DewT(P, MezclaPlato_i(iter));
                                    [T2(iter), Xs(iter,:), K, flag3]= self.MEdE.DewT(P, MezclaPlato_i(iter), sup_x, [], [], Tr);
                                    for i = 1:num_sust
                                        alfa2(iter, i) = K(i)./K(HKey);
                                    end

                                    for i = 1:num_sust
                                        yfondo(iter+1, i) = ratioLsVs.*Xs(iter,i) + ratioVsB.*xB(i) ;
                                    end
                                    if any(Xs(iter,:) < 0)
                                        indice = find(Xs(iter,:) < 0);
                                        for indic = min(indice):max(indice)
                                            Xs(iter, indic) = 1e-15;
                                        end
                                        Xs(iter, :) = Xs(iter, :)./sum(Xs(iter, :));
                                    end
                                    check_terminador2 = self.Bflow.*Xs(iter, self.lk);
                                    fid = fopen('Output_LM.dat', 'a+');
                                    fprintf(fid, 'El plato %d bajando desde la alimentación tiene una composición: T=%f K\n      X       Y\n\n', iter, T2(iter));
                                    fprintf(1, 'El plato %d bajando desde la alimentación tiene una composición: T=%f K\n      X       Y\n\n', iter, T2(iter));
                                    fclose('all');
                                    matrix_output = [Xs(iter,:); yfondo(iter,:)];
                                    self.printmfract(matrix_output);    

                                    if check_terminador2 < terminador
                                        PlatosStrip = iter;
                                        fid = fopen('Output_LM.dat', 'a+');
                                        fprintf(fid, 'El plato %d es el fondo\n', iter);
                                        fprintf(1, 'El plato %d es el fondo\n', iter);
                                        fclose('all');
                                        break
                                    end
                                end
                        end
                    PlatosTot = (PlatosStrip) + PlatosRect; %Obtengo núm de platos totales
                    self.PlatosTot = PlatosTot;
                end
                
                yfondoaux = yfondo;
                T2aux = T2;
                Xsaux = Xs;
                for h = 1:length(T2)
                    T2(end - h + 1) = T2aux(h);
                end
                
                tamano = size(yfondo);
                for h = 1:tamano(1)
                    yfondo(end - h + 1, :) = yfondoaux(h, :);
                end
                
                tamano = size(Xs);
                for h = 1:tamano(1)
                    Xs(end - h + 1, :) = Xsaux(h, :);
                end
                PlatosTot = (PlatosStrip) + PlatosRect + 1; %Obtengo núm de platos totales
            elseif self.startbottom
                if self.lk == 1 && self.hk == self.num_sust
                    for  kui = 1:2
                        try 
                            self.Dconc = X1(end,:);
                            self.Dconc = self.Dconc./(sum(self.Dconc));
                            REPESADO = (self.reco_top.*self.concFeed.*self.flujo)./sum((self.reco_top.*self.concFeed.*self.flujo));
                            REPESADO = REPESADO(self.hk);
                            if REPESADO > self.Dconc(self.hk)
                                RevisionPesado = (self.concFeed.*self.flujo - self.Dflow.*self.Dconc)./(sum(self.concFeed.*self.flujo - self.Dflow.*self.Dconc));
                                RevisionPesado = RevisionPesado(self.hk);
                                self.Bconc(self.hk) = RevisionPesado;
                                self.Bconc = self.Bconc./(sum(self.Bconc));
                            end
                        catch ME
                            self.Dconc = (self.reco_top.*self.concFeed.*self.flujo)./sum((self.reco_top.*self.concFeed.*self.flujo));
                            self.Bconc = (self.reco_bottom.*self.concFeed.*self.flujo)./sum((self.reco_bottom.*self.concFeed.*self.flujo));
                        end
                        ratioLV = (self.reflujo + 1)./(self.reflujo);
                        ratioBV = 1./(self.reflujo);
                        Ls = self.cfeed(1).molF.*self.cfeed(1).q + L;
                        Vs = Ls - self.Bflow;
                        ratioLsVs = Vs./Ls;
                        ratioVsB = self.Bflow ./ Ls;
                        fid = fopen('Output_LM.dat','a+');    
                        fprintf(fid,'Para el plato de fondo, los flujos molares de liquido y vapor son: Ls=%f Vs=%f\n', Ls , Vs );
                        fprintf(1,'Para el plato de fondo, los flujos molares de liquido y vapor son: Ls=%f Vs=%f\n', Ls , Vs );
                        fprintf(fid, 'Con ello, se calcula las relaciones Vs/Ls = %f y W/Vs = %f\n\n', ratioLsVs, ratioVsB);
                        fprintf(1, 'Con ello, se calcula las relaciones Vs/Ls = %f y W/Vs = %f\n\n', ratioLsVs, ratioVsB);
                        X = self.Bconc;
                        xB = self.Bconc;
                        x = self.concFeed;
                        xD = self.Dconc;
                        xB = self.Bconc;
                        terminador = (x(LKey) - ((xD(LKey).*(1 - q))/( 1+reflux)))/(x(HKey)- ((xD(HKey).*(1-q))./(reflux + 1)));
                        check_terminador1 = terminador - 1;            
                        fprintf(fid,'En la alimentacion xF(LKey)/xF(HKey)=%f\n', terminador);
                        fprintf(1,'En la alimentacion xF(LKey)/xF(HKey)=%f\n', terminador);
                        fclose('all');
                        iter = 0;
                        Xs = X;
                        while check_terminador1 <= terminador
                            iter = iter + 1;
                            MezclaPlato_i(iter) = Mezcla(comp, Xs(iter, :), kij);
                            [Tb, sup_y, ~, flag3]= self.MEdE.BubbleT(P, MezclaPlato_i(iter));
                            [T2(iter), yfondo(iter,:), K, flag3]= self.MEdE.BubbleT(P, MezclaPlato_i(iter), sup_y, [], [], Tb);
                            for i = 1:num_sust
                                alfa2(iter, i) = K(i)./K(HKey);
                            end

                            for i = 1:num_sust
                                Xs(iter+1, i) = ratioLsVs.*yfondo(iter,i) + ratioVsB.*xB(i) ;
                            end
                            check_terminador2 = Xs(iter + 1, self.lk)./(Xs(iter + 1, self.hk));
                            fid = fopen('Output_LM.dat', 'a+');
                            fprintf(fid, 'El plato %d subiendo del calderín tiene una composición: T=%f K\n      X       Y\n\n', iter, T2(iter));
                            fprintf(1, 'El plato %d subiendo del calderín tiene una composición: T=%f K\n      X       Y\n\n', iter, T2(iter));
                            fclose('all');
                            matrix_output = [Xs(iter,:); yfondo(iter,:)];
                            self.printmfract(matrix_output);    

                            if check_terminador2 > terminador
                                PlatosStrip = iter;
                                fid = fopen('Output_LM.dat', 'a+');
                                fprintf(fid, 'El plato %d se descarta y la alimentación es fijada en el anterior. Se continúa con la zona de rectificación\n', iter);
                                fprintf(1, 'El plato %d se descarta y la alimentación es fijada en el anterior. Se continúa con la zona de rectificación\n', iter);
                                fclose('all');
                                break
                            end
                        end
                            %%%%%%%% Se continúa operando en la región de
                            %%%%%%%% enriquecimiento

                            X1 = Xs(end, :);
                            xD = self.Dconc;
                            x = self.concFeed;
                            terminador = self.flujo.*(1 - self.reco_hk).*self.concFeed(self.hk);
                            check_terminador1 = terminador - 1;        
                            fid = fopen('Output_LM.dat', 'a+');
                            fprintf(fid,'Pasando el plato de alimentación con la fracción de liviano a pesado de %f \n', terminador);
                            fprintf(1,'Pasando el plato de alimentación con la fracción de liviano a pesado de %f \n', terminador);
                            fclose('all');
                            iter = 0;
                            while check_terminador1 <= terminador
                                iter = iter + 1;
                                MezclaPlato_i(iter) = Mezcla(comp, X1(iter, :), kij);
                                [Tb, sup_y, ~, flag3]= self.MEdE.BubbleT(P, MezclaPlato_i(iter));
                                [T1(iter), Y(iter,:), K, flag3]= self.MEdE.BubbleT(P, MezclaPlato_i(iter), sup_y, [], [], Tb);
                                for i = 1:num_sust
                                    alfa2(iter, i) = K(i)./K(HKey);
                                end

                                for i = 1:num_sust
                                    X1(iter+1, i) = ratioLV.*Y(iter,i) - ratioBV.*xD(i) ;
                                end
                                if any(X1(iter+1,:) < 0)
                                    indice = find(X1(iter+1,:) < 0);
                                    for indic = min(indice):max(indice)
                                        X1(iter+1, indic) = 1e-15;
                                    end
                                    X1(iter+1, :) = X1(iter+1, :)./sum(X1(iter+1, :));
                                end
    %                             iter = iter + 1;
    %                             for i = 1:num_sust
    %                                 yfondo(iter, i) = ratioLsVs.*Xs(iter,i) + ratioVsB.*xB(i) ;
    %                             end
    %                             MezclaPlato_i(iter+1) = Mezcla(comp, yfondo(iter, :), kij);
    %                             [Tr, sup_x, ~, flag3]= self.MEdE.DewT(P, MezclaPlato_i(iter+1));
    %                             [T2(iter), Xs(iter+1,:), K, flag3]= self.MEdE.DewT(P, MezclaPlato_i(iter+1), sup_x, [], [], Tr);
    %                             for i = 1:num_sust
    %                                 alfa2(iter+1, i) = K(i)./K(HKey);
    %                             end
                                check_terminador2 = self.Dflow.*X1(iter, self.hk);
                                fid = fopen('Output_LM.dat', 'a+');
                                fprintf(fid, 'El plato %d subiendo del calderín tiene una composición: T=%f K\n      X       Y\n\n', iter, T1(iter));
                                fprintf(1, 'El plato %d subiendo del calderín tiene una composición: T=%f K\n      X       Y\n\n', iter, T1(iter));
                                fclose('all');
                                matrix_output = [X1(iter,:); Y(iter,:)];
                                self.printmfract(matrix_output);    

                                if check_terminador2 < terminador
                                    PlatosRect = iter;
                                    fid = fopen('Output_LM.dat', 'a+');
                                    fprintf(fid, 'El plato %d se descarta y la alimentación es fijada en el anterior. Se continúa con la zona de rectificación\n', iter);
                                    fprintf(1, 'El plato %d se descarta y la alimentación es fijada en el anterior. Se continúa con la zona de rectificación\n', iter);
                                    fclose('all');
                                    break
                                end
                            end
                    end
                end
                tamano = size(Xs);
                Xsaux = zeros(tamano(1)+1, tamano(2));
                Xsaux(1,:) = self.Bconc;
                Xsaux(2:end, :) = Xs(1:end, :);
                Xs = Xsaux;
                tamano = size(Y);
                Yaux = zeros(tamano(1)+1, tamano(2));
                Yaux(1,:) = self.Dconc;
                Yaux(2:end, :) = Y(1:end, :);
                Y = Yaux;
                PlatosTot = (PlatosStrip) + PlatosRect; %Obtengo núm de platos totales
                
                Yaux = Y;
                T1aux = T1;
                X1aux = X1;
                for h = 1:length(T1)
                    T1(end - h + 1) = T1aux(h);
                end
                
                tamano = size(Y);
                for h = 1:tamano(1)
                    Y(end - h + 1, :) = Yaux(h, :);
                end
                
                tamano = size(X1);
                for h = 1:tamano(1)
                    X1(end - h + 1, :) = X1aux(h, :);
                end
            end
            %Recalculo composiciones cercanas a la alimentación
            
            no_seguir=0; %Si ya se ha corregido suficiente las composiciones no sigue
            fid = fopen('Output_LM.dat','a+'); 
            fprintf(fid,'Corrigiendo las composiciones para la alimentación Plato %d\n     x       y    \n', PlatosRect);
            fprintf(1,'Corrigiendo las composiciones para la alimentación Plato %d\n     x       y    \n', PlatosRect);
            fclose('all');
            
            tamano = size(X1)-1;
            for i = 1: tamano(1) + 1
                for j = 1:num_sust
                    if X1(i,j) == 0
                        X1(i,j) = 1e-15;
                    end
                end
            end
            tamano2 = size(Xs) - 1;
            for i = 1: tamano2(1) + 1
                for j = 1:num_sust
                    if Xs(i,j) == 0
                        Xs(i,j) = 1e-15;
                    end
                end
            end
            if num_sust > HKey
            %Fuerzo la composición de pesado de alimentación a ser la de zona Striping
            plato = 0;
            while (plato < length(T1)) && (no_seguir == 0)
                plato = plato + 1;
                for iter_pesados=HKey+1:num_sust
                    save = X1(end-plato+1,iter_pesados);
                        if plato == 1
                            X1(end-plato+1,iter_pesados) = Xs(end,iter_pesados);
                        elseif plato > 1
                            X1(end-plato+1,iter_pesados) = X1(end-plato+2, iter_pesados).*alfa1(end-plato+2, iter_pesados).*X1(end - plato +1, HKey)./(X1(end-plato+2, HKey));
                        end
                        if plato == 1
                            comp_pesados = (Xs(end,iter_pesados)-save);
                        elseif plato > 1
                            comp_pesados = X1(end -plato+1,iter_pesados)-save;
                        end
                        suma = sum(X1(end-plato+1,1:iter_pesados-1));
                        if (abs(X1(end-plato+1,iter_pesados)-save)<1e-5) && iter_pesados == HKey + 1
                            no_seguir = 1;
                        end
                        for iter_livianos = 1:HKey
                            X1(end-plato+1,iter_livianos) = X1(end-plato+1,iter_livianos)*((suma-comp_pesados)/suma);
                        end
                        fid = fopen('Output_LM.dat','a+');
                        if PlatosRect - plato + 1 >0
                            fprintf(fid,'Corrigiendo la composición del plato %d \n', PlatosRect - plato+1);
                            fprintf(1,'Corrigiendo la composición del plato %d \n', PlatosRect - plato+1);
                        else 
                            fprintf(fid,'Corrigiendo la composición del condensador %d \n', PlatosRect - plato+1);
                             fprintf(1,'Corrigiendo la composición del condensador %d \n', PlatosRect - plato+1);
                        end
                        MezclaPlato_i(end-plato+1).conc = X1(end - plato+1,:)./sum(X1(end - plato+1,:));
                        [Tb, sup_y, ~, flag5]= self.MEdE.BubbleT(P, MezclaPlato_i(end-plato+1));
                        [Treval(tamano(1)-plato+2), Y(end-plato+1,:), K, flag5]= self.MEdE.BubbleT(P, MezclaPlato_i(end-plato +1), sup_y, [], [], Tb);
                        matrix_output = [X1(end-plato+1,:);Y(end-plato+1,:)];
                        self.printmfract(matrix_output);                
                        fid = fopen('Output_LM.dat','a+'); 
                        fprintf(fid,'Y la temperatura del plato es %f\n\n',Treval(tamano(1)-plato+2));
                        fprintf(1,'Y la temperatura del plato es %f\n\n',Treval(tamano(1)-plato+2));
                        fclose('all');
                end
                
            end
                TRectreval = T1;
                Treval = Treval(Treval~=0);
                presionreval = zeros(length(T1),1);
                for i = length(Treval)-1:-1:0
                    TRectreval(end-i+1) = Treval(end-i);
                end
                TRectreval = TRectreval(2:end);
            end
            %%%%%%%%%%%%%%%%%%%%%%
            if LKey > 1
            plato = 0;
            while (plato < length(T2)) && (no_seguir == 0)
                plato = plato + 1;
                for iter_liviano=1:LKey-1
                    if plato == 1
                        Xs(end,:) = X1(end,:);
                    elseif plato > 1
                        save = Xs(end-plato+1,iter_liviano);
                        Xs(end-plato+1,iter_liviano) = Xs(end-plato+2, iter_liviano).*alfa2(end-plato+2, LKey).*Xs(end - plato +1, LKey)./(alfa2(end-plato+2, iter_liviano)*Xs(end-plato+2, LKey));
                        if plato > 1
                            comp_livianos = Xs(end-plato+1,iter_liviano)-save;
                        end
                        suma = sum(Xs(end-plato+1,iter_liviano:end));
                        if abs(Xs(end-plato+1,iter_liviano)-save)<1e-5 && (iter_liviano == LKey - 1)
                            no_seguir = 1;
                        end
                        for iter_pesado = iter_liviano+1:num_sust
                            Xs(end-plato+1,iter_pesado) = Xs(end-plato+1,iter_pesado)*((suma-comp_livianos)/suma);
                        end
                        fid = fopen('Output_LM.dat','a+');
                        if PlatosStrip - plato + 1 >0
                            fprintf(fid,'Corrigiendo la composición del plato %d \n', PlatosStrip - plato+1);
                            fprintf(1,'Corrigiendo la composición del plato %d \n', PlatosStrip - plato+1);
                        else 
                            fprintf(fid,'Corrigiendo la composición del rehervidor %d \n', PlatosStrip - plato+1);
                            fprintf(1,'Corrigiendo la composición del rehervidor %d \n', PlatosStrip - plato+1);
                        end
                        LiquidoPlato_i(end-plato+1).conc = Xs(end - plato+1,:)./sum(Xs(end - plato+1,:));
                        [Tr, sup_y, ~, flag5]= self.MEdE.BubbleT(P, LiquidoPlato_i(end-plato+1));
                        [T2reval(tamano2(1)-plato+2), yfondo(end - plato +1,:), K, flag5]= self.MEdE.BubbleT(P, LiquidoPlato_i(end-plato +1), sup_y, [], [], Tr);
                        matrix_output = [Xs(end-plato+1,:);yfondo(end-plato+1,:)];
                        self.printmfract(matrix_output);                
                        fid = fopen('Output_LM.dat','a+'); 
                        fprintf(fid,'Y la temperatura del plato es %f\n\n',T2reval(tamano2(1)-plato+2));
                        fprintf(1,'Y la temperatura del plato es %f\n\n',T2reval(tamano2(1)-plato+2));
                        fclose('all');
                    end
                end
                
            end
                try
                    TStripreval = T2;
                    T2reval = T2reval(T2reval~=0);
                    presionreval = zeros(length(T1),1);
                    
                    for i = length(T2reval)-1:-1:0
                        TStripreval(end-i) = T2reval(end-i);
                    end
                    TStripreval = TStripreval(2:end);
                    T2Stripreval = TStripreval;
                    for i = 1:length(TStripreval)
                        TStripreval(i) = T2Stripreval(end-i+1);
                    end
                    Xsaux = Xs;
                    yfondoaux = yfondo;
                    
                    for i =1:length(Xsaux)
                        Xs(i,:) = Xsaux(end-i+1,:);
                    end
                    for i = 1:length(yfondoaux)
                        yfondo(i,:) = yfondoaux(end-i + 1, :);
                    end
                catch
                    TStripreval = T2;
                    T2Stripreval = TStripreval;
                    for i = 1:length(TStripreval)
                        TStripreval(i) = T2Stripreval(end-i+1);
                    end
                    Xsaux = Xs;
                    yfondoaux = yfondo;
                    
                    for i =1:length(Xsaux)
                        Xs(i,:) = Xsaux(end-i+1,:);
                    end
                    for i = 1:length(yfondoaux)
                        yfondo(i,:) = yfondoaux(end-i + 1, :);
                    end
                end
            end
            self.TPlatoPlato = zeros(1, PlatosTot+1);
            try
                for i = 1:length(TRectreval)
                    self.TPlatoPlato(i) = TRectreval(i);
                end
            catch
                for i = 1:length(T1)
                    self.TPlatoPlato(i) = T1(i);
                end
            end
            try
                for i = length(TRectreval)+1:1:length(TRectreval)+length(TStripreval)
                    self.TPlatoPlato(i) = TStripreval(end - i + length(TRectreval)+1);
                end
            catch
                try
                    for i = length(T1)+1:1:length(T1)+length(TStripreval)
                        self.TPlatoPlato(i) = T2(end - i + length(T1)+1);
                    end
                catch
                    for i = length(T1)+1:1:length(T1)+length(T2)
                        self.TPlatoPlato(i) = T2(end - i + length(T1)+1);
                    end
                end
            end
            self.concLiqPlato = zeros(PlatosTot+1, num_sust);
            self.concVapPlato = zeros(PlatosTot+1, num_sust);
            try
                for i = 1:length(TRectreval)
                    self.concLiqPlato(i,:) = X1(i, :)./sum(X1(i, :));
                    self.concVapPlato(i,:) = Y(i,:)./sum(Y(i,:));
                end
            catch
                for i = 1:length(T1)
                    self.concLiqPlato(i,:) = X1(i, :)./sum(X1(i, :));
                    self.concVapPlato(i,:) = Y(i,:)./sum(Y(i,:));
                end
            end
            try
                for i = length(TRectreval)+1:1:length(TRectreval)+length(TStripreval)-1
                    self.concLiqPlato(i,:) = Xs(end - i + length(TRectreval),:)./sum(Xs(end - i + length(TRectreval),:));
                    self.concVapPlato(i,:) = yfondo(end -i + length(TRectreval),:)./sum(yfondo(end -i + length(TRectreval),:));                   
                end
            catch
                try
                    for i = length(T1)+1:1:length(T1)+length(TStripreval)-1
                        self.concLiqPlato(i,:) = Xs(end - i + length(T1),:)./sum(Xs(end - i + length(T1),:));
                        self.concVapPlato(i,:) = yfondo(end -i + length(T1),:)./sum(yfondo(end -i + length(T1),:));

                    end
                catch
                    for i = length(T1)+1:1:length(T1)+length(T2)-1
                        self.concLiqPlato(i,:) = Xs(end - i + length(T1),:)./sum(Xs(end - i + length(T1),:));
                        self.concVapPlato(i,:) = yfondo(end -i + length(T1),:)./sum(yfondo(end -i + length(T1),:));

                    end
                end
            end
                self.PlatosTot = PlatosTot;
            self.etapaoptima = PlatosRect;
            indice = 0;
            for variableforauxiliar = 1:length(self.TPlatoPlato) 
                if self.TPlatoPlato(end + 1 - variableforauxiliar) > 1e-8
                    if indice > 0
                        self.TPlatoPlato = self.TPlatoPlato(1:length(self.TPlatoPlato)-indice);
                    end
                    break
                else
                    indice = indice + 1;
                end
            end
            indice = 0;
            tamano = size(self.concVapPlato);
            for variableforauxiliar = 1:tamano(1)
                if sum(self.concVapPlato(end + 1 - variableforauxiliar,:)) > 1e-8
                    if indice > 0
                        self.concVapPlato = self.concVapPlato(1:tamano(1)-indice, :);
                    end
                    break
                else
                    indice = indice + 1;
                end
            end
            indice = 0;
            tamano = size(self.concLiqPlato);
            for variableforauxiliar = 1:tamano(1)
                if sum(self.concLiqPlato(end + 1 - variableforauxiliar,:)) > 1e-8
                    if indice > 0
                        self.concLiqPlato = self.concLiqPlato(1:tamano(1)-indice, :);
                    end
                    break
                else
                    indice = indice + 1;
                end
            end
            D = self.Dflow;
            B = self.Bflow;
            MezclaY = Mezcla(comp, self.Dconc, kij);
            Vapor1 = Corriente(MezclaY, self.TPlatoPlato(2),'T', 1, 'P', self.Dflow.*(reflux+1), 'm', self.MEdE); 
            VaporD =  Corriente(MezclaY, P,'P', 1, 'x', self.Dflow, 'm', self.MEdE); 
            self.TPlatoPlato(2) = Vapor1.T;
            LiquidoD = Corriente(MezclaY, 0, 'x', self.TPlatoPlato(1), 'T', self.Dflow, 'm', self.MEdE);
            Liquido0 = LiquidoD;
            Liquido0.molF = self.Dflow.*reflux;
            Liquido0.entalpia();
            Qc = D.*((reflux+1).*(Vapor1.H) - reflux.*Liquido0.H - LiquidoD.H);
            MezclaB = Mezcla(comp, self.Bconc, kij);
            LiquidoW = Corriente(MezclaB, 0, 'x', P, 'P', self.Bflow, 'm', self.MEdE);
            F = self.flujo;
            Qb = D.*(VaporD.H) + B.*LiquidoW.H + Qc + self.heat_loss - F.*self.cfeed.H;
            self.Qc = Qc;   
            self.Qb = -Qb;
            Lim1 = D*reflux;
            Gi = D*(reflux+1);
            self.ratioLV = reflux;
            self.ratioLV(2) = (D.*(reflux + 1) - D)./Gi;
            for i = 1:PlatosRect
                Liquido_i = Mezcla(comp, self.concLiqPlato(i, :), kij);
                Vapor_i = Mezcla(comp, self.concVapPlato(i+1,:), kij);
                Liquidoi(i) = Corriente(Liquido_i, 0, 'x', self.TPlatoPlato(i+1), 'T', Lim1(i), 'm', self.MEdE);
                Vapori(i) = Corriente(Vapor_i, 1, 'x', P, 'P', Lim1(i), 'm', self.MEdE);
                Gi(i+1) = ((Qc + D.*(LiquidoD.H - Liquidoi(i).H))./(Vapori(i).H - Liquidoi(i).H) );
                if i == PlatosRect
                    break
                end
                Lim1(i+1) = Gi(i+1) - D;
                self.ratioLV(i+2) = Lim1(i+1)./(Gi(i+1));
            end
            self.ratioLsVs = zeros(1, PlatosStrip-1);
            self.ratioLsVs(1) = (self.flujo.*self.q + Lim1(end))./(Gi(end) + self.flujo.*(1-self.q));
            for i = 0:PlatosStrip-3
                Liquido_is = Mezcla(comp, self.concLiqPlato(end-i-1,:), kij);
                Vapor_is = Mezcla(comp, self.concVapPlato(end-i,:), kij);
                Liquidois(i+1) = Corriente(Liquido_is, 0, 'x', self.TPlatoPlato(end - i - 1), 'T', B, 'm', self.MEdE);
                Vaporis(i+1) = Corriente(Vapor_is, 1, 'x', P, 'P', B, 'm', self.MEdE);
                Gis(i+1) = (Qb + B.*(Liquidois(i+1).H - LiquidoW.H))/(Vaporis(i+1).H - Liquidois(i+1).H);
                Lis(i+1) = Gis(i+1) + B;
                self.ratioLsVs(end - i) = Lis(i+1)./Gis(i+1);
            end
            salida = self;
%             
%             no_seguir=0; %Si ya se ha corregido suficiente las composiciones no sigue
%             fid = fopen('Output_LM.dat','a+'); 
%             fprintf(fid,'Corrigiendo las composiciones para la alimentación Plato %d\n     x       y    \n', PlatosRect);
%             fclose('all');
%             if num_sust > HKey
%             Fuerzo la composición de pesado de alimentación a ser la de zona Striping
%             plato = 0;
%             tamano = size(X1)-1;
%             for i = 1: tamano(1) + 1
%                 for j = 1:num_sust
%                     if X1(i,j) == 0
%                         X1(i,j) = 1e-15;
%                     end
%                 end
%             end
%             tamano2 = size(Xs) - 1;
%             for i = 1: tamano2(1) + 1
%                 for j = 1:num_sust
%                     if Xs(i,j) == 0
%                         Xs(i,j) = 1e-15;
%                     end
%                 end
%             end
%             while (plato < length(T1)) && (no_seguir == 0)
%                 plato = plato + 1;
%                 for iter_pesados=HKey+1:num_sust
%                     save = X1(end-plato+1,iter_pesados);
%                         if plato == 1
%                             X1(end-plato+1,iter_pesados) = Xs(end,iter_pesados);
%                         elseif plato > 1
%                             X1(end-plato+1,iter_pesados) = X1(end-plato+2, iter_pesados).*alfa1(end-plato+2, iter_pesados).*X1(end - plato +1, HKey)./(X1(end-plato+2, HKey));
%                         end
%                         if plato == 1
%                             comp_pesados = (Xs(end,iter_pesados)-save);
%                         elseif plato > 1
%                             comp_pesados = X1(end -plato+1,iter_pesados)-save;
%                         end
%                         suma = sum(X1(end-plato+1,1:iter_pesados-1));
%                         if (abs(X1(end-plato+1,iter_pesados)-save)<1e-5) && iter_pesados == HKey + 1
%                             no_seguir = 1;
%                         end
%                         for iter_livianos = 1:HKey
%                             X1(end-plato+1,iter_livianos) = X1(end-plato+1,iter_livianos)*((suma-comp_pesados)/suma);
%                         end
%                         fid = fopen('Output_LM.dat','a+');
%                         if PlatosRect - plato + 1 >0
%                             fprintf(fid,'Corrigiendo la composición del plato %d \n', PlatosRect - plato+1);
%                         else 
%                             fprintf(fid,'Corrigiendo la composición del condensador %d \n', PlatosRect - plato+1);
%                         end
%                         MezclaPlato_i(end-plato+1).conc = X1(end - plato+1,:);
%                         [Tb, sup_y, ~, flag5]= self.MEdE.BubbleT(P, MezclaPlato_i(end-plato+1));
%                         [Treval(tamano(1)-plato+2), Y(end-plato+1,:), K, flag5]= self.MEdE.BubbleT(P, MezclaPlato_i(end-plato +1), sup_y, [], [], Tb);
%                         matrix_output = [X1(end-plato+1,:);Y(end-plato+1,:)];
%                         self.printmfract(matrix_output);                
%                         fid = fopen('Output_LM.dat','a+'); 
%                         fprintf(fid,'Y la temperatura del plato es %f\n\n',Treval(tamano(1)-plato+2));
%                         fclose('all');
%                 end
%                 
%             end
%                 TRectreval = T1;
%                 Treval = Treval(Treval~=0);
%                 presionreval = zeros(length(T1),1);
%                 for i = length(Treval)-1:-1:0
%                     TRectreval(end-i+1) = Treval(end-i);
%                 end
%                 TRectreval = TRectreval(2:end);
%             end
%             %%%%%%%%%%%%%%%%%%%%%
%             if LKey > 1
%             plato = 0;
%             while (plato < length(T2)) && (no_seguir == 0)
%                 plato = plato + 1;
%                 for iter_liviano=1:LKey-1
%                     if plato == 1
%                         Xs(end,:) = X1(end,:);
%                     elseif plato > 1
%                         save = Xs(end-plato+1,iter_liviano);
%                         Xs(end-plato+1,iter_liviano) = Xs(end-plato+2, iter_liviano).*alfa2(end-plato+2, LKey).*Xs(end - plato +1, LKey)./(alfa2(end-plato+2, iter_liviano)*Xs(end-plato+2, LKey));
%                         if plato > 1
%                             comp_livianos = Xs(end-plato+1,iter_liviano)-save;
%                         end
%                         suma = sum(Xs(end-plato+1,iter_liviano:end));
%                         if abs(Xs(end-plato+1,iter_liviano)-save)<1e-5 && (iter_liviano == LKey - 1)
%                             no_seguir = 1;
%                         end
%                         for iter_pesado = iter_liviano+1:num_sust
%                             Xs(end-plato+1,iter_pesado) = Xs(end-plato+1,iter_pesado)*((suma-comp_livianos)/suma);
%                         end
%                         fid = fopen('Output_LM.dat','a+');
%                         if PlatosStrip - plato + 1 >0
%                             fprintf(fid,'Corrigiendo la composición del plato %d \n', PlatosStrip - plato+1);
%                         else 
%                             fprintf(fid,'Corrigiendo la composición del rehervidor %d \n', PlatosStrip - plato+1);
%                         end
%                         LiquidoPlato_i(end-plato+1).conc = Xs(end - plato+1,:);
%                         [Tr, sup_y, ~, flag5]= self.MEdE.BubbleT(P, LiquidoPlato_i(end-plato+1));
%                         [T2reval(tamano2(1)-plato+2), yfondo(end - plato +1,:), K, flag5]= self.MEdE.BubbleT(P, LiquidoPlato_i(end-plato +1), sup_y, [], [], Tr);
%                         matrix_output = [Xs(end-plato+1,:);yfondo(end-plato+1,:)];
%                         self.printmfract(matrix_output);                
%                         fid = fopen('Output_LM.dat','a+'); 
%                         fprintf(fid,'Y la temperatura del plato es %f\n\n',T2reval(tamano2(1)-plato+2));
%                         fclose('all');
%                     end
%                 end
%                 
%             end
%                 TStripreval = T2;
%                 T2reval = T2reval(T2reval~=0);
%                 presionreval = zeros(length(T1),1);
%                 for i = length(T2reval)-1:-1:0
%                     TStripreval(end-i) = T2reval(end-i);
%                 end
%             end
%             self.TPlatoPlato = zeros(1, PlatosTot+1);
%             
%             for i = 1:length(TRectreval)
%                 self.TPlatoPlato(i) = TRectreval(i);
%             end
%             for i = length(TRectreval)+1:1:length(TRectreval)+length(TStripreval)
%                 self.TPlatoPlato(i) = TStripreval(end - i + length(TRectreval)+1);
%             end
%             self.concLiqPlato = zeros(PlatosTot+1, num_sust);
%             self.concVapPlato = zeros(PlatosTot+1, num_sust);
%             for i = 1:length(TRectreval)
%                 self.concLiqPlato(i,:) = X1(i, :);
%                 self.concVapPlato(i,:) = Y(i,:);
%             end
%             for i = length(TRectreval)+1:1:length(TRectreval)+length(TStripreval)-1
%                 self.concLiqPlato(i,:) = Xs(end - i + length(TRectreval),:);
%                 self.concVapPlato(i,:) = yfondo(end -i + length(TRectreval),:);
%             end
%             self.PlatosTot = PlatosTot+1;
%             self.etapaoptima = PlatosRect+1;
%             Liquido1 = Mezcla(comp, self.concLiqPlato(1,:), kij);
%             [Tb , self.Dconc, ~, flag] = self.MEdE.BubbleT(P, Liquido1);
%             [~, self.Dconc, ~, flag] = self.MEdE.BubbleT(P, Liquido1, self.Dconc, [], [], Tb);
%             self.Dflow = self.reco_lk.*self.flujo.*self.concFeed(self.lk)./(self.Dconc(self.lk));
%             self.Bflow = 1-self.Dflow;
%             
%             D = self.Dflow;
%             B = self.Bflow;
%             MezclaY = Mezcla(comp, self.Dconc, kij);
%             Vapor1 = Corriente(MezclaY, self.TPlatoPlato(1),'T', P, 'P', self.Dflow.*(reflux+1), 'm', self.MEdE); 
%             if abs(Vapor1.beta - 1) <5e-4
%                 Vapor1 = Corriente(MezclaY, 1,'x', P, 'P', self.Dflow.*(reflux+1), 'm', self.MEdE); 
%             end
%             LiquidoD = Corriente(MezclaY, 0, 'x', P, 'P', self.Dflow, 'm', self.MEdE);
%             Liquido0 = LiquidoD;
%             Liquido0.molF = self.Dflow.*reflux;
%             Liquido0.entalpia();
%             Qc = D.*((reflux+1).*(Vapor1.H) - reflux.*Liquido0.H - LiquidoD.H);
%             MezclaB = Mezcla(comp, self.Bconc, kij);
%             LiquidoW = Corriente(MezclaB, 0, 'x', P, 'P', self.Bflow, 'm', self.MEdE);
%             F = self.flujo;
%             Qb = D.*(LiquidoD.H) + B.*LiquidoW.H + Qc + self.heat_loss - F.*self.cfeed.H;
%             self.Qc = Qc;
%             self.Qb = Qb;
%             Lim1 = D*reflux;
%             Gi = D*(reflux+1);
%             self.ratioLV = Lim1./Gi;
%             for i = 1:PlatosRect
%                 Liquido_i = Mezcla(comp, self.concVapPlato(i+1, :), kij);
%                 Vapor_i = Mezcla(comp, self.concLiqPlato(i,:), kij);
%                 Liquidoi(i) = Corriente(Liquido_i, 0, 'x', P, 'P', Gi(i), 'm', self.MEdE);
%                 Vapori(i) = Corriente(Vapor_i, 1, 'x', P, 'P', Lim1(i), 'm', self.MEdE);
%                 Gi(i+1) = ((Qc + D.*(LiquidoD.H - Liquidoi(i).H))./(Vapori(i).H - Liquidoi(i).H) );
%                 Lim1(i+1) = Gi(i+1) - D;
%                 self.ratioLV(i+1) = Lim1(i+1)./(Gi(i+1));
%             end
%             self.ratioLsVs = zeros(1, PlatosStrip-1);
%             for i = 0:PlatosStrip-2
%                 Liquido_is = Mezcla(comp, self.concVapPlato(end-i,:), kij);
%                 Vapor_is = Mezcla(comp, self.concLiqPlato(end-i,:), kij);
%                 Liquidois(i+1) = Corriente(Liquido_is, 0, 'x', P, 'P', B, 'm', self.MEdE);
%                 Vaporis(i+1) = Corriente(Vapor_is, 1, 'x', P, 'P', B, 'm', self.MEdE);
%                 Gis(i+1) = (Qb + B.*(Liquidois(i+1).H - LiquidoW.H))/(Vaporis(i+1).H - Liquidois(i+1).H);
%                 Lis(i+1) = Gis(i+1) + B;
%                 self.ratioLsVs(end - i) = Lis(i+1)./Gis(i+1);
%             end
%             salida = self;
        end
        function salida = reset(self) %resetea cada propiedad conservando nombre
            
            self.platos_reales = false;
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
            self.reco_bottom = false;
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
            self.gdlibertad = 1;
            self.volat = false;
            self.volat_top = false;
            self.volat_bot = false;
            self.concFeed = false;
            self.concliqtope = false;
            self.concvapfondo = false;
            self.beta_feed = false;
            self.heat_loss = 0;
            self.ratioLV = false;
            self.ratioLsVs = false;
            self.ratioLG = false;
            self.volat_top = false;
            self.volat_bot = false;
            self.volat = false;
            salida = self;
        end
    end    
end

