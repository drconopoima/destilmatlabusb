classdef Hengstebeck < handle
    %Hengstebeck Resolución del método de Hengstebeck para sistemas
    %pseudo-binarios con aproximación del flujo de componentes no claves
    %por el método propuesto por Jenny (1939)
    %   Luis Jesús Díaz Manzo
    
    properties
        id
        etapas
        perfil_p
        perfil_t
        perfil_t_bin
        qc
        qb
        dflujo
        bflujo
        reflujo
        num_sust
        alimentaciones
        cdest
        cfond
        salidas
        entradas
        gdlibertad
        comp
        cmpbinario
        sust
        MEdE
        dsubi
        bsubi
        concfeed
        lk
        hk
        Fi
        Di
        Bi
        xDi
        xBi
        UWi
        reco_lk
        reco_hk
        volat_top
        volat_bot
        volat
        E_x
        E_y
        Tbin
        Rmin
        Nmin
        metodo = 1
        Murphree
        Mur_x
        Mur_y
        klk = 0;
        rectaq
        TPlato
        TPlatobin
        intersectx
        intersecty
        intersectequilx 
        intersectequil
        rectarectmin
        rectastripmin
        vivo
        platosmin
        platosaefficmin
        etapasemin
        condparcial = false
        rectaoperacion
        rectaoper_min
        vap_min
        rectasminima
        rectas
    end
    
    methods
        function self = Hengstebeck(alim, lightk, recup_top, heavyk, recup_bot, salid, reflux, dest, bot, Pperfil, effic, vivo, parcial, nombre)
            self.comp = cell.empty(0,0);
            self.sust = Sustancia.empty(0,0);
            if nargin > 0 && ~isempty(alim) && isa(alim, 'Corriente') && isa(alim(1), 'Corriente')
                long = length(alim);
                self.entradas = zeros(1, long*3);        
                for i = 1:1:long
                    corr = alim(i);
                    susta = corr.comp;
                    for j = 1:length(susta)
                        if ~any(strcmp(self.comp, susta(j).id))
                            self.comp{length(self.comp) + 1} = susta(j).id;
                            self.sust(length(self.sust) + 1) = susta(j);
                        end
                        try    
                            if lightk == j && strcmpi(susta(j).id, self.sust(lightk).id)
                                self.entradas(3.*i - 2) = corr.conc(j)/(corr.conc(j)+corr.conc(j+1));
                                if isempty(self.cmpbinario)
                                    self.cmpbinario = Sustancia.empty(0,2);
                                    self.cmpbinario(1) = susta(j);
                                end
                            elseif strcmp(lightk, susta(j).id)
                                self.entradas(3.*i - 2) = corr.conc(j)/(corr.conc(j)+corr.conc(j+1));
                                if isempty(self.cmpbinario)
                                    self.cmpbinario = Sustancia.empty(0,2);
                                    self.cmpbinario(1) = susta(j);
                                end
                            end
                        catch
                        end
                        try
                            if heavyk == j
                                if length(self.cmpbinario) == 1
                                    self.cmpbinario(2) = susta(j);
                                    klkhk = strcat('k', int2str(lightk), int2str(heavyk));
                                    self.klk = alim(1).mezcla.kij(find(strcmpi(alim(1).mezcla.kbinario, klkhk),1,'first'));
                                end
                            elseif strcmp(heavyk, susta(j).id)
                                if length(self.cmpbinario) == 1
                                    self.cmpbinario(2) = susta(j);
                                    klkhk = strcat('k', int2str(lightk), int2str(heavyk));
                                    self.klk = alim(1).mezcla.kij(find(strcmpi(alim(1).mezcla.kbinario, klkhk),1,'first'));
                                end
                            end
                        catch
                            continue
                        end
                    end
                    self.entradas(3*(i-1) + 2) = alim(i).q;
                    self.entradas(3*(i-1) + 3) = alim(i).molF;
                end
                self.num_sust = length(self.comp);
                self.concfeed = zeros(long, length(self.comp));
                self.dsubi = zeros(1, length(self.comp));
                self.bsubi = zeros(1, length(self.comp));
                for i = 1:1:long
                    corr = alim(i);
                    susta = corr.comp;
                    for j = 1:length(susta)
                        indice = find(strcmp(susta(j).id, self.comp),1, 'first');
                        self.concfeed(i, indice) = corr.conc(j);
                    end
                end
                self.Fi = zeros(self.num_sust, 1);
                for l = 1:length(self.comp)
                    for u = 3:3:length(self.entradas)
                        self.Fi(l) = self.Fi(l) + self.entradas(u).*self.concfeed(u/3, l);
                    end
                end
                self.alimentaciones = alim;
                self.MEdE = alim(1).MEdE;
                self.rectaq = zeros(2, (length(self.entradas)/3)*2);
                iter = 0;
                for i = 2:3:length(self.entradas)
                    iter = iter + 1;
                    self.rectaq(:, (iter)*2-1) = [self.entradas(i - 1), self.entradas(i - 1)];
                    if abs(self.entradas(i) - 1) < 1e-5
                        self.rectaq(:, (iter)*2) = [self.entradas(i - 1), 2];
                    elseif abs(self.entradas(i)) < 1e-5
                        self.rectaq(:, (iter)*2) = [-2, self.entradas(i - 1)];
                    else
                        self.rectaq(:, (iter)*2) = [0, -(self.entradas(i - 1)/(self.entradas(i)-1))];
                    end
                end
            end
            
            if nargin > 1 && ~isempty(lightk) 
                self.lk = lightk;
            end
            if nargin > 2   
                if ~isempty(recup_top) && ~isempty(lightk)
                    self.set_reco_top(lightk, recup_top);
                end
            end
            if nargin > 3 && ~isempty(heavyk) 
                self.hk = heavyk;
            end
            if nargin >  4
                if ~isempty(recup_bot) && ~isempty(heavyk)
                    self.set_reco_bot(heavyk, recup_bot);
                end
            end
            if nargin > 5 && ~isempty(salid) && isa(salid, 'double')
                self.salidas = salid;
            end
            if nargin > 6 && ~isempty(reflux)
                self.reflujo = reflux;
            end
            if nargin > 7 && isa(dest, 'Corriente')
                self.cdest = dest;
            else
                self.cdest = Corriente.empty(1,0);
            end
            if nargin > 8 && isa(bot, 'Corriente')
                self.cfond = bot;
            else
                self.cfond = Corriente.empty(1,0);
            end
            if nargin > 9 && isempty(Pperfil)
                self.perfil_p = ones(1, 2)*alim(1).P;
            elseif nargin > 9 && length(Pperfil) == 2
                self.perfil_p = Pperfil;
            elseif nargin > 9 && length(Pperfil) == 1
                self.perfil_p = ones(1, 2)*Pperfil;
            elseif nargin > 9 
                self.perfil_p(1) = Pperfil(1);
                self.perfil_p(2) = Pperfil(end);
            elseif nargin > 0 && isa(alim(1), 'Corriente')
                self.perfil_p = zeros(1,2);
                self.perfil_p(:) = alim(1).P;
            end 
            if nargin > 10 && ~isempty(effic) && effic < 1 && effic > 0
                 self.Murphree = effic;
            else
                self.Murphree = 1;
            end
            if nargin > 11 && ~isempty(vivo) && vivo == 1
                self.vivo = true;
            else 
                self.vivo = false;
            end
            if nargin > 12 && ~isempty(parcial) && parcial
                self.condparcial = 1;
            end
            if nargin > 13 && ~isempty(nombre) && isa(nombre, 'char')
                self.id = nombre;
            end
            if ~isempty(self.reco_lk) && ~isempty(self.lk)
                self.Di(self.lk) = self.reco_lk*self.Fi(self.lk);
                if isempty(self.salidas)
                    self.Bi(self.lk) = (1 - self.reco_lk)*self.Fi(self.lk);
                else
                    for k = 1:3:length(self.salidas)
                        self.Bi(self.lk) = (1 - self.reco_lk)*self.Fi(self.lk);
                    end
                end
            end
            if ~isempty(self.reco_hk) && ~isempty(self.hk)
                self.Bi(self.hk) = self.reco_hk*self.Fi(self.hk);
                if isempty(self.salidas)
                    self.Di(self.hk) = (1 - self.reco_hk)*self.Fi(self.hk);
                else
                    self.Di(self.hk) = (1 - self.reco_hk)*self.Fi(self.hk);
                end
            end
            if nargin>10
                try
                    self.xDi = self.Di(1)./(sum(self.Di));
                    self.xBi = self.Bi(1)./(sum(self.Bi));
                catch
                end
            end
            
            if ~isempty(self.cmpbinario)
                try
                    self.equilibrio;
                catch ME
                end
            end
        end
        function salida = set_reco_top(self, clave_liviano, recup_lk)
            if ischar(clave_liviano) || length(clave_liviano)== 1
                if ischar(clave_liviano)
                    %si el usuario introdujo el clave como un string
                    %con el id del clave deseado.
                    for i = 1:self.num_sust
                        if strcmpi(clave_liviano, self.comp(i))
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
                %Si usuario introduce varios valores numéricos
                self.lk = zeros(1, length(clave_liviano));
                for j=1:length(clave_liviano)
                    if ischar(clave_liviano{j})      
                        for i = 1:self.num_sust
                            if strcmpi(clave_liviano{j}, self.comp(i))
                                self.lk(j) = i; %se identifica ubicación del clave
                            end
                        end
                    else
                        self.lk = clave_liviano;
                        break
                    end
                end                          
                self.reco_lk = recup_lk;
            end
            salida = self;
        end
        function salida = set_reco_bot(self, clave_pesado, recup_hk)
            if ischar(clave_pesado) || length(clave_pesado)== 1
                if ischar(clave_pesado)
                    %si el usuario introdujo el clave como un string
                    %con el id del clave deseado.
                    for i = 1:self.num_sust
                        if strcmpi(clave_pesado, self.comp(i))
                            self.hk = i; %se identifica ubicación del clave
                        end
                    end
                else                        
                    %Si usuario introduce valor numérico que refleja la
                    %ubicación del clave liviano entre los
                    %ides de compuestos.
                    self.hk = clave_pesado;
                end
                self.reco_hk = recup_hk; %Se identifica la recuperación del clave
            else 
                %Si usuario introduce varios valores numéricos
                self.hk = zeros(1, length(clave_pesado));
                for j=1:length(clave_pesado)
                    if ischar(clave_pesado{j})      
                        for i = 1:self.num_sust
                            if strcmpi(clave_pesado{j}, self.comp(i))
                                self.hk(j) = i; %se identifica ubicación del clave
                            end
                        end
                    else
                        self.hk = clave_pesado;
                        break
                    end
                end                          
                self.reco_hk = recup_hk;
            end
            salida = self;
        end
        function salida = intervalos(self, N) %divide equilibrio en N intervalos
            salida = linspace(0, 1, N+1);
            self.equilibrio(salida);
        end
        function salida = equilibrio(self, equil)
            if nargin > 1
                self.E_x = equil;
                self.Tbin = self.E_x;
            end
            if isempty(self.E_x)
                self.E_x = linspace(0,1, 81);
                self.E_y = zeros(2, length(self.E_x));
                self.E_y(:,end) = [1;1];
                self.Tbin = self.E_x;
            end
            
            self.E_y = zeros(2, length(self.E_x));
            self.E_y(:,end) = [1;1];
            try
                self.Tbin(1) = self.MEdE.EdE.isofugT(self.perfil_p(2), self.cmpbinario(2), 1e-7);
                self.Tbin(end) = self.MEdE.EdE.isofugT(self.perfil_p(1), self.cmpbinario(1), 1e-7);
            catch ME
                tsat = self.cmpbinario(2).tsat;
                tsat = tsat{1};
                options = optimset('Display', 'none');
                self.Tbin(1) = fzero(@(T) tsat(T, self.perfil_p(end)), (self.cmpbinario(2).tcri - 1), options);
                tsat = self.cmpbinario(1).tsat;
                tsat = tsat{1};
                self.Tbin(end) = fzero(@(T) tsat(T, self.perfil_p(1)), (self.cmpbinario(1).tcri - 1), options);
            end
            %consigo la curva de equilibrio para 40 puntos
            for i = 2:2:length(self.E_x)-1
                mezcla = Mezcla(self.cmpbinario, self.E_x(i), self.klk);
                self.Tbin(i) = (self.Tbin(end) - self.Tbin(1))*self.E_x(i) + self.Tbin(1);
                [self.Tbin(i), self.E_y(:,i)] = self.MEdE.BubbleT(sqrt(self.perfil_p(1)*self.perfil_p(2)), mezcla, [], [], [], self.Tbin(i));
            end
            self.E_y = self.E_y(1,:);
            salida = self;
            % interpolo el resto de los puntos mediante Lagrange cuadrático
            for i = 3:2:length(self.E_x)-1
                self.E_y(i) = self.E_y(i-2)*(((self.E_x(i) - self.E_x(i-1))*(self.E_x(i)-self.E_x(i+1)))/((self.E_x(i-2) - self.E_x(i-1))*(self.E_x(i-2) - self.E_x(i+1)))) + self.E_y(i-1)*(((self.E_x(i) - self.E_x(i-2))*(self.E_x(i)-self.E_x(i+1)))/((self.E_x(i-1) - self.E_x(i-2))*(self.E_x(i-1) - self.E_x(i+1)))) + self.E_y(i+1)*(((self.E_x(i) - self.E_x(i-2))*(self.E_x(i)-self.E_x(i-1)))/((self.E_x(i+1) - self.E_x(i-2))*(self.E_x(i+1) - self.E_x(i-1))));
                self.Tbin(i) = self.Tbin(i-2)*(((self.E_x(i) - self.E_x(i-1))*(self.E_x(i)-self.E_x(i+1)))/((self.E_x(i-2) - self.E_x(i-1))*(self.E_x(i-2) - self.E_x(i+1)))) + self.Tbin(i-1)*(((self.E_x(i) - self.E_x(i-2))*(self.E_x(i)-self.E_x(i+1)))/((self.E_x(i-1) - self.E_x(i-2))*(self.E_x(i-1) - self.E_x(i+1)))) + self.Tbin(i+1)*(((self.E_x(i) - self.E_x(i-2))*(self.E_x(i)-self.E_x(i-1)))/((self.E_x(i+1) - self.E_x(i-2))*(self.E_x(i+1) - self.E_x(i-1))));
            end
            if self.Murphree ~= 1
                self.Mur_x = self.E_x;
                for i = 1:length(self.E_y)
                    self.Mur_y(i) = self.Murphree*(self.E_y(i) - self.E_x(i)) + self.E_x(i);
                end
            else
                self.Mur_x = self.E_x;
                self.Mur_y = self.E_y;
            end
        end
        function distributemin(self)
            iter = 0;

            self.intersectx = zeros(1, floor(length(self.rectaq)/2));
            self.intersectequilx = zeros(1, floor(length(self.rectaq)/2));
            self.intersecty = self.intersectx;
            self.intersectequil = self.intersectx;
            for i = 1:floor(length(self.rectaq)/2)
                iter = iter + 2;
                A = self.rectaq(1:2, iter - 1);
                B = self.rectaq(1:2, iter);
                A = A'; B = B';
                if abs(A(1) - B(1)) > 1e-3
                    options = optimset('Display', 'none');
                    max = find(self.Mur_x > self.entradas(floor(iter/2-1)*3+1), 1, 'first');
                    self.intersectx(floor(iter / 2), 1) = fsolve(@(x) interpolLagrange2(x, self.Mur_x, self.Mur_y, A, B), self.Mur_x(max), options);
                    self.intersectequilx(floor(iter / 2), 1) = fsolve(@(x) interpolLagrange2(x, self.E_x, self.E_y, A, B), self.E_x(max), options);
                    self.intersecty(floor(iter / 2), 1) = (A(2) - B(2))/(A(1) - B(1))*self.intersectx(floor(iter/2)) - (A(2) - B(2))/(A(1) - B(1))*A(1) + A(2);
                    self.intersectequil(floor(iter / 2), 1) = (A(2) - B(2))/(A(1) - B(1))*self.intersectequilx(floor(iter/2)) - (A(2) - B(2))/(A(1) - B(1))*A(1) + A(2);
                else
                    max = find(self.Mur_x > self.entradas(floor(iter/2-1)*3+1), 1, 'first');
                    self.intersectx(1, floor(iter / 2)) = A(1);
                    self.intersectequilx(1, floor(iter / 2)) = A(1);
                    %Interpolacion cuadrática de Lagrange
                    try
                        self.intersecty(1, floor(iter / 2)) = self.Mur_y(max-2)*(((A(1) - self.Mur_x(max-1))*(A(1) - self.Mur_x(max)))/((self.Mur_x(max - 2) - self.Mur_x(max - 1))*(self.Mur_x(max - 2) - self.Mur_x(max)))) + self.Mur_y(max-1)*(((A(1) - self.Mur_x(max-2))*(A(1) - self.Mur_x(max)))/((self.Mur_x(max-1) - self.Mur_x(max-2))*(self.Mur_x(max-1) - self.Mur_x(max)))) + self.Mur_y(max)*(((A(1) - self.Mur_x(max-2))*(A(1)-self.Mur_x(max-1)))/((self.Mur_x(max) - self.Mur_x(max-2))*(self.Mur_x(max) - self.Mur_x(max-1))));
                        self.intersectequil(1, floor(iter / 2)) = self.E_y(max-2)*(((A(1) - self.E_x(max-1))*(A(1) - self.E_x(max)))/((self.E_x(max - 2) - self.E_x(max - 1))*(self.E_x(max - 2) - self.E_x(max)))) + self.E_y(max-1)*(((A(1) - self.E_x(max-2))*(A(1) - self.E_x(max)))/((self.E_x(max-1) - self.E_x(max-2))*(self.E_x(max-1) - self.E_x(max)))) + self.E_y(max)*(((A(1) - self.E_x(max-2))*(A(1)-self.E_x(max-1)))/((self.E_x(max) - self.E_x(max-2))*(self.E_x(max) - self.E_x(max-1))));
                    catch
                        try 
                            self.intersecty(1, floor(iter / 2)) = self.Mur_y(max-1)*(((A(1) - self.Mur_x(max))*(A(1) - self.Mur_x(max+1)))/((self.Mur_x(max - 1) - self.Mur_x(max ))*(self.Mur_x(max - 1) - self.Mur_x(max+1)))) + self.Mur_y(max)*(((A(1) - self.Mur_x(max-1))*(A(1) - self.Mur_x(max+1)))/((self.Mur_x(max) - self.Mur_x(max-1))*(self.Mur_x(max) - self.Mur_x(max+1)))) + self.Mur_y(max+1)*(((A(1) - self.Mur_x(max-1))*(A(1)-self.Mur_x(max+1)))/((self.Mur_x(max+1) - self.Mur_x(max-1))*(self.Mur_x(max+1) - self.Mur_x(max))));
                            self.intersectequil(1, floor(iter / 2)) = self.E_y(max-1)*(((A(1) - self.E_x(max))*(A(1) - self.E_x(max+1)))/((self.E_x(max - 1) - self.E_x(max ))*(self.E_x(max - 1) - self.E_x(max+1)))) + self.E_y(max)*(((A(1) - self.E_x(max-1))*(A(1) - self.E_x(max+1)))/((self.E_x(max) - self.E_x(max-1))*(self.E_x(max) - self.E_x(max+1)))) + self.E_y(max+1)*(((A(1) - self.E_x(max-1))*(A(1)-self.E_x(max+1)))/((self.E_x(max+1) - self.E_x(max-1))*(self.E_x(max+1) - self.E_x(max))));
                        catch
                            try
                                self.intersecty(1, floor(iter / 2)) = self.Mur_y(max-3)*(((A(1) - self.Mur_x(max-2))*(A(1) - self.Mur_x(max-1)))/((self.Mur_x(max - 3) - self.Mur_x(max - 2))*(self.Mur_x(max - 3) - self.Mur_x(max-1)))) + self.Mur_y(max-2)*(((A(1) - self.Mur_x(max-2))*(A(2) - self.Mur_x(max-1)))/((self.Mur_x(max-2) - self.Mur_x(max-3))*(self.Mur_x(max-2) - self.Mur_x(max-1)))) + self.Mur_y(max-1)*(((A(1) - self.Mur_x(max-3))*(A(1)-self.Mur_x(max-2)))/((self.Mur_x(max-1) - self.Mur_x(max-3))*(self.Mur_x(max-1) - self.Mur_x(max-2))));
                                self.intersectequil(1, floor(iter / 2)) = self.E_y(max-3)*(((A(1) - self.E_x(max-2))*(A(1) - self.E_x(max-1)))/((self.E_x(max - 3) - self.E_x(max - 2))*(self.E_x(max - 3) - self.E_x(max-1)))) + self.E_y(max-2)*(((A(1) - self.E_x(max-2))*(A(2) - self.E_x(max-1)))/((self.E_x(max-2) - self.E_x(max-3))*(self.E_x(max-2) - self.E_x(max-1)))) + self.E_y(max-1)*(((A(1) - self.E_x(max-3))*(A(1)-self.E_x(max-2)))/((self.E_x(max-1) - self.E_x(max-3))*(self.E_x(max-1) - self.E_x(max-2))));
                            catch
                                self.intersecty(1, floor(iter / 2)) = self.Mur_y(max)*(((A(1) - self.Mur_x(max+1))*(A(1) - self.Mur_x(max+2)))/((self.Mur_x(max) - self.Mur_x(max+1 ))*(self.Mur_x(max ) - self.Mur_x(max+2)))) + self.Mur_y(max+1)*(((A(1) - self.Mur_x(max))*(A(1) - self.Mur_x(max+2)))/((self.Mur_x(max+1) - self.Mur_x(max))*(self.Mur_x(max+1) - self.Mur_x(max+2)))) + self.Mur_y(max+2)*(((A(1) - self.Mur_x(max))*(A(1)-self.Mur_x(max+2)))/((self.Mur_x(max+2) - self.Mur_x(max))*(self.Mur_x(max+2) - self.Mur_x(max+1))));
                                self.intersectequil(1, floor(iter / 2)) = self.E_y(max)*(((A(1) - self.E_x(max+1))*(A(1) - self.E_x(max+2)))/((self.E_x(max) - self.E_x(max+1 ))*(self.E_x(max ) - self.E_x(max+2)))) + self.E_y(max+1)*(((A(1) - self.E_x(max))*(A(1) - self.E_x(max+2)))/((self.E_x(max+1) - self.E_x(max))*(self.E_x(max+1) - self.E_x(max+2)))) + self.E_y(max+2)*(((A(1) - self.E_x(max))*(A(1)-self.E_x(max+2)))/((self.E_x(max+2) - self.E_x(max))*(self.E_x(max+2) - self.E_x(max+1))));
                            end
                        end
                    end
                end
            
            if ~isempty(self.salidas)
                for i=1:length(self.salidas)/3
                    if self.salidas(i*3-1) == 1;
                        indice = find(self.salidas(i*3-2) > self.intersectx, 1, 'first');
                        if ~isempty(indice);
                            self.rectasminimas();
                            iterint = 0;
                            for j = length(self.rectaq)+1:-2:2*(indice)
                                iterint = iterint + 1;
                                self.rectaq(:,j:j+1) = self.rectaq(:,j-2:j-1);
                            end
                            self.intersectx(end + 1) = self.intersectx(end + 1 - 1);
                            self.intersectequilx(end + 1) = self.intersectx(end + 1 - 1);
                            self.intersecty(end + 1) = self.intersecty(end + 1 - 1);
                            self.intersectequil(end + 1) = self.intersectequil(end + 1 - 1);
                            self.intersectx(end - 1) = self.intersectx(end - 2);
                            self.intersectequilx(end - 1) = self.intersectequilx(end - 2);
                            self.intersecty(end - 1) = self.intersecty(end - 2);
                            self.intersectequil(end - 1) = self.intersectequil(end - 2);
                            max = find(self.Mur_x > self.salidas(floor(i)*3-2), 1, 'first');
                            self.rectaq(:,indice*2-1:indice*2) = [self.salidas(i*3-2), self.salidas(i*3-2); self.salidas(i*3-2), 2];
                            self.intersectx(1,indice) = self.salidas(i*3-2);
                            self.intersectequilx(1,indice) = self.salidas(i*3-2);
                            self.intersecty(1,indice) = self.Mur_y(max-2)*(((self.intersectx(1,indice) - self.Mur_x(max-1))*(self.intersectx(1,indice) - self.Mur_x(max)))/((self.Mur_x(max - 2) - self.Mur_x(max - 1))*(self.Mur_x(max - 2) - self.Mur_x(max)))) + self.Mur_y(max-1)*(((self.intersectx(1,indice) - self.Mur_x(max-2))*(self.intersectx(1,indice) - self.Mur_x(max)))/((self.Mur_x(max-1) - self.Mur_x(max-2))*(self.Mur_x(max-1) - self.Mur_x(max)))) + self.Mur_y(max)*(((self.intersectx(1,indice) - self.Mur_x(max-2))*(self.intersectx(1,indice)-self.Mur_x(max-1)))/((self.Mur_x(max) - self.Mur_x(max-2))*(self.Mur_x(max) - self.Mur_x(max-1))));
                            self.intersectequil(1,indice) = self.E_y(max-2)*(((self.intersectequilx(1,indice) - self.E_x(max-1))*(self.intersectequilx(1,indice) - self.E_x(max)))/((self.E_x(max - 2) - self.E_x(max - 1))*(self.E_x(max - 2) - self.E_x(max)))) + self.E_y(max-1)*(((self.intersectequilx(1,indice) - self.E_x(max-2))*(self.intersectequilx(1,indice) - self.E_x(max)))/((self.E_x(max-1) - self.E_x(max-2))*(self.E_x(max-1) - self.E_x(max)))) + self.E_y(max)*(((self.intersectequilx(1,indice) - self.E_x(max-2))*(self.intersectequilx(1,indice)-self.E_x(max-1)))/((self.E_x(max) - self.E_x(max-2))*(self.E_x(max) - self.E_x(max-1))));
                        end
                    end
                end
            else
                self.rectasminimas();
            end
            reflujoguardado = self.reflujo;
            self.reflujo = 1e308;
            platosmin = self.platosmin();
            self.distribute([], 'off');
            self.Nmin = self.platosmin();
            self.platosmin = platosmin;
            self.reflujo = reflujoguardado;
            end
        end
        function mccabethiele(self, reflux)
            if nargin == 1
                %Si se desea hallar el mínimo de etapas
                yequil = self.xDi(1);
                xequil = self.xDi(1);
                if self.intersectx(end) ~= self.xBi(1)
                    self.intersectx(end + 1) = self.xBi(1);
                    self.intersectequilx(end + 1) = self.xBi(1);
                    self.intersecty(end + 1) = self.xBi(1);
                    self.intersectequil(end + 1) = self.xBi(1);
                end
                plato = 0;
                self.etapasemin = zeros(1, length(self.intersectx) - 1);
                if self.Murphree == 1
                    self.Mur_x = self.E_x;
                    self.Mur_y = self.E_y;
                end
                for i = 1:length(self.intersectx)
                    platoinicia = plato;
                    while yequil > self.intersecty(i)
                        plato = plato + 1;
                        A = [xequil, yequil];
                        B = [0, yequil];
                        max = find(self.Mur_x < xequil, 1, 'last');
                        options = optimset('Display', 'none');
                        if plato == 1 && self.condparcial
                            xequilattempt = fsolve(@(x) interpolLagrange2(x, self.E_x, self.E_y, A, B), self.E_x(max), options);
                            B = [xequilattempt, yequil];
                        else
                            xequilattempt = fsolve(@(x) interpolLagrange2(x, self.Mur_x, self.Mur_y, A, B), self.Mur_x(max), options);
                        
                        if xequilattempt < self.xBi(1)
                            xequilattempt = fsolve(@(x) interpolLagrange2(x, self.E_x, self.E_y, A, [0, yequil]), self.E_x(max), options);
                            if plato ~= 1
                                B = [xequilattempt, yequil];
                            end
                        end                        
                        xequil = xequilattempt;
                        B = [xequil, yequil];
                            end
                        if plato == 2
                            if self.Murphree == 1
                                
                                plot(B, [B(2), A(2)], 'Color', [.12,.12,.12], 'LineStyle', ':')
                            else
                                plot(B, [B(2), A(2)], 'Color', [.12,.12,.12], 'LineStyle', ':')
                            end
                        else
                            if self.Murphree == 1
                                plot(B, A, 'Color', [.12,.12,.12], 'LineStyle', ':')
                            else
                                plot(B, A, 'Color', [.12,.12,.12], 'LineStyle', ':')
                            end
                        end
                        srt = strcat(int2str(plato),'\rightarrow');
                        if plato == 1
                            if self.Murphree == 1
                                text(B(1), yequil, srt, 'HorizontalAlignment','right', 'Color', [.12,.12,.12], 'LineStyle', ':')
                            else
                                text(B(1), yequil, srt, 'HorizontalAlignment','right', 'Color', [.12,.12,.12], 'LineStyle', ':')
                            end
                        else
                            if self.Murphree == 1
                                text(xequil, yequil, srt, 'HorizontalAlignment','right', 'Color', [.12,.12,.12], 'LineStyle', ':')
                            else
                                text(xequil, yequil, srt, 'HorizontalAlignment','right', 'Color', [.12,.12,.12], 'LineStyle', ':')
                            end
                        end
                        tequil = find(self.Mur_y < xequil, 1, 'last');
                        tequil = self.Tbin(tequil);
                        mezclabinaria = Mezcla(self.cmpbinario, xequil, self.klk);
                        [self.TPlatobin(plato), ~] = self.MEdE.DewT(sqrt(self.perfil_p(1)*self.perfil_p(2)), mezclabinaria, [], [],[], tequil);
                        
                        if self.Murphree == 1
                            plot([xequil, xequil],[xequil,yequil], 'Color', [.12,.12,.12], 'LineStyle', ':')
                        else
                            plot([xequil, xequil],[xequil,yequil], 'Color', [.12,.12,.12], 'LineStyle', ':')
                        end
                        yequil = B(1);
                        if plato - platoinicia > 300
                            break
                        end
                    end
                    if plato - platoinicia > 300
                        error('No se pudo alcanzar una resolución porque se superaron 300 platos')
                    end
                    if i ~=length(self.intersectx)
                        if i == 1
                            self.etapasemin(i) = plato;
                        else
                            self.etapasemin(i) = plato - self.etapasemin(i-1);
                        end
                    end
                end
                if self.Murphree == 1
                    self.platosmin = plato;
                else
                    self.platosaefficmin = plato;
                end
                self.Nmin = plato;
            else    
                %Si se desea hallar el número real de etapas
                if self.Murphree==1
                    self.grafica('off');
                end
                yequil = self.xDi(1);
                xequil = self.xDi(1);
                if self.intersectx(end) ~= self.xBi(1)
                    self.intersectx(end + 1) = self.xBi(1);
                    self.intersecty(end + 1) = self.xBi(1);
                    self.intersectequil(end + 1) = self.xBi(1);
                end
                plato = 0;
                self.etapasemin = zeros(1, length(self.intersectx) - 1);
                self.rectas = zeros(length(self.intersectx) - 1, 2);
                for i = 1:length(self.intersectx)
                    if i == 1
                        liquido = reflux.*sum(self.Di);
                        vapor = (reflux + 1).*sum(self.Di);
                        b = self.xDi(1)/(reflux + 1);                        
                    end
                    if self.vivo && i == length(self.intersectx)
                        if abs(self.entradas((i-1)*3-1)-1) > 1e-3
                            x45 = (-(self.entradas((i-1)*3-2))/(self.entradas((i-1)*3-1)-1) -self.rectas(2,i-1))/(self.rectas(1,i-1) - ((self.entradas((i-1)*3-1))/(self.entradas((i-1)*3-1)-1)));
                        else
                            x45 = (self.entradas((i-1)*3-2));
                        end
                        y45 = self.rectas(2,i-1)+self.rectas(1,i-1)*x45;
                        m = (y45)/(x45 - self.intersectx(i));
                        b = 0 - m*self.xBi(1);
                    else
                        m = liquido / vapor;
                    end
                    self.rectas(1:2, i) = [m;b];
                    x = 0:0.5:1;
                    y = x*m + b;
                    plot(x, y);
                    platoinicia = plato;
                    while yequil > self.intersecty(i)
                        plato = plato + 1;
                        A = [xequil, yequil];
                        if plato == 1
                            B = [0, yequil];
                        end
                        max = find(self.Mur_x < xequil, 1, 'last');
                        options = optimset('Display', 'none');
                        if plato == 1 &&  self.condparcial
                            xequilattempt = fsolve(@(x) interpolLagrange2(x, self.E_x, self.E_y, A, B), self.E_x(max), options);
                        else
                            xequilattempt = fsolve(@(x) interpolLagrange2(x, self.Mur_x, self.Mur_y, A, [0, yequil]), self.Mur_x(max), options);
                        end
                        if plato == 1
                            B = [xequilattempt, yequil];
                        end
                        y1equil = yequil;
                        yequilattempt = xequilattempt*m + b;
                        if abs(self.Murphree-1)>1e-4
                            if yequilattempt < self.xBi(1)
                                xequilattempt = fsolve(@(x) interpolLagrange2(x, self.E_x, self.E_y, A, [0, yequil]), self.E_x(max), options);
                                if plato == 1
                                    B = [xequilattempt, yequil];
                                end
                                y1equil = yequil;
                                yequilattempt = xequil*m + b;
                            end     
                        end
                        xequil = xequilattempt;
                        yequil = yequilattempt;
                        if self.Murphree == 1
                            plot([A(1), xequil], [B(2),B(2)], 'Color', [.12,.12,.12], 'LineStyle', ':')
                        else
                            plot([A(1), xequil], [B(2),B(2)], 'Color', [.12,.12,.12], 'LineStyle', ':')
                        end
                        
                        B = [xequil, yequil];
                        srt = strcat(int2str(plato),'\rightarrow');
                        if self.Murphree == 1
                            text(xequil, y1equil, srt, 'HorizontalAlignment','right', 'Color', [.12,.12,.12], 'LineStyle', ':')
                        else
                            text(xequil, y1equil, srt, 'HorizontalAlignment','right', 'Color', [.12,.12,.12], 'LineStyle', ':')
                        end
                        tequil = find(self.Mur_y < xequil, 1, 'last');
                        tequil = self.Tbin(tequil);
                        mezclabinaria = Mezcla(self.cmpbinario, xequil, self.klk);
                        [self.TPlatobin(plato), ~] = self.MEdE.DewT(sqrt(self.perfil_p(1)*self.perfil_p(2)), mezclabinaria, [], [],[], tequil);
                        if self.Murphree == 1
                            plot([B(1), B(1)], [yequil, y1equil], 'Color', [.12,.12, .12], 'LineStyle', ':')
                        else
                            plot([B(1), B(1)], [yequil, y1equil], 'Color', [.12,.12, .12], 'LineStyle', ':')
                        end
                        
                        
                        if yequil < self.xBi(1) || yequilattempt < self.xBi(1) || xequil < self.xBi(1)
                            break
                        end
                        if plato - platoinicia > 300
                            break
                        end
                    end
                    if plato - platoinicia > 300
                        error('No se pudo alcanzar una resolución porque se superaron 300 platos')
                    end
                    if yequil < self.xBi(1)
                        break
                    end
                    if i ~=length(self.intersectx)
                        if i == 1
                            self.etapasemin(i) = plato;
                        else
                            self.etapasemin(i) = plato - self.etapasemin(i-1);
                        end
                    end
                    display(yequil)
                    if i ~= length(self.intersectx)-1 && i ~= length(self.intersectx)
                        try
                            indice = find(self.intersectx(i) == self.salidas, 1, 'first');
                            if isempty(indice)
                                error('No hay salidas')
                            end
                            liquido = liquido - self.salidas(indice + 1)*self.salidas(indice + 2);
                            vapor = vapor +  (1 - self.salidas(indice + 1))*self.salidas(indice + 2);
                            b = self.salidas(indice)*m + b - liquido / vapor * self.salidas(indice);
                        catch 
                            try 
                                entradai = entradai + 3;
                                entradaj = entradaj + 1;
                            catch
                                entradai = 1;
                                entradaj = 1;
                            end
                            indice = entradai;
                            
                            liquido = liquido + self.entradas(indice + 1)*self.entradas(indice + 2)*(self.alimentaciones(entradaj).conc(self.lk) + self.alimentaciones(entradaj).conc(self.hk));
                            vapor = vapor -  (1 - self.entradas(indice + 1))*self.entradas(indice + 2)*(self.alimentaciones(entradaj).conc(self.lk) + self.alimentaciones(entradaj).conc(self.hk));
                            xu = (-(self.entradas(indice + 1)*b - b + self.entradas(indice)))/((self.entradas(indice + 1)*(m - 1))-m);
                            if abs(self.entradas(indice + 1) - 1) < 1e-3
                                b = self.entradas(indice)*m + b - liquido/vapor*self.entradas(indice);
                            else
                                b = xu*m + b -  liquido/vapor*xu;
                            end
                        end
                    elseif  i == length(self.intersectx)-1
                        get_intery_point = self.rectastripmin(1);
                        try
                            indice = find(self.intersectx(i) == self.salidas, 1, 'first');
                            if isempty(indice)
                                error('No hay salidas')
                            end
                                liquido = liquido - (self.salidas(indice + 1))*self.salidas(indice + 2);
                                vapor = vapor + ((1 - self.salidas(indice + 1))*self.entradas(indice + 2));
                                m1 = liquido / vapor;
                            if self.vivo 
                                b = 0 - m1.*self.xBi(1);
                            else
                                b = self.xBi(1)-m1*self.xBi(1);
                            end
                            m = m1;
                        catch 
                            try 
                                entradai = entradai + 3;
                                entradaj = entradaj + 1;
                            catch
                                entradai = 1;
                                entradaj = 1;
                            end
                                
                                liquido = liquido + ((self.entradas(entradai + 1))*self.entradas(entradai + 2))*(self.alimentaciones(entradaj).conc(self.lk) + self.alimentaciones(entradaj).conc(self.hk));
                                vapor = vapor - ((1-self.entradas(entradai + 1))*self.entradas(entradai + 2))*(self.alimentaciones(entradaj).conc(self.lk) + self.alimentaciones(entradaj).conc(self.hk));
                            	m1 = liquido / vapor;
                            if self.vivo 
                                b = 0 - m1*self.xBi(1);
                            else
                                b = self.xBi(1)-m1*self.xBi(1);
                            end
                            m = m1;
                        end
                    end
                    if yequil < self.xBi(1)
                        break
                    end
                end
                if self.Murphree == 1
                    self.platosmin = plato;
                else
                    self.etapas = plato;
                end
            end
                
        end
        function rectasminimas(self)
            if isempty(self.xDi)
                self.xDi = self.Di./(sum(self.Di));
            end
            if isempty(self.xBi)
                self.xBi = self.Bi./(sum(self.Bi));
            end
            self.rectarectmin = [self.intersectx(1), self.xDi(1); self.intersecty(1), self.xDi(1)];
            m = (self.rectarectmin(2, 1) - self.rectarectmin(2, 2))/(self.rectarectmin(1,1) - self.rectarectmin(1,2));
            b = self.intersecty(1) - m*self.intersectx(1);
            self.rectarectmin = [0, 1; b, m*1 + b];
            if isempty(self.Rmin)
                self.Rmin = -m /(m - 1);
            end
            self.intersectx = self.intersectx( self.intersectx~=0)';
            self.intersecty = self.intersecty( self.intersecty~=0)';
            
            if self.intersectx(end) ~= self.xBi(1)
                self.intersectx(end + 1) = self.xBi(1);
                self.intersecty(end + 1) = self.xBi(1);
            end
            if self.intersectx(end) ~= self.xBi(1)
                self.intersecty(end + 1) = self.xBi(1);
            end
            self.rectasminima(1:2,1) = [m; b];
            if self.vivo == false
                self.rectastripmin = [ self.intersectx(end-1), self.xBi(1); self.intersecty(end-1), self.xBi(1)];
                m = (self.rectastripmin(2, 1) - self.rectastripmin(2, 2))/(self.rectastripmin(1,1) - self.rectastripmin(1,2));
                b = self.xBi(1) - m*self.intersectx(end);
            else 
                self.rectastripmin = [self.intersectx(end-1), self.xBi(1); self.intersecty(end-1), 0];
                m = (self.rectastripmin(2, 1) - self.rectastripmin(2, 2))/(self.rectastripmin(1,1) - self.rectastripmin(1,2));
                b = 0 - m*self.intersectx(end);
            end
            self.rectastripmin = [0, 1; b, m*1 + b];
            self.rectasminima(1:2,2) = [m; b];
            if isempty(self.vap_min)
                self.vap_min = m./(m-1);
            end
            if length(self.alimentaciones) > 1
                Vfondo = self.vap_min*self.bflujo;
                Vtope = self.dflujo*(self.Rmin + 1);
                for i = 2 :3:length(self.entradas)
                    Vfondo = Vfondo + (1-self.entradas{i}).*self.entradas{i+1};
                end
                if Vfondo > Vtope
                    self.Rmin = (Vfondo / self.dflujo)-1;
                    fprintf(1, 'Tuvo que ser reajustado el reflujo de tope mínimo por ser\n menor al necesario según recta mínima inferior');
                end
            end
        end
        function reflux_remin(self, reflujo2rmin)
            if ~isempty(self.Rmin) && nargin > 1
                self.reflujo = self.Rmin.*reflujo2rmin;
            end
        end
        function distribute(self, reflux, onoff)
            if nargin > 1 && ~isempty(reflux)
                self.reflujo = reflux;
            end
            if nargin < 3 || ~strcmpi(onoff, 'off')
                self.distributemin();
            end
            if self.metodo == 1
                if self.Murphree ~= 1
                    effic = self.Murphree;
                    self.Murphree = 1;
                    mur_x = self.Mur_x;
                    mur_y = self.Mur_y;
                    if self.intersectequil(end) ~= self.xBi(1)
                        self.intersectequil(end + 1) = self.xBi(1);
                    end
                    intersecty = self.intersecty;
                    self.intersecty = self.intersectequil;
                    self.Mur_x = self.E_x;
                    self.Mur_y = self.E_y;
                    self.mccabethiele(self.reflujo);
                    self.Murphree = effic;
                    self.Mur_x = mur_x;
                    self.Mur_y = mur_y;
                    self.intersecty = intersecty;
                    self.grafica()
                end
                self.mccabethiele(self.reflujo);
            else
                self.ponchonsavarit();
            end
        end
        function grafica(self, off)
            figure            
            hold on
            if nargin > 1 && ~isempty(off)
                hold off
                grid off
            end
            try 
                xi = 0:0.0005:1;
                yi=interp1(self.E_x,self.E_y,xi,'spline');
                plot([self.E_x(1), self.E_x(end)], [self.E_x(1), self.E_x(end)],'Color', 'black','LineWidth',0.5, 'LineStyle', '-')
                hold on
                plot(xi, yi,'Color', 'red', 'LineWidth',1, 'LineStyle', '-.');
                legend('', 'Equilibrio binario de LK y HK','Equilibrio según Eficiencia de Murphree','Location', 'southeast');
            catch ME
                error('No se tiene la curva de equilibrio construida')
            end
            hold on
            if self.Murphree ~= 1
                yi=interp1(self.Mur_x,self.Mur_y,xi,'spline');
                plot(xi, yi, 'Color', [0.75, 0.375, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
                legend('' , 'Equilibrio binario de LK y HK','Equilibrio según Eficiencia de Murphree', 'Location', 'southeast')                
            end
            try 
                for i = 1:2:length(self.rectaq)
                    plot(self.rectaq(1, i:i+1), self.rectaq(2, i:i+1), ':', 'Color', [0.1875, 0.5625, 0.5], 'MarkerEdgeColor', 'black',   'LineWidth', 1.5);
                end
                legend('', 'Equilibrio binario de LK y HK','Equilibrio según Eficiencia de Murphree','Alimentaciones', 'Location', 'southeast')
            catch ME 
                error('No se tiene la recta de la alimentación construida')
            end
            axis([-0.0002 1.0002 -0.0002 1.0002])
        end
        function graficamin(self, off)
            figure
            hold on
            grid off
            if nargin > 1 && ~isempty(off)
                hold off
                grid off
            end
            try 
                xi = 0:0.0005:1;
                yi=interp1(self.E_x,self.E_y,xi,'spline');
                plot([self.E_x(1), self.E_x(end)], [self.E_x(1), self.E_x(end)],'Color', 'black','LineWidth',0.5, 'LineStyle', '-')
                hold on
                plot(xi, yi,'Color', 'red', 'LineWidth',1, 'LineStyle', '-.');
                legend('', 'Equilibrio binario de LK y HK','Location', 'southeast');
            catch ME
                error('No se tiene la curva de equilibrio construida')
            end
            hold on
            if isempty(self.Murphree) || self.Murphree == 1
                yi=interp1(self.Mur_x,self.Mur_y,xi,'spline');
                plot(xi, yi, 'Color', [0.5, 0.375, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
                legend('' , 'Equilibrio binario de LK y HK', 'Location', 'southeast')                
            else
                yi=interp1(self.Mur_x,self.Mur_y,xi,'spline');
                plot(xi, yi, 'Color', [0.5, 0.375, 0.5], 'LineWidth', 0.5, 'LineStyle', '--');
                legend('' , 'Equilibrio binario de LK y HK','Equilibrio según Eficiencia de Murphree', 'Location', 'southeast')                
            end
            try 
                for i = 1:2:length(self.rectaq)
                    plot(self.rectaq(1, i:i+1), self.rectaq(2, i:i+1), ':', 'Color', [0.1875, 0.5625, 0.5], 'MarkerEdgeColor', 'black', 'LineWidth', 1.5);
                end
                if isempty(self.Murphree) || self.Murphree == 1
                    legend('', 'Equilibrio binario de LK y HK','Alimentaciones', 'Location', 'southeast')
                else
                    legend('', 'Equilibrio binario de LK y HK','Eficiencia según equilibrio de Murprhee', 'Alimentaciones', 'Location', 'southeast')
                end
            catch ME 
                error('No se tiene la recta de la alimentación construida')
            end
            axis([-0.0002 1.0002 -0.0002 1.0002])
            try 
                plot(self.rectarectmin(1,1:2), self.rectarectmin(2,1:2))
            catch ME
                error('No se tiene la recta de rectificación mínima construida')
            end
            try 
                plot(self.rectastripmin(1,1:2), self.rectastripmin(2,1:2))
            catch ME
                error('No se tiene la recta de despojo mínima inferior construida')
            end
            if self.Murphree ~= 1
                effic = self.Murphree;
                self.Murphree = 1;
                self.mccabethiele();
                self.Murphree = effic;
                self.mccabethiele();
            else 
                self.mccabethiele();
            end
        end
        function alter_define(self, alim, lightk, recup_top, heavyk, recup_bot, equil)
            long = length(alim);
            self.entradas = zeros(1, long*3);     
            for i = 1:1:long
                corr = alim(i);
                susta = corr.comp;
                for j = 1:length(susta)
                    if ~any(strcmp(self.comp, susta(j).id))
                        self.comp{length(self.comp) + 1} = susta(j).id;
                        self.sust(length(self.sust) + 1) = susta(j);
                    end
                    try    
                        if lightk == j && strcmpi(susta(j).id, self.sust(lightk).id)
                            self.entradas(3.*i - 2) = corr.conc(j)/(corr.conc(j)+corr.conc(j+1));
                            if isempty(self.cmpbinario)
                                self.cmpbinario = Sustancia.empty(0,2);
                                self.cmpbinario(1) = susta(j);
                            end
                        elseif strcmp(lightk, susta(j).id)
                            self.entradas(3.*i - 2) = corr.conc(j)/(corr.conc(j)+corr.conc(j+1));
                            if isempty(self.cmpbinario)
                                self.cmpbinario = Sustancia.empty(0,2);
                                self.cmpbinario(1) = susta(j);
                            end
                        end
                    catch
                        
                    end
                    try
                        if heavyk == j
                            if length(self.cmpbinario) == 1
                                self.cmpbinario(2) = susta(j);
                                if isempty self.klk
                                    klkhk = strcat('k', int2str(lightk), int2str(heavyk));
                                    self.klk = alim(1).mezcla.kij(find(strcmpi(alim(1).mezcla.kbinario, klkhk),1,'first'));
                                end
                            end
                        elseif strcmp(heavyk, susta(j).id)
                            if length(self.cmpbinario) == 1
                                self.cmpbinario(2) = susta(j);
                                if isempty self.klk
                                    klkhk = strcat('k', int2str(lightk), int2str(heavyk));
                                    self.klk = alim(1).mezcla.kij(find(strcmpi(alim(1).mezcla.kbinario, klkhk),1,'first'));
                                end
                            end
                        end
                    catch
                        continue
                    end
                end
                self.entradas(3*(i-1) + 2) = alim(i).q;
                self.entradas(3*(i-1) + 3) = alim(i).molF;
            end
            if isempty(self.num_sust)
                self.num_sust = length(self.comp);
            end
            self.concfeed = zeros(long, length(self.comp));
            self.dsubi = zeros(1, length(self.comp));
            self.bsubi = zeros(1, length(self.comp));
            for i = 1:1:long
                corr = alim(i);
                susta = corr.comp;
                for j = 1:length(susta)
                    indice = find(strcmp(susta(j).id, self.comp),1, 'first');
                    self.concfeed(i, indice) = corr.conc(j);
                end
            end
            self.Fi = zeros(self.num_sust, 1);
            for l = 1:length(self.comp)
                for u = 3:3:length(self.entradas)
                    self.Fi(l) = self.Fi(l) + self.entradas(u).*self.concfeed(u/3, l);
                end
            end
            if isempty(self.alimentaciones)
                self.alimentaciones = alim;
            end
            if isempty(self.MEdE)
                self.MEdE = alim(1).MEdE;
            end
            if isempty(self.rectaq)
                self.rectaq = zeros(2, (length(self.entradas)/3)*2);
            end
            iter = 0;
            for i = 2:3:length(self.entradas)
                iter = iter + 1;
                self.rectaq(:, (iter)*2-1) = [self.entradas(i - 1), self.entradas(i - 1)];
                if abs(self.entradas(i) - 1) < 1e-5
                    self.rectaq(:, (iter)*2) = [self.entradas(i - 1), 2];
                elseif abs(self.entradas(i)) < 1e-5
                    self.rectaq(:, (iter)*2) = [-2, self.entradas(i - 1)];
                else
                    self.rectaq(:, (iter)*2) = [0, -(self.entradas(i - 1)/(self.entradas(i)-1))];
                end
            end
            if nargin > 3   
                if ~isempty(recup_top) && ~isempty(lightk)
                    self.set_reco_top(lightk, recup_top);
                end
            end
            if nargin >  5
                if ~isempty(recup_bot) && ~isempty(heavyk)
                    self.set_reco_bot(heavyk, recup_bot);
                end
            end
            if nargin > 6 && ~isempty(equil)
                self.equilibrio(equil)
            end
        end
    end
    
end

