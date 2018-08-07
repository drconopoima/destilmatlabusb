classdef BubblePoint < handle
% %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENGLISH EXPLANATION
%
%    BUBBLEPOINT Rigorous Bubble Point method for solving complex distillation columns. Its use is best between
%       components that have close volatilities between them
%
%       References:
%           - Seader, Henley, Roper. "Separation Process Principles. 3rd
%           Edition.
%       Luis Jesús Díaz Manzo
%
% Separation of a binary mixture by rigorous methods
%
% A binary mixture of Butane, Pentane, 11% molar and with a binary interaction parameter of 0.0011
%
% mixture1 = Mixture ([Substance('Butane'), Substance('Pentane')],[0.11,0.89],[1.1e-3]);
%
% A 100 PSIA, saturated liquid
%
% current1 = current (mix1, 0, 'x', 689.475728, 'P', 100, 'm', RMVdW(PREdE));
%
% BUBBLE POINT
%
% By the Bubble Point method, only the reflux and the flow rate of distillate or waste is available as a single set of variables.
%
% In this example, it is considered that with 8 stages, reflux of 3 and 60 moles of background flow a residue is obtained.
% high in pentane purity
%
% Tower2 = BubblePoint(8, {4, current1}, Current.empty(0,1), Current.empty(0,1), {1,1,40}, 689.475728, 0, 3, 60, `B', 100, `C-001');
% This represents the following entry conditions
%         BubblePoint(Plates, Feed, Distilled currents and residue, Stage 1 output, q=1 liquid 40 moles, pressure, profile_q = 0 adiabatic, reflux=3, 60 = `B' background flow, 100 maximum iterations, name )
% Tower2.wanghenke();
%
% The method converges to the results:
%
% Distillate molar flows
% Tower2.dsubi = 10.4444, 29.5556
% Bottom residual molar flows
% Tower2.bsubi = 0.5575, 59.4425
%
% Temperature profile:
%
% Tower2.profile_t = 365.2696, 371.7052, 374.7180, 375.9715, 377.3161, 378.2996, 378.9864, 379.4471
%
% Steam profile:
%
% Tower2.profile_v = 0 , 160 , 160.6258, 158.4190, 155.0443, 162.8477, 171.4716 180.6595
%
% Liquid profile:
%
% Tower2.profile_l = 120, 120.6258, 118.4190, 215.0443, 222.8477, 231.4716, 240.6595, 60
%
% Profiles by individual components: Tower2.profile_vi and Tower2.profile_li
%
% The weakness of the method is that it does not converge energy balances, it only meets mass balances.
% %
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPLICACIÓN EN ESPAÑOL
%
% BUBBLEPOINT Método riguroso de Bubble Point para resolución de columnas de destilación complejas. Su uso es mejor entre
%       componentes que tienen volatilidades cercanas entre ellas
%
%       Referencias:
%           - Seader, Henley, Roper. "Separation Process Principles". 3rd
%           Edition.
%       Luis Jesús Díaz Manzo
%
% Separación de una mezcla binaria por métodos rigurosos
%
% Una mezcla binaria de Butano, Pentano, 11% molar y con un parámetro de interacción binaria de 0.0011
%
% mezcla1 = Mezcla ([Sustancia('Butane'), Sustancia('Pentane')], [0.11,0.89], [1.1e-3]);
%
% A 100 PSIA, en líquido saturado
%
% corriente1 = Corriente(mezcla1, 0, 'x', 689.475728, 'P', 100, 'm', RMVdW(PREdE));
%
% BUBBLE POINT
%
% Por el método de Bubble Point sólo existe como único conjunto de variables el reflujo, y el caudal de destilado o de residuo
%
% En este ejemplo se considera que con 8 etapas, reflujo de 3 y 60 moles de flujo de fondo se obtiene un residuo
% alto en pureza de pentano
%
% Torre2 = BubblePoint(8, {4, corriente1}, Corriente.empty(0,1), Corriente.empty(0,1), {1,1,40}, 689.475728, 0, 3, 60, 'B', 100, 'C-001');
% Lo que representa las siguientes condiciones de entrada
%         BubblePoint(Platos,  Alimentacion, Corrientes destilado y residuo,Salida etapa1, q=1 líquida 40 moles, presion, perfil_q = 0 adiabático, reflujo=3, 60 = 'B' flujo fondo, 100 maximo iteraciones, nombre   )
% Torre2.wanghenke();
%
% El método converge a los resultados:
%
% Flujos molares de destilado
% Torre2.dsubi = 10.4444, 29.5556
% Flujos molares de residuo de fondo
% Torre2.bsubi = 0.5575,  59.4425
%
% Perfil de temperaturas:
%
% Torre2.perfil_t = 365.2696, 371.7052, 374.7180, 375.9715, 377.3161, 378.2996, 378.9864, 379.4471
%
% Perfil de vapor:
%
% Torre2.perfil_v = 0 , 160 , 160.6258, 158.4190, 155.0443,  162.8477, 171.4716 180.6595
%
% Perfil de líquido:
%
% Torre2.perfil_l = 120, 120.6258, 118.4190, 215.0443, 222.8477, 231.4716, 240.6595, 60
%
% Perfiles por componentes individuales: Torre2.perfil_vi y Torre2.perfil_li
%
% La debilidad del método es que no converge los balances de energía, únicamente cumple los balances de masa.
% %

    properties
        id
        etapas
        perfil_p
        perfil_q
        perfil_t
        perfil_v
        perfil_l
        perfil_vi
        perfil_li
        perfil_k
        tol
        qc
        qb
        dflujo
        bflujo
        reflujo
        num_sust
        alimentaciones
        iteraciones = 100;
        cdest
        cfond
        salidas
        entradas
        gdlibertad
        comp
        sust
        MEdE
        dsubi
        bsubi
        ysubj
        xsubj
        concfeed
        matriz_tridiag
        lk
        hk
        recover
        Fi
        UWi
        platos
        respaldoiter = 0
        laststep = logical(0)
        diferencia
        actualiter
        respaldovalI
        respaldoV
        respaldoL
        respaldovi
        respaldoli
        respaldot
        respaldoxi
        respaldoyi
        respaldoqc
        respaldoqb
        actualvalI
    end

    methods
        function self = BubblePoint(platos, alim, dest, bot, salid, Pperfil, Qperfil, reflux, dobflow, d_o_b, iteraciones,  nombre)
            if nargin > 11&& ~isempty(nombre) && isa(nombre, 'char')
                self.id = nombre;
            end
            if nargin > 0 && ~isempty(platos)
                self.etapas = platos;
            end
            if nargin > 10 && ~isempty(iteraciones)
                self.iteraciones = iteraciones;
            else
                self.iteraciones = 100;
            end

            self.comp = cell.empty(0,0);
            self.sust = Sustancia.empty(0,0);
            if nargin > 1 && ~isempty(alim) && isa(alim, 'cell') && isa(alim{2}, 'Corriente')
                long = length(alim);
                self.entradas = cell.empty(0,0);
                for i = 2:2:long
                    corr = alim{i};
                    susta = corr.comp;
                    for j = 1:length(susta)
                        if ~any(strcmp(self.comp, susta(j).id))
                            self.comp{length(self.comp) + 1} = susta(j).id;
                            self.sust(length(self.sust) + 1) = susta(j);
                        end
                    end
                    self.entradas{length(self.entradas) + 1} = alim{i-1};
                    self.entradas{length(self.entradas) + 1} = alim{i}.q;
                    self.entradas{length(self.entradas) + 1} = alim{i}.molF;
                end
                self.num_sust = length(self.comp);
                self.concfeed = zeros(floor(long/2), length(self.comp));
                self.dsubi = zeros(1, length(self.comp));
                self.bsubi = zeros(1, length(self.comp));
                for i = 2:2:long
                    corr = alim{i};
                    susta = corr.comp;
                    for j = 1:length(susta)
                        indice = find(strcmp(susta(j).id, self.comp),1, 'first');
                        self.concfeed(floor(i/2), indice) = corr.conc(j);
                    end
                end
                self.xsubj = zeros(platos, length(self.comp));
                self.ysubj = zeros(platos, length(self.comp));
                self.alimentaciones = alim;
                self.MEdE = alim{2}.MEdE;
                if length(self.entradas) > 3
                    for i = 4:3:length(self.entradas)
                        savedstage = self.entradas{i};
                        savedq = self.entradas{i+1};
                        savedf = self.entradas{i+2};
                        savedalim = self.alimentaciones{floor((i-1)/3)*2+2};
                        if savedstage < self.entradas{i-3}
                            self.alimentaciones{floor((i-1)/3)*2+1} = self.entradas{i-3};
                            self.alimentaciones{floor((i-1)/3)*2-1} = savedstage;
                            self.entradas{i} = self.entradas{i-3};
                            self.entradas{i+1} = self.entradas{i-2};
                            self.entradas{i+2} = self.entradas{i-1};
                            self.entradas{i-3} = savedstage;
                            self.entradas{i-2} = savedq;
                            self.entradas{i-1} = savedf;
                            self.alimentaciones{floor((i-1)/3)*2+2} = self.alimentaciones{floor((i-1)/3)*2};
                            self.alimentaciones{floor((i-1)/3)*2} = savedalim;
                        end
                    end
                end
                self.Fi = zeros(1, self.num_sust);
                for l = 2:2:length(self.alimentaciones)
                    self.Fi = self.Fi + self.alimentaciones{l}.molF * self.alimentaciones{l}.conc;
                end
            end
            if nargin > 2 && isa(dest, 'Corriente')
                self.cdest = dest;
            end
            if nargin > 3 && isa(bot, 'Corriente')
                self.cfond = bot;
            end
            if nargin > 4 && ~isempty(salid) && isa(salid, 'cell')
                self.salidas = salid;
            end
            if nargin > 5 && isempty(Pperfil)
                self.perfil_p = ones(1, platos).*alim{2}.P;
            elseif nargin > 5 && isa(Pperfil, 'cell')
                %Si se provee una celda, se supone que cada 2do elemento contiene presi�n y cada elem. 1ero
            %contiene un n�mero de plato, se tiene que proveer necesariamente presi�n del condensador y el rehervidor
                self.perfil_p = zeros(1, platos);
                self.perfil_p(1) = Pperfil{2};
                self.perfil_p(end) = Pperfil{end};
                for i = 2:platos - 1
                    if any(cell2mat(Pperfil) == i)
                        self.perfil_p(floor(find(cell2mat(Pperfil)==1, 1, 'first')/2)) = Pperfil{i};
                    end
                end
                deltaP = Pperfil{end} - Pperfil{2};
                for i = 2 : platos - 1
                    if ~isempty(self.perfil_p(i))
                        self.perfil_p(i) = ((i-1)*(deltaP))/(platos - 1) + Pperfil{2};
                    end
                end
            elseif nargin > 5 && length(Pperfil) == platos
                self.perfil_p = Pperfil;
            elseif nargin > 5 && length(Pperfil) == 1
                self.perfil_p = ones(1, platos).*Pperfil;
            end
            if nargin > 6 && isempty(Qperfil)
                self.perfil_q = zeros(1, platos-2);
            elseif nargin > 6 && length(Qperfil) == 1
                self.perfil_q = ones(1,platos-2).*Qperfil;
            elseif nargin > 6 && length(Qperfil) == platos - 2
                self.perfil_q = Qperfil;
            end
            if nargin > 7 && ~isempty(reflux)
                self.reflujo = reflux;
            end
            if nargin > 8 && ~isempty(dobflow) && strcmpi(d_o_b, 'D')
                self.dflujo = dobflow;
            elseif nargin > 8 && ~isempty(dobflow) && strcmpi(d_o_b, 'B')
                self.bflujo = dobflow;
            end
        end
        function self = balanmasa(self)
            if ~any(self.bsubi ~= 0) && ~any(self.dsubi ~= 0)
                long = length(self.entradas);
                entran = zeros(1, floor(long/3));
                j = 0;
                for i = 3:3:long
                    j = j + 1;
                    entran(j) = self.entradas{i};
                end
                %Balance de masa global
                if ~isempty(self.dflujo) && isempty(self.bflujo) && isempty(self.salidas)
                    self.bflujo =  sum(entran) - self.dflujo;
                elseif ~isempty(self.bflujo) && isempty(self.dflujo) && isempty(self.salidas)
                    self.dflujo = sum(entran) - self.bflujo;
                elseif ~isempty(self.dflujo) && isempty(self.bflujo) && ~isempty(self.salidas)
                    long = length(self.salidas);
                    salen = zeros(1, floor(long/3));
                    j = 0;
                    for i = 3:3:long
                        j = j + 1;
                        salen(j) = self.salidas{i};
                    end
                    self.bflujo =  sum(entran) - self.dflujo - sum(salen);
                elseif ~isempty(self.bflujo) && isempty(self.dflujo) && ~isempty(self.salidas)
                    long = length(self.salidas);
                    salen = zeros(1, floor(long/3));
                    j = 0;
                    for i = 3:3:long
                        j = j + 1;
                        salen(j) = self.salidas{i};
                    end
                    self.dflujo =  sum(entran) - self.bflujo - sum(salen);
                end
                %Balance de masa por componente
                iter = 0;
                valI = 100;
                tamano = size(self.entradas);
                    iter = iter + 1;
                    salidasi = zeros(1, length(self.comp));
                    self.UWi = salidasi;
                    self.xsubj(2:self.etapas - 1,:) = zeros(self.etapas - 2, self.num_sust);
                    self.ysubj(2:self.etapas - 1,:) = zeros(self.etapas - 2, self.num_sust);
                    dsubi = zeros(1, length(self.dsubi));
                    bsubi = zeros(1, length(self.dsubi));
                    self.lk = [];
                    self.hk = [];
                    for l = 1:length(self.comp)
                        for m = 3:3:tamano(2)
                            try
                                if abs(self.dflujo/self.entradas{m}) > 1e-3
                                    %Suponiendo la distribuci�n de los
                                    %claves m�s livianos no ocurre y salen
                                    %completos por el tope
                                    if (sum((dsubi(l) + self.concfeed(m/3, l)*self.entradas{m}))+sum(dsubi(1:l-1))+sum(dsubi(l+1:end))) < self.dflujo
                                        dsubi(l) = dsubi(l) + self.concfeed(m/3, l).*self.entradas{m};
                                    else
                                        self.lk = l;
                                        if self.lk == self.num_sust
                                            self.lk = l-1;
                                            dsubi(l) = 0;
                                        elseif m~=3
                                            dsubi(l) = 0;
                                        end
                                        Flk = 0;
                                        Fhk = 0;
                                        for i = 3:3:tamano(2)
                                            Flk = Flk + self.concfeed(i/3, self.lk).*self.entradas{i};
                                        end
                                        Flk = Flk - salidasi(self.lk);
                                        for i = 3:3:tamano(2)
                                            Fhk = Fhk + self.concfeed(i/3, self.lk+1).*self.entradas{i};
                                        end
                                        Fhk = Fhk - salidasi(self.lk+1);
                                        break
                                    end
                                else
                                    if (sum((dsubi(l) + self.concfeed(m/3, l)*self.entradas{m}))+sum(dsubi(1:l-1))+sum(dsubi(l+1:end))) < self.salidas{3}
                                        dsubi(l) = dsubi(l) + self.concfeed(m/3, l).*self.entradas{m};
                                    else
                                        self.lk = l;
                                        if self.lk == self.num_sust
                                            self.lk = l-1;
                                            dsubi(l) = 0;
                                        elseif m~=3
                                            dsubi(l) = 0;
                                        end
                                        Flk = 0;
                                        Fhk = 0;
                                        for i = 3:3:tamano(2)
                                            Flk = Flk + self.concfeed(i/3, self.lk).*self.entradas{i};
                                        end
                                        Flk = Flk - salidasi(self.lk);
                                        for i = 3:3:tamano(2)
                                            Fhk = Fhk + self.concfeed(i/3, self.lk+1).*self.entradas{i};
                                        end
                                        Fhk = Fhk - salidasi(self.lk+1);
                                        break
                                    end
                                end
                            catch
                                self.lk = l-1;
                                dsubi(l-1) = 0;
                                Flk = 0;
                                Fhk = 0;
                                for i = 3:3:tamano(2)
                                    Flk = Flk + self.concfeed(i/3, self.lk).*self.entradas{i};
                                end
                                Flk = Flk - salidasi(self.lk);
                                for i = 3:3:tamano(2)
                                    Fhk = Fhk + self.concfeed(i/3, self.lk+1).*self.entradas{i};
                                end
                                Fhk = Fhk - salidasi(self.lk+1);
                                break
                            end
                        end
                        if ~isempty(self.lk)
                            break
                        end
                        dsubi(l) = dsubi(l) - salidasi(l);
                    end
                    self.hk = self.lk + 1;
                    if isempty(self.salidas) || ~isa(self.salidas, 'cell')
                        dsubi(self.lk) = self.dflujo - sum(dsubi);
                    elseif isempty(self.dflujo) && self.salidas{1} == 1
                        dsubi(self.lk) = self.salidas{3} - sum(dsubi);
                        self.dflujo = 0;
%                     end
                    elseif ~isempty(self.dflujo) && self.salidas{1} == 1
                        if dsubi(self.lk) ~= 0
                            dsubi(self.hk) = self.salidas{3}  + self.dflujo - sum(dsubi);
                        else
                            dsubi(self.lk) = self.salidas{3}  + self.dflujo - sum(dsubi);
                        end
                    end
%                     bsubi(self.hk) = Fhk - dsubi(self.hk);
                    bsubi = self.Fi - dsubi;

                    %dsubi = [160, 365.39, 4.61, 1e-10, 1e-14];
                    %bsubi = [1e-11, 4.61, 235.39, 25, 5 ];
                    if isempty(self.perfil_t) || all(self.perfil_t == 0)
                        if all(bsubi == 0)
                            bsubi = self.Fi - dsubi;
                        end
                        self.xsubj(end, :) = bsubi./sum(bsubi);
                        self.perfil_t = zeros(1, self.etapas);
                        if ~isempty(self.salidas) && isa(self.salidas, 'cell')
                            if self.salidas{1} == 1 && self.salidas{2} == 1
                                if ~isempty(self.salidas) && isa(self.salidas, 'cell')
                                    if self.salidas{1} == 1 && self.salidas{2} == 1
                                        if self.dflujo ~= 0 && self.dflujo/(self.dflujo + self.salidas{3}) > 1e-4
                                            zsubj = dsubi./sum(dsubi);
                                            mezcla = Mezcla(self.sust, zsubj);
                                            [Tb, self.ysubj(1,:), Ktope2, flag] = self.MEdE.BubbleT(self.perfil_p(1), mezcla);
                                            [Tr, self.xsubj(1,:), Ktope1, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                                            self.perfil_t(1) = ((self.dflujo)*Tr + (self.salidas{3})*Tb)/(self.dflujo + self.salidas{3});
                                        elseif self.dflujo == 0 || self.dflujo/(self.dflujo + self.salidas{3}) < 1e-4
                                            self.xsubj(1,:) = dsubi./sum(dsubi);
                                            mezcla = Mezcla(self.sust, self.xsubj(1, :));
                                            [Tb, self.ysubj(1,:), Ktope2, flag] = self.MEdE.BubbleT(self.perfil_p(1), mezcla);
                                            mezcla = Mezcla(self.sust, self.ysubj(1, :));
                                            [Tr, ~, Ktope1, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                                            self.perfil_t(1) = ((self.dflujo)*Tr + (self.salidas{3})*Tb)/(self.dflujo + self.salidas{3});
                                        end
                                        Ktope = (self.dflujo .* Ktope2 + self.salidas{3} .* Ktope1)/(self.dflujo + self.salidas{3});
                                    elseif self.salidas{1} == 1 && self.salidas{2} == 0
                                        self.ysubj(1,:) = dsubi./sum(dsubi);
                                        mezcla = Mezcla(self.sust, self.ysubj(1, :));
                                        [self.perfil_t(1), self.xsubj(1,:), Ktope1, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                                        mezcla = Mezcla(self.sust, self.xsubj(1, :));
                                        [Tb, ~, Ktope2, flag] = self.MEdE.BubbleT(self.perfil_p(1), mezcla);
                                        Ktope = Ktope1;
                                    else
                                        self.ysubj(1,:) = dsubi./sum(dsubi);
                                        mezcla = Mezcla(self.sust, self.ysubj(1, :));
                                        [self.perfil_t(1), self.xsubj(1,:), Ktope1, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                                        mezcla = Mezcla(self.sust, self.xsubj(1, :));
                                        [Tb, ~, Ktope2, flag] = self.MEdE.BubbleT(self.perfil_p(1), mezcla);
                                        Ktope = Ktope1;
                                    end
                                else
                                    self.ysubj(1,:) = dsubi./sum(dsubi);
                                    mezcla = Mezcla(self.sust, self.ysubj(1, :));
                                    [self.perfil_t(1), self.xsubj(1,:), Ktope1, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                                    mezcla = Mezcla(self.sust, self.xsubj(1, :));
                                    [Tb, ~, Ktope2, flag] = self.MEdE.BubbleT(self.perfil_p(1), mezcla);
                                    Ktope = Ktope1;
                                end
                            else
                                self.ysubj(1,:) = dsubi./sum(dsubi);
                                mezcla = Mezcla(self.sust, self.ysubj(1, :));
                                [self.perfil_t(1), self.ysubj(1, :), Ktope, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                            end
                        else
                            self.ysubj(1,:) = dsubi./sum(dsubi);
                            mezcla = Mezcla(self.sust, self.ysubj(1, :));
                            [self.perfil_t(1), self.xsubj(1, :), Ktope, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                        end
                    elseif ~any(self.perfil_t == 0)
                        if all(bsubi == 0)
                            bsubi = self.Fi - dsubi;
                        end
                            self.xsubj(end, :) = bsubi./sum(bsubi);
                        if ~isempty(self.salidas) && isa(self.salidas, 'cell')
                            if self.salidas{1} == 1 && self.salidas{2} == 1
                                if ~isempty(self.salidas) && isa(self.salidas, 'cell')
                                    if self.salidas{1} == 1 && self.salidas{2} == 1
                                        if self.dflujo ~= 0 && self.dflujo/(self.dflujo + self.salidas{3}) > 1e-4
                                            zsubj = dsubi./sum(dsubi);
                                            mezcla = Mezcla(self.sust, zsubj);
                                            [Tb, self.ysubj(1,:), Ktope2, flag] = self.MEdE.BubbleT(self.perfil_p(1), mezcla);
                                            [Tr, self.xsubj(1,:), Ktope1, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                                        elseif self.dflujo == 0 || self.dflujo/(self.dflujo + self.salidas{3}) < 1e-4
                                            self.xsubj(1,:) = dsubi./sum(dsubi);
                                            mezcla = Mezcla(self.sust, self.xsubj(1, :));
                                            [Tb, self.ysubj(1,:), Ktope2, flag] = self.MEdE.BubbleT(self.perfil_p(1), mezcla);
                                            mezcla = Mezcla(self.sust, self.ysubj(1, :));
                                            [Tr, ~, Ktope1, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                                        end
                                        Ktope = (self.dflujo .* Ktope2 + self.salidas{3} .* Ktope1)/(self.dflujo + self.salidas{3});
                                    elseif self.salidas{1} == 1 && self.salidas{2} == 0
                                        self.ysubj(1,:) = dsubi./sum(dsubi);
                                        mezcla = Mezcla(self.sust, self.ysubj(1, :));
                                        [~, self.xsubj(1,:), Ktope1, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                                        mezcla = Mezcla(self.sust, self.xsubj(1, :));
                                        [Tb, ~, Ktope2, flag] = self.MEdE.BubbleT(self.perfil_p(1), mezcla);
                                        Ktope = Ktope1;
                                    else
                                        self.ysubj(1,:) = dsubi./sum(dsubi);
                                        mezcla = Mezcla(self.sust, self.ysubj(1, :));
                                        [~, self.xsubj(1,:), Ktope1, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                                        mezcla = Mezcla(self.sust, self.xsubj(1, :));
                                        [Tb, ~, Ktope2, flag] = self.MEdE.BubbleT(self.perfil_p(1), mezcla);
                                        Ktope = Ktope1;
                                    end
                                else
                                    self.ysubj(1,:) = dsubi./sum(dsubi);
                                    mezcla = Mezcla(self.sust, self.ysubj(1, :));
                                    [~, self.xsubj(1,:), Ktope1, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                                    mezcla = Mezcla(self.sust, self.xsubj(1, :));
                                    [Tb, ~, Ktope2, flag] = self.MEdE.BubbleT(self.perfil_p(1), mezcla);
                                    Ktope = Ktope1;
                                end
                            else
                                self.ysubj(1,:) = dsubi./sum(dsubi);
                                mezcla = Mezcla(self.sust, self.ysubj(1, :));
                                [~, self.ysubj(1, :), Ktope, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                            end
                        else
                            self.ysubj(1,:) = dsubi./sum(dsubi);
                            mezcla = Mezcla(self.sust, self.ysubj(1, :));
                            [~, self.xsubj(1, :), Ktope, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                        end
                    end
                    mezcla = Mezcla(self.sust, self.xsubj(end, :));
                    [self.perfil_t(end), self.ysubj(end,:), Kfondo, flag] = self.MEdE.BubbleT(self.perfil_p(end), mezcla);
                    subj = (self.etapas - 1);
                    deltaXsubj = self.xsubj(end,:) - self.xsubj(1,:);
                    deltaTsubj = self.perfil_t(end) - self.perfil_t(1);
                    if all(self.ysubj(1,:)==0)
                        self.ysubj(1,:) = self.xsubj(1,:);
                    end
                    deltaYsubj = self.ysubj(end,:) - self.ysubj(1,:);
                    for i = 2:self.etapas - 1
                        self.xsubj(i, :) = (deltaXsubj)./subj*(i-1) + self.xsubj(1,:);
                        %Normalizo la cantidad xsubj
                        self.xsubj(i,:) = self.xsubj(i,:)./sum(self.xsubj(i,:));
                        self.ysubj(i, :) = (deltaYsubj)./subj*(i-1) + self.ysubj(1,:);
                        %Normalizo la cantidad ysubj
                        self.ysubj(i,:) = self.ysubj(i,:)./sum(self.ysubj(i,:));
                    end
                    for i=2:self.etapas - 1
                        self.perfil_t(i) = (deltaTsubj)./subj.*(i - 1) + self.perfil_t(1);
                    end
                    valI = 0;
                    for i = 1:length(dsubi)
                        if dsubi(i)~=0
                            valI = valI + ((dsubi(i)-self.dsubi(i))/(dsubi(i)))^2;
                        end
                    end
                    for i = 1:length(bsubi)
                        if bsubi(i)~=0
                            valI = valI + ((bsubi(i)-self.bsubi(i))/(bsubi(i)))^2;
                        end
                    end
                    %display(valI);
                    self.dsubi = dsubi;
                    self.bsubi = bsubi;
                self.recover = (self.dsubi(self.lk))/(Flk);
                self.perfil_k = zeros(self.num_sust, self.etapas);
                self.perfil_k(:,1) = Ktope;
                self.perfil_k(:,end) = Kfondo;
                for n = 1: self.num_sust
                    deltaK = -Ktope(n) + Kfondo(n);
                    for i = 2:self.etapas - 1
                        self.perfil_k(n,i) = ((deltaK)./(self.etapas - 1)).*(i-1) + Ktope(n);
                    end
                end
                %if self.dflujo/self.entradas{3} > 1e-3 && isempty(self.perfil_v)
                if isempty(self.perfil_v)
                    self.perfil_v = zeros(1, self.etapas);
                end
                if ~isempty(self.salidas) && isa(self.salidas, 'cell') && all(self.perfil_v == zeros(1, self.etapas))
                    if self.salidas{1} == 1 && self.salidas{2} == 1
                        if abs(self.dflujo./self.salidas{3}) < 1e-3
                            self.perfil_v(1, 2:end) = self.perfil_v(1, 2:end) + ones(1, self.etapas-1).*(self.reflujo + 1).*self.salidas{3};
                        else
                            self.perfil_v(1, 2:end) = self.perfil_v(1, 2:end) + ones(1, self.etapas-1).*(self.reflujo + 1).*(self.salidas{3} + self.dflujo);
                            self.perfil_v(1,1) = self.dflujo;
                        end
                    elseif self.salidas{1} == 1 && self.salidas{2} == 0
                        self.perfil_v(1, 2:end) = self.perfil_v(1, 2:end) +  ones(1, self.etapas-1).*(self.reflujo+1).*(self.dflujo + self.salidas{3});
                        self.perfil_v(1,1) = self.dflujo;
                    else
                        self.perfil_v(1, 2:end) = self.perfil_v(1, 2:end) +  ones(1, self.etapas-1).*(self.reflujo+1).*self.dflujo;
                        self.perfil_v(1,1) = self.dflujo;
                    end
                elseif all(self.perfil_v == zeros(1, self.etapas))
                    self.perfil_v(1, 2:end) = self.perfil_v(1, 2:end) +  ones(1, self.etapas-1).*(self.reflujo+1).*self.dflujo;
                    self.perfil_v(1,1) = self.dflujo;
                end
            elseif any(self.dsubi ~= 0)
                for irrw = 2:2:length(self.alimentaciones)
                    for irre = 1:self.num_sust
                        self.bsubi(irre) = self.alimentaciones{irrw}.molF * self.alimentaciones{irrw}.conc(irre);
                    end
                end
                for irre = 1:self.num_sust
                    self.bsubi(irre) = self.bsubi(irre) - self.dsubi(irre);
                end
                self.xsubj(end, :) = self.bsubi./sum(self.bsubi);
                    self.perfil_t = zeros(1, self.etapas);
                    if ~isempty(self.salidas) && isa(self.salidas, 'cell')
                        if self.salidas{1} == 1 && self.salidas{2} == 1
                            if ~isempty(self.salidas) && isa(self.salidas, 'cell')
                                if self.salidas{1} == 1 && self.salidas{2} == 1
                                    if self.dflujo ~= 0 && self.dflujo/(self.dflujo + self.salidas{3}) > 1e-4
                                        zsubj = self.dsubi./sum(self.dsubi);
                                        mezcla = Mezcla(self.sust, zsubj);
                                        [Tb, self.ysubj(1,:), Ktope2, flag] = self.MEdE.BubbleT(self.perfil_p(1), mezcla);
                                        [Tr, self.xsubj(1,:), Ktope1, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                                        self.perfil_t(1) = ((self.dflujo)*Tr + (self.salidas{3})*Tb)/(self.dflujo + self.salidas{3});
                                    elseif self.dflujo == 0 || self.dflujo/(self.dflujo + self.salidas{3}) < 1e-4
                                        self.xsubj(1,:) = self.dsubi./sum(self.dsubi);
                                        mezcla = Mezcla(self.sust, self.xsubj(1, :));
                                        [Tb, self.ysubj(1,:), Ktope2, flag] = self.MEdE.BubbleT(self.perfil_p(1), mezcla);
                                        mezcla = Mezcla(self.sust, self.ysubj(1, :));
                                        [Tr, ~, Ktope1, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                                        self.perfil_t(1) = ((self.dflujo)*Tr + (self.salidas{3})*Tb)/(self.dflujo + self.salidas{3});
                                    end
                                    Ktope = (self.dflujo .* Ktope2 + self.salidas{3} .* Ktope1)/(self.dflujo + self.salidas{3});
                                elseif self.salidas{1} == 1 && self.salidas{2} == 0
                                    self.ysubj(1,:) = self.dsubi./sum(self.dsubi);
                                    mezcla = Mezcla(self.sust, self.ysubj(1, :));
                                    [self.perfil_t(1), self.xsubj(1,:), Ktope1, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                                    mezcla = Mezcla(self.sust, self.xsubj(1, :));
                                    [Tb, ~, Ktope2, flag] = self.MEdE.BubbleT(self.perfil_p(1), mezcla);
                                    Ktope = Ktope1;
                                else
                                    self.ysubj(1,:) = self.dsubi./sum(self.dsubi);
                                    mezcla = Mezcla(self.sust, self.ysubj(1, :));
                                    [self.perfil_t(1), self.xsubj(1,:), Ktope1, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                                    mezcla = Mezcla(self.sust, self.xsubj(1, :));
                                    [Tb, ~, Ktope2, flag] = self.MEdE.BubbleT(self.perfil_p(1), mezcla);
                                    Ktope = Ktope1;
                                end
                            else
                                self.ysubj(1,:) = self.dsubi./sum(self.dsubi);
                                mezcla = Mezcla(self.sust, self.ysubj(1, :));
                                [self.perfil_t(1), self.xsubj(1,:), Ktope1, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                                mezcla = Mezcla(self.sust, self.xsubj(1, :));
                                [Tb, ~, Ktope2, flag] = self.MEdE.BubbleT(self.perfil_p(1), mezcla);
                                Ktope = Ktope1;
                            end
                        else
                            self.ysubj(1,:) = self.dsubi./sum(self.dsubi);
                            mezcla = Mezcla(self.sust, self.ysubj(1, :));
                            [self.perfil_t(1), self.ysubj(1, :), Ktope, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                        end
                    else
                        self.ysubj(1,:) = self.dsubi./sum(self.dsubi);
                        mezcla = Mezcla(self.sust, self.ysubj(1, :));
                        [self.perfil_t(1), self.xsubj(1, :), Ktope, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                    end
                    mezcla = Mezcla(self.sust, self.xsubj(end, :));
                    [self.perfil_t(end), self.ysubj(end,:), Kfondo, flag] = self.MEdE.BubbleT(self.perfil_p(end), mezcla);
                    subj = (self.etapas - 1);
                    deltaXsubj = self.xsubj(end,:) - self.xsubj(1,:);
                    deltaTsubj = self.perfil_t(end) - self.perfil_t(1);
                    if all(self.ysubj(1,:)==0)
                        self.ysubj(1,:) = self.xsubj(1,:);
                    end
                    deltaYsubj = self.ysubj(end,:) - self.ysubj(1,:);
                    for i = 2:self.etapas - 1
                        self.xsubj(i, :) = (deltaXsubj)./subj*(i-1) + self.xsubj(1,:);
                        %Normalizo la cantidad xsubj
                        self.xsubj(i,:) = self.xsubj(i,:)./sum(self.xsubj(i,:));
                        self.ysubj(i, :) = (deltaYsubj)./subj*(i-1) + self.ysubj(1,:);
                        %Normalizo la cantidad ysubj
                        self.ysubj(i,:) = self.ysubj(i,:)./sum(self.ysubj(i,:));
                    end
                    for i=2:self.etapas - 1
                        self.perfil_t(i) = (deltaTsubj)./subj.*(i - 1) + self.perfil_t(1);
                    end
                    self.perfil_k = zeros(self.num_sust, self.etapas);
                self.perfil_k(:,1) = Ktope;
                self.perfil_k(:,end) = Kfondo;
                for n = 1: self.num_sust
                    deltaK = -Ktope(n) + Kfondo(n);
                    for i = 2:self.etapas - 1
                        self.perfil_k(n,i) = ((deltaK)./(self.etapas - 1)).*(i-1) + Ktope(n);
                    end
                end
                %if self.dflujo/self.entradas{3} > 1e-3 && isempty(self.perfil_v)
                if isempty(self.perfil_v)
                    self.perfil_v = zeros(1, self.etapas);
                end
                if ~isempty(self.salidas) && isa(self.salidas, 'cell') && all(self.perfil_v == zeros(1, self.etapas))
                    if self.salidas{1} == 1 && self.salidas{2} == 1
                        if abs(self.dflujo./self.salidas{3}) < 1e-3
                            self.perfil_v(1, 2:end) = self.perfil_v(1, 2:end) + ones(1, self.etapas-1).*(self.reflujo + 1).*self.salidas{3};
                        else
                            self.perfil_v(1, 2:end) = self.perfil_v(1, 2:end) + ones(1, self.etapas-1).*(self.reflujo + 1).*(self.salidas{3} + self.dflujo);
                            self.perfil_v(1,1) = self.dflujo;
                        end
                    elseif self.salidas{1} == 1 && self.salidas{2} == 0
                        self.perfil_v(1, 2:end) = self.perfil_v(1, 2:end) +  ones(1, self.etapas-1).*(self.reflujo+1).*(self.dflujo + self.salidas{3});
                        self.perfil_v(1,1) = self.dflujo;
                    else
                        self.perfil_v(1, 2:end) = self.perfil_v(1, 2:end) +  ones(1, self.etapas-1).*(self.reflujo+1).*self.dflujo;
                        self.perfil_v(1,1) = self.dflujo;
                    end
                elseif all(self.perfil_v == zeros(1, self.etapas))
                    self.perfil_v(1, 2:end) = self.perfil_v(1, 2:end) +  ones(1, self.etapas-1).*(self.reflujo+1).*self.dflujo;
                    self.perfil_v(1,1) = self.dflujo;
                end
            elseif any(self.bsubi ~= 0)
                for irrw = 2:2:length(self.alimentaciones)
                    for irre = 1:self.num_sust
                        self.dsubi(irre) = self.alimentaciones{irrw}.molF * self.alimentaciones{irrw}.conc(irre);
                    end
                end
                for irre = 1:self.num_sust
                    self.dsubi(irre) = self.dsubi(irre) - self.bsubi(irre);
                end
                self.xsubj(end, :) = self.bsubi./sum(self.bsubi);
                    self.perfil_t = zeros(1, self.etapas);
                    if ~isempty(self.salidas) && isa(self.salidas, 'cell')
                        if self.salidas{1} == 1 && self.salidas{2} == 1
                            if ~isempty(self.salidas) && isa(self.salidas, 'cell')
                                if self.salidas{1} == 1 && self.salidas{2} == 1
                                    if self.dflujo ~= 0 && self.dflujo/(self.dflujo + self.salidas{3}) > 1e-4
                                        zsubj = self.dsubi./sum(self.dsubi);
                                        mezcla = Mezcla(self.sust, zsubj);
                                        [Tb, self.ysubj(1,:), Ktope2, flag] = self.MEdE.BubbleT(self.perfil_p(1), mezcla);
                                        [Tr, self.xsubj(1,:), Ktope1, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                                        self.perfil_t(1) = ((self.dflujo)*Tr + (self.salidas{3})*Tb)/(self.dflujo + self.salidas{3});
                                    elseif self.dflujo == 0 || self.dflujo/(self.dflujo + self.salidas{3}) < 1e-4
                                        self.xsubj(1,:) = self.dsubi./sum(self.dsubi);
                                        mezcla = Mezcla(self.sust, self.xsubj(1, :));
                                        [Tb, self.ysubj(1,:), Ktope2, flag] = self.MEdE.BubbleT(self.perfil_p(1), mezcla);
                                        mezcla = Mezcla(self.sust, self.ysubj(1, :));
                                        [Tr, ~, Ktope1, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                                        self.perfil_t(1) = ((self.dflujo)*Tr + (self.salidas{3})*Tb)/(self.dflujo + self.salidas{3});
                                    end
                                    Ktope = (self.dflujo .* Ktope2 + self.salidas{3} .* Ktope1)/(self.dflujo + self.salidas{3});
                                elseif self.salidas{1} == 1 && self.salidas{2} == 0
                                    self.ysubj(1,:) = self.dsubi./sum(self.dsubi);
                                    mezcla = Mezcla(self.sust, self.ysubj(1, :));
                                    [self.perfil_t(1), self.xsubj(1,:), Ktope1, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                                    mezcla = Mezcla(self.sust, self.xsubj(1, :));
                                    [Tb, ~, Ktope2, flag] = self.MEdE.BubbleT(self.perfil_p(1), mezcla);
                                    Ktope = Ktope1;
                                else
                                    self.ysubj(1,:) = self.dsubi./sum(self.dsubi);
                                    mezcla = Mezcla(self.sust, self.ysubj(1, :));
                                    [self.perfil_t(1), self.xsubj(1,:), Ktope1, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                                    mezcla = Mezcla(self.sust, self.xsubj(1, :));
                                    [Tb, ~, Ktope2, flag] = self.MEdE.BubbleT(self.perfil_p(1), mezcla);
                                    Ktope = Ktope1;
                                end
                            else
                                self.ysubj(1,:) = self.dsubi./sum(self.dsubi);
                                mezcla = Mezcla(self.sust, self.ysubj(1, :));
                                [self.perfil_t(1), self.xsubj(1,:), Ktope1, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                                mezcla = Mezcla(self.sust, self.xsubj(1, :));
                                [Tb, ~, Ktope2, flag] = self.MEdE.BubbleT(self.perfil_p(1), mezcla);
                                Ktope = Ktope1;
                            end
                        else
                            self.ysubj(1,:) = self.dsubi./sum(self.dsubi);
                            mezcla = Mezcla(self.sust, self.ysubj(1, :));
                            [self.perfil_t(1), self.ysubj(1, :), Ktope, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                        end
                    else
                        self.ysubj(1,:) = self.dsubi./sum(self.dsubi);
                        mezcla = Mezcla(self.sust, self.ysubj(1, :));
                        [self.perfil_t(1), self.xsubj(1, :), Ktope, flag] = self.MEdE.DewT(self.perfil_p(1), mezcla);
                    end
                    mezcla = Mezcla(self.sust, self.xsubj(end, :));
                    [self.perfil_t(end), self.ysubj(end,:), Kfondo, flag] = self.MEdE.BubbleT(self.perfil_p(end), mezcla);
                    subj = (self.etapas - 1);
                    deltaXsubj = self.xsubj(end,:) - self.xsubj(1,:);
                    deltaTsubj = self.perfil_t(end) - self.perfil_t(1);
                    if all(self.ysubj(1,:)==0)
                        self.ysubj(1,:) = self.xsubj(1,:);
                    end
                    deltaYsubj = self.ysubj(end,:) - self.ysubj(1,:);
                    for i = 2:self.etapas - 1
                        self.xsubj(i, :) = (deltaXsubj)./subj*(i-1) + self.xsubj(1,:);
                        %Normalizo la cantidad xsubj
                        self.xsubj(i,:) = self.xsubj(i,:)./sum(self.xsubj(i,:));
                        self.ysubj(i, :) = (deltaYsubj)./subj*(i-1) + self.ysubj(1,:);
                        %Normalizo la cantidad ysubj
                        self.ysubj(i,:) = self.ysubj(i,:)./sum(self.ysubj(i,:));
                    end
                    for i=2:self.etapas - 1
                        self.perfil_t(i) = (deltaTsubj)./subj.*(i - 1) + self.perfil_t(1);
                    end
                    self.perfil_k = zeros(self.num_sust, self.etapas);
                self.perfil_k(:,1) = Ktope;
                self.perfil_k(:,end) = Kfondo;
                for n = 1: self.num_sust
                    deltaK = -Ktope(n) + Kfondo(n);
                    for i = 2:self.etapas - 1
                        self.perfil_k(n,i) = ((deltaK)./(self.etapas - 1)).*(i-1) + Ktope(n);
                    end
                end
                %if self.dflujo/self.entradas{3} > 1e-3 && isempty(self.perfil_v)
                if isempty(self.perfil_v)
                    self.perfil_v = zeros(1, self.etapas);
                end
                if ~isempty(self.salidas) && isa(self.salidas, 'cell') && all(self.perfil_v == zeros(1, self.etapas))
                    if self.salidas{1} == 1 && self.salidas{2} == 1
                        if abs(self.dflujo./self.salidas{3}) < 1e-3
                            self.perfil_v(1, 2:end) = self.perfil_v(1, 2:end) + ones(1, self.etapas-1).*(self.reflujo + 1).*self.salidas{3};
                        else
                            self.perfil_v(1, 2:end) = self.perfil_v(1, 2:end) + ones(1, self.etapas-1).*(self.reflujo + 1).*(self.salidas{3} + self.dflujo);
                            self.perfil_v(1,1) = self.dflujo;
                        end
                    elseif self.salidas{1} == 1 && self.salidas{2} == 0
                        self.perfil_v(1, 2:end) = self.perfil_v(1, 2:end) +  ones(1, self.etapas-1).*(self.reflujo+1).*(self.dflujo + self.salidas{3});
                        self.perfil_v(1,1) = self.dflujo;
                    else
                        self.perfil_v(1, 2:end) = self.perfil_v(1, 2:end) +  ones(1, self.etapas-1).*(self.reflujo+1).*self.dflujo;
                        self.perfil_v(1,1) = self.dflujo;
                    end
                elseif all(self.perfil_v == zeros(1, self.etapas))
                    self.perfil_v(1, 2:end) = self.perfil_v(1, 2:end) +  ones(1, self.etapas-1).*(self.reflujo+1).*self.dflujo;
                    self.perfil_v(1,1) = self.dflujo;
                end
            end
            if ~isempty(self.salidas) && isa(self.salidas, 'cell')
                if self.salidas{1} == 1 && self.salidas{2} == 1
                    self.perfil_l(1) = self.perfil_v(2)  - self.dflujo - self.salidas{3};
                elseif self.salidas{2} == 0 && self.salidas{1} == 1
                    self.perfil_l(1) = self.perfil_v(2)  - self.dflujo - self.salidas{3};
                else
                    self.perfil_l(1) = self.perfil_v(2)  - self.dflujo;
                end
            else
                self.perfil_l(1) = self.perfil_v(2)  - self.dflujo;
            end
            self.generar_etapas()
            for iteru = 1:self.etapas
                if iteru == 1
                    self.platos(iteru).salidaV = self.dflujo;
                    if ~isempty(self.salidas)
                        if self.salidas{1} == 1 && self.salidas{2} == 1
                            self.platos(iteru).salidaL = self.salidas{3};
                        end
                    end
                end
                self.platos(iteru).setV(self.perfil_v(iteru));
                self.platos(iteru).setyi(self.ysubj(iteru,:));
                self.platos(iteru).K = self.perfil_k(:, iteru);
                if iteru == self.etapas
                    self.platos(iteru).salidaL = self.bflujo;
                end
            end
            if ~isempty(self.salidas) && isa(self.salidas, 'cell')
                if self.salidas{1} == 1 && self.salidas{2} == 1
                    self.dflujo = self.dflujo + self.salidas{3};
                elseif self.salidas{1} == 1 && self.salidas{2} == 0
                    self.dflujo = self.dflujo + self.salidas{3};
                end
            end

            for itert = 2:self.etapas
                if itert >= 2
                    self.perfil_l(itert) = -self.perfil_v(1);
                end
                if itert < self.etapas
                    self.perfil_l(itert) = self.perfil_l(itert) + self.perfil_v(itert + 1);
                end
                for iters = 1:itert
                    if ~isempty(self.platos(iters).aliment)
                        if self.platos(iters).aliment.q >= 1
                            self.perfil_l(itert) = self.perfil_l(itert) + self.platos(iters).aliment.molF;
                        elseif self.platos(iters).aliment.q < 1 && self.platos(iters).aliment.q > 0
                            self.perfil_l(itert) = self.perfil_l(itert) + (self.platos(iters).aliment.q).*self.platos(iters).aliment.molF;
                        end
                    end
                end
                for iters = 1:itert
                    if iters ~= 1
                        if ~isempty(self.platos(iters).salidaV)
                            self.perfil_l(itert) = self.perfil_l(itert) - self.platos(iters).salidaV;
                        end
                    end
                    if iters ~= self.etapas
                        if ~isempty(self.platos(iters).salidaL)
                            self.perfil_l(itert) = self.perfil_l(itert) - self.platos(iters).salidaL;
                        end
                    end
                end
            end
            self.perfil_l(end) = self.bflujo;

            for op = 2:2:length(self.alimentaciones)
                if self.alimentaciones{op}.q <= 1 && self.alimentaciones{op}.q  >= 0
                    self.perfil_v(self.alimentaciones{op-1}:end) = self.perfil_v(self.alimentaciones{op-1}:end) - self.alimentaciones{op}.molF*(1-self.alimentaciones{op}.q);
                    if any(self.perfil_v < 0)
                        warning('La condicion de reflujo espeficificada es muy peque�a, pruebe una mayor')
                    end
                elseif self.alimentaciones{op}.q < 0
                    mezclaX = Mezcla(self.sust, self.xsubj(self.alimentaciones{op-1},:), self.alimentaciones{2}.mezcla.kij);
                    mezclaY = Mezcla(self.sust, self.ysubj(self.alimentaciones{op-1},:), self.alimentaciones{2}.mezcla.kij);
                    Hvap = self.MEdE.entalpia(self.perfil_t(self.alimentaciones{op-1}), self.perfil_p(self.alimentaciones{op-1}), mezclaX, 'liq') - self.MEdE.entalpia(self.perfil_t(self.alimentaciones{op-1}), self.perfil_p(self.alimentaciones{op-1}), mezclaY, 'liq');
                    self.perfil_l(self.alimentaciones{op-1}:end-1) =  self.perfil_l(self.alimentaciones{op-1}:end-1) + ( self.alimentaciones{op}.q*self.alimentaciones{op}.molF*Hvap) / Hvap;
                    self.perfil_v(self.alimentaciones{op-1}:end) = self.perfil_v(self.alimentaciones{op-1}:end) - self.alimentaciones{op}.molF + ( self.alimentaciones{op}.q*self.alimentaciones{op}.molF*Hvap) / Hvap;

                else
                    mezclaX = Mezcla(self.sust, self.xsubj(self.alimentaciones{op-1},:), self.alimentaciones{2}.mezcla.kij);
                    mezclaY = Mezcla(self.sust, self.ysubj(self.alimentaciones{op-1},:), self.alimentaciones{2}.mezcla.kij);
                    Hvap = self.MEdE.entalpia(self.perfil_t(self.alimentaciones{op-1}), self.perfil_p(self.alimentaciones{op-1}), mezclaX, 'liq') - self.MEdE.entalpia(self.perfil_t(self.alimentaciones{op-1}), self.perfil_p(self.alimentaciones{op-1}), mezclaY, 'liq');
                    self.perfil_l(self.alimentaciones{op-1}:end-1) =  self.perfil_l(self.alimentaciones{op-1}:end-1) + ( (self.alimentaciones{op}.q - 1)*self.alimentaciones{op}.molF*Hvap) / Hvap;
                    self.perfil_v(self.alimentaciones{op-1}:end) = self.perfil_v(self.alimentaciones{op-1}:end) + ( (self.alimentaciones{op}.q - 1)*self.alimentaciones{op}.molF*Hvap) / Hvap;
                end
            end
            self.perfil_l(end) = self.bflujo;
            if isempty(self.perfil_vi)
                for o = 1:self.etapas
                    self.perfil_vi(o,:) = self.ysubj(o,:).*self.perfil_v(o);
                end
            end
            if isempty(self.perfil_li)
                for o = 1:self.etapas
                    self.perfil_li(o,:) = self.xsubj(o,:).*self.perfil_l(o);
                end
            end
        end
        function self = wanghenke(self)
            if isempty(self.perfil_k) || isempty(self.perfil_v) || isempty(self.perfil_t)
                self = self.balanmasa();
            end
            iter = 0;
            if isempty(self.respaldovalI)
                valI = uint32(2^32);
            else
                valI = self.respaldovalI;
            end
            for iteru=1:self.etapas
                self.platos(iteru).setL(self.perfil_l(iteru));
                self.platos(iteru).setxi(self.xsubj(iteru,:));
            end
            if isempty(self.tol)
                self.tol = 1e-8.*self.etapas;
            end
            try
                respaldoV = self.perfil_v;
                respaldoL = self.perfil_l;
                respaldovi = self.perfil_vi;
                respaldoli = self.perfil_li;
                respaldot = self.perfil_t;
                respaldovalI = valI;
                respaldok = self.perfil_k;
                respaldoqc = self.qc;
                respaldoqb = self.qb;
                respaldoplatos = self.platos;
                respaldoxi = self.xsubj;
                respaldoyi = self.ysubj;
            if isempty(self.respaldoV)
                self.respaldoV = zeros(2, self.etapas);
                self.respaldoV(1, :) = self.perfil_v;
            else
                self.respaldoV((self.actualiter+1)+1, :) = zeros(1, self.etapas);
            end
            if isempty(self.respaldoL)
                self.respaldoL = zeros(2, self.etapas);
                self.respaldoL = self.perfil_l;
            else
                self.respaldoL((self.actualiter+1)+1, :) = zeros(1, self.etapas);
            end
            if isempty(self.respaldovi)
                self.respaldovi = zeros(2.*self.etapas, self.num_sust);
                self.respaldovi(1:self.etapas, :) = self.perfil_vi;
            else
                self.respaldovi((self.actualiter+1).*self.etapas+1:(self.actualiter+1).*self.etapas + self.etapas, :) = zeros(self.etapas, self.num_sust);
            end
            if isempty(self.respaldoli)
                self.respaldoli = zeros(2.*self.etapas, self.num_sust);
                self.respaldoli(1:self.etapas, :) = self.perfil_li;
            else
                self.respaldoli((self.actualiter+1).*self.etapas+1:(self.actualiter+1).*self.etapas + self.etapas, :) = zeros(self.etapas, self.num_sust);
            end
            if isempty(self.respaldot)
                self.respaldot = zeros(2, self.etapas);
                self.respaldot(1, :) = self.perfil_t;
            else
                self.respaldot((self.actualiter+1)+1, :) = zeros(1, self.etapas);
            end
            if isempty(self.respaldoqc)
                self.respaldoqc = zeros(2,1);
            else
                self.respaldoqc((self.actualiter+1)+1,1) = 0;
            end
            if isempty(self.respaldoqc)
                self.respaldoqb = zeros(2,1);
            else
                self.respaldoqb((self.actualiter+1)+1) = 0;
            end
            while valI > self.tol && iter < self.iteraciones
                iter = iter + 1;
                %%%%%% Voy componente a componente resolviendo Thomas
                tamano = size(self.entradas);
                for ierto =1: self.num_sust
                    Dij = zeros(1, self.etapas);
                    Cij = zeros(1, self.etapas - 1);
                    Bij = zeros(1, self.etapas);
                    Aij = zeros(1, self.etapas - 1);
                    Cij(1,:) = Cij(1, :) + self.perfil_v(2:end).*(self.perfil_k(ierto,2:end));
                    Bij(:) = self.perfil_v(1);
                    Aij(:) = -self.perfil_v(1);
                    for k = 1:self.etapas
                        for j = 1:3:tamano(2)
                            if self.entradas{j} == k
                                Dij(k) = -self.platos(k).aliment.molF.*self.concfeed((j+2)/3, ierto);
                            end
                        end
                        if k < self.etapas
                            Bij(k) = Bij(k) - (self.platos(k+1).V);
                        end
                        if k > 1
                            Aij(k-1) = Aij(k-1) + (self.platos(k).V);
                        end
                        for h = 1:k
                            if ~isempty(self.platos(h).aliment)
                                Bij(k) = Bij(k) - self.platos(h).aliment.molF;
                            end
                        end
                        if k > 1
                            for h = 1:k-1
                                if ~isempty(self.platos(h).aliment)
                                    Aij(k-1) = Aij(k-1) + self.platos(h).aliment.molF;
                                end
                            end
                        end
                        for h = 1:k
                            if ~isempty(self.platos(h).salidaV)
                                if h ~= 1
                                    Bij(k) = Bij(k) + self.platos(h).salidaV;
                                end
                            end
                            if ~isempty(self.platos(h).salidaL)
                                Bij(k) = Bij(k) + self.platos(h).salidaL;
                            end
                        end
                        if k > 1
                            for h = 1:k-1
                                if ~isempty(self.platos(h).salidaL)
                                    if h~= self.etapas
                                        Aij(k-1) = Aij(k-1) - self.platos(h).salidaL;
                                    end
                                end
                                if ~isempty(self.platos(h).salidaV)
                                    if h~= 1
                                        Aij(k-1) = Aij(k-1) - self.platos(h).salidaV;
                                    end
                                end
                            end
                        end
                        if ~isempty(self.platos(k).salidaV)
                            if k > 1
                                Bij(k) = Bij(k) - (self.perfil_v(k) + self.platos(k).salidaV)*(self.perfil_k(ierto,k)) ;
                            else
                                Bij(k) = Bij(k) - (self.perfil_v(k))*(self.perfil_k(ierto,k)) ;
                            end
                        else
                            Bij(k) = Bij(k) - (self.perfil_v(k))*(self.perfil_k(ierto,k)) ;
                        end
                        if ~isempty(self.platos(k).salidaL)
                            Bij(k) = Bij(k) - self.platos(k).salidaL;
                        end
                    end
                    self.xsubj(:,ierto) = tridiagThomas(Bij, Aij, Cij, Dij);
                    self.xsubj(:,ierto) = abs(self.xsubj(:,ierto));
                end
                for i = 1:self.etapas
                    self.xsubj(i,:) = self.xsubj(i,:)./(sum(self.xsubj(i,:))); %Normalizo
                end
                Tbackup = self.perfil_t;
                Vbackup = self.perfil_v;
                for i = 1:self.etapas
                    mezcla = Mezcla(self.sust, self.xsubj(i, :));
                    mezcla = mezcla.definek(self.alimentaciones{2}.mezcla.kij);
                    [self.perfil_t(i), self.ysubj(i,:), self.perfil_k(:,i), flag] = self.MEdE.BubbleT(self.perfil_p(i), mezcla);
                end
                self.bsubi = self.xsubj(end,:).*(self.perfil_l(end));
                if ~isempty(self.salidas) && isa(self.salidas, 'cell')
                    if self.salidas{1} == 1 && self.salidas{2} == 1
                        self.dsubi = (self.dflujo - self.salidas{3}).* self.ysubj(1,:) + (self.salidas{3}).*self.xsubj(1,:);
                    else
                        self.dsubi = (self.dflujo).* self.ysubj(1,:);
                    end
                else
                    self.dsubi = (self.dflujo).* self.ysubj(1,:);
                end
                %V2 is obtained from (10-6) Seader & Henley & Roper:
                %Separation Process Principles 3rd Edition
                %V2 se obtiene de (10-6) Seader & Henley & Roper

                %Condenser Duty
                %Calor del condensador
                Fhf = 0;
                if self.alimentaciones{1} == 1 % Primer plato
                    Fhf = Fhf + self.alimentaciones{2}.molF*self.alimentaciones{2}.H;
                end
                Uh = 0; %Salidas Laterales Liquidas 1
                Wh = 0; %Salidas Laterales Vapor 1
                Lh = 0; %Liquido que va al plato 2
                Vh = 0; %Pudiera haber vapor fuga del plato 1
                Vh2 = 0; %Vapor que entra al plato 1 del plato 2
                cpi_gi = cell.empty(0,self.num_sust);
                Href = zeros(1, self.num_sust);


                for i = 1:self.num_sust
                    Href(i) = self.sust(i).href;
                    cpi_gi{i} = self.sust(i).cp_gi{1};
                end
                mezclaL1 = Mezcla(self.sust, self.xsubj(1, :), self.alimentaciones{2}.mezcla.kij);
                HL1gi = 0;
                HV1gi = 0;
                HV2gi = 0;
                deltaHL1 = zeros(1, self.num_sust);
                for i = 1:self.num_sust
                    deltaHL1(i) = integral(cpi_gi{i}, 273.15, self.perfil_t(1));
                    HL1gi = HL1gi  + deltaHL1(i) * self.xsubj(1,i) + Href(i) * self.xsubj(1,i);
                end
                HdepL1 = self.MEdE.entalpia(self.platos(1).T,  self.perfil_p(1), mezclaL1, 'liq');
                Hdepref = self.MEdE.entalpia(273.15, 101.325, mezclaL1, 'liq');
                HL1 = HL1gi - HdepL1 + Hdepref;
                %HL1 = self.MEdE.entalpia(TL1, PL1, mezcla, 'liq');

                if ~isempty(self.platos(1).salidaL)
                    Uh = Uh + self.platos(1).salidaL;
                end
                Lh = Lh + self.perfil_l(1);
                Vh = Vh + self.perfil_v(1);
                mezclaV1 = Mezcla(self.sust, self.ysubj(1, :), self.alimentaciones{2}.mezcla.kij);
                for i = 1:self.num_sust
                    HV1gi = HV1gi  + deltaHL1(i) * self.ysubj(1, i) + Href(i) * self.ysubj(1, i);
                end
                HdepV1 = self.MEdE.entalpia(self.platos(1).T,  self.perfil_p(1), mezclaV1, 'vap');
                Hdepref = self.MEdE.entalpia(273.15, 101.325, mezclaV1, 'liq');
                HV1 = HV1gi - HdepV1 + Hdepref;
                if ~isempty(self.platos(1).salidaV)
                    Wh = Wh + self.platos(1).salidaV;
                end

                concV2 = self.ysubj(2, :);
                flujoV2 = self.perfil_v(2);
                TV2 = self.perfil_t(2);
                mezclaV2 = Mezcla(self.sust, concV2, self.alimentaciones{2}.mezcla.kij);
                deltaHV2 = zeros(1 , self.num_sust);
                for i = 1:self.num_sust
                    deltaHV2(i) = integral(cpi_gi{i}, 273.15, TV2);
                    HV2gi = HV2gi  + deltaHV2(i) *  self.ysubj(2, i) + Href(i) *  self.ysubj(2, i);
                end
                HdepV2 = self.MEdE.entalpia(self.platos(2).T, self.perfil_p(2), mezclaV2, 'vap' );
                Hdepref = self.MEdE.entalpia(273.15, 101.325, mezclaV2, 'liq');
                HV2 = HV2gi - HdepV2 + Hdepref;
                Vh2 = Vh2 + flujoV2;
                self.qc = Vh2*HV2 + Fhf - (Lh + Uh)*HL1 - (Vh + (Wh - Vh))*HV1;
                Fh = 0;
                Uh = 0;
                Vh = 0;
                Vh1 = 0;
                LhEND = 0;
                Wh = 0;
                Lh = 0;
                for iterx = 1:self.etapas
                    if iterx == 1 && self.alimentaciones{1} == 1
                        Fh = Fh + self.alimentaciones{2}.molF * self.alimentaciones{2}.H;
                    elseif iterx == self.alimentaciones{1}
                        Fh = Fh + self.alimentaciones{2}.molF * self.alimentaciones{2}.H;
                    elseif length(self.alimentaciones)>2
                        if iterx == self.alimentaciones{3}
                            Fh = Fh + self.alimentaciones{4}.molF * self.alimentaciones{4}.H;
                        end
                    elseif length(self.alimentaciones) > 4
                        if iterx == self.alimentaciones{5}
                            Fh = Fh + self.alimentaciones{6}.molF * self.alimentaciones{6}.H;
                        end
                    elseif length(self.alimentaciones ) > 6
                        if iterx == self.alimentaciones{7}
                            Fh = Fh + self.alimentaciones{8}.molF * self.alimentaciones{8}.H;
                        end
                    elseif length(self.alimentaciones) > 8
                        if iterx == self.alimentaciones{9}
                            Fh = Fh + self.alimentaciones{10}.molF * self.alimentaciones{10}.H;
                        end
                    elseif length(self.alimentaciones ) > 10
                        if iterx == self.alimentaciones{11}
                            Fh = Fh + self.alimentaciones{12}.molF * self.alimentaciones{12}.H;
                        end
                    elseif length(self.alimentaciones) > 12
                        if iterx == self.alimentaciones{13}
                            Fh = Fh + self.alimentaciones{14}.molF * self.alimentaciones{14}.H;
                        end
                    end
                    deltaHgi = zeros(1 , self.num_sust);
                    HLgi = 0;
                    HVgi = 0;

                    if ~isempty(self.platos(iterx).salidaV)
                        concVk = self.ysubj(iterx,:);
                        mezclaVk = Mezcla(self.sust, concVk, self.alimentaciones{2}.mezcla.kij);
                        PVk = self.perfil_p(iterx);
                        TVk = self.perfil_t(iterx);

                        for i = 1:self.num_sust
                            deltaHgi(i) = integral(cpi_gi{i}, 273.15, TVk);
                            HVgi = HVgi  + deltaHgi(i) * concVk(i) + Href(i) * concVk(i);
                        end
                        Hdep_Vk_ref = self.MEdE.entalpia(273.15, 101.325, mezclaVk, 'liq');
                        Hdep_Vk = self.MEdE.entalpia(TVk, PVk, mezclaVk, 'vap');
                        HVk = HVgi - Hdep_Vk + Hdep_Vk_ref;
                        if iterx == 1
                            Vh1 = Vh1 +  self.perfil_v(1) * HVk;  % que no se sume
                                %2 veces el destilado vapor
                        else
                            Wh = Wh + self.platos(iterx).salidaV * HVk;
                        end
                    end
                    if ~isempty(self.platos(iterx).salidaL)
                        concLk = self.xsubj(iterx,:);
                        mezclaLk = Mezcla(self.sust, concLk,self.alimentaciones{2}.mezcla.kij);
                        PLk = self.perfil_p(iterx);
                        TLk = self.perfil_t(iterx);

                        for i = 1:self.num_sust
                            if deltaHgi(i) == 0
                                deltaHgi(i) = integral(cpi_gi{i}, 273.15, TLk);
                            end
                            HLgi = HLgi  + deltaHgi(i) * self.xsubj(iterx,i) + Href(i) * self.xsubj(iterx,i) ;
                        end
                        Hdep_Lk = self.MEdE.entalpia(TLk, PLk, mezclaLk, 'liq');
                        Hdep_Lk_ref = self.MEdE.entalpia(273.15, 101.325, mezclaLk, 'liq');
                        HLk = HLgi - Hdep_Lk + Hdep_Lk_ref;
                        if iterx == self.etapas
                            LhEND = LhEND  + self.platos(iterx).salidaL * HLk;
                            %que no se sume dos veces el fondo l�quido
                        else
                            Uh = Uh + self.platos(iterx).salidaL * HLk;
                        end
                    end
                end
                self.qb = Fh - Uh  - Wh - Vh1 - LhEND - sum(self.perfil_q) - self.qc;

                %Matriz didiagonal de platos 3 a N de donde se despeja Vj

                Dij = zeros(1, self.etapas - 2);
                Cij = zeros(1, self.etapas - 3);
                Bij = zeros(1, self.etapas - 2);
                Aij = zeros(1, self.etapas - 3);
                alfa2 = HL1 - HV2;
                Dij = Dij - self.perfil_v(1);
                deltaHgi = zeros(1, self.num_sust);

                for iterz = 3:self.etapas
                    NL = iterz - 1;
                    PL = self.perfil_p(NL);
                    TL = self.perfil_t(NL);
                    mezclaNL = Mezcla(self.sust, self.xsubj(NL, :));
                    mezclaNL.definek(self.alimentaciones{2}.mezcla.kij);
                    HLgijm1 = 0;
                    for i = 1:self.num_sust
                        deltaHgi(i) = integral(cpi_gi{i}, 273.15, TL);
                        HLgijm1 = HLgijm1  + deltaHgi(i) * self.xsubj(NL, i) + Href(i);
                    end
                    HLdep_jm1 = self.MEdE.entalpia(TL, PL, mezclaNL, 'liq');
                    Hdep_Ljm1_ref = self.MEdE.entalpia(273.15, 101.325, mezclaNL, 'liq');
                    HLjm1 = HLgijm1 - HLdep_jm1 + Hdep_Ljm1_ref;
                    NV = iterz;
                    PV = self.perfil_p(NV);
                    TV = self.perfil_t(NV);

                    mezclaNV = Mezcla(self.sust, self.ysubj(NV, :));
                    mezclaNV.definek(self.alimentaciones{2}.mezcla.kij);
                    HVgi = 0;
                    for i = 1:self.num_sust
                        deltaHgi(i) = integral(cpi_gi{i}, 273.15, TV);
                        HVgi = HVgi  + deltaHgi(i) * self.ysubj(NV, i) + Href(i);
                    end
                    HVdep_j = self.MEdE.entalpia(TV, PV, mezclaNV, 'vap');
                    Hdep_Vj_ref = self.MEdE.entalpia(273.15, 101.325, mezclaNV, 'liq');

                    HVj =  HVgi - HVdep_j + Hdep_Vj_ref;

                    mezclaNVm1 = Mezcla(self.sust, self.ysubj(NL, :),self.alimentaciones{2}.mezcla.kij);
                    HVgi = 0;
                    for i = 1:self.num_sust
                        deltaHgi(i) = integral(cpi_gi{i}, 273.15, TL);
                        HVgi = HVgi  + deltaHgi(i) * self.ysubj(NL, i) + Href(i) * self.ysubj(NL, i);
                    end
                    HVdep_jm1 = self.MEdE.entalpia(TL, PL, mezclaNVm1, 'vap');
                    Hdep_Vjm1_ref = self.MEdE.entalpia(273.15, 101.325, mezclaNVm1, 'liq');
                    HVjm1 = HVgi - HVdep_jm1 + Hdep_Vjm1_ref;

                    if iterz < self.etapas
                        Aij(NL-1) = (HLjm1 - HVj);     %alpha j
                    end
                    for iteraa = 1 : NL-1
                        if ~isempty(self.platos(iteraa).aliment)
                            Dij(NL-1) = Dij(NL-1) +  self.platos(iteraa).aliment.molF;
                        end
                        if ~isempty(self.platos(iteraa).salidaV)
                            Dij(NL-1) = Dij(NL-1) -  self.platos(iteraa).salidaV;
                        end
                        if ~isempty(self.platos(iteraa).salidaL)
                            Dij(NL-1) = Dij(NL-1) -  self.platos(iteraa).salidaL;
                        end
                    end
                    %corriente = Corriente(mezclaNL,  TV, 'T', 0, 'x', self.perfil_l(NV),'m' , self.MEdE);

                    %HLj = self.MEdE.entalpia(TV, PV, mezcla, 'liq');
                    Dij(NL-1) = Dij(NL-1)*(HLjm1 - HL1);
                    if ~isempty(self.platos(NL).aliment)
                        Dij(NL-1) = Dij(NL-1) + self.platos(NL).aliment.molF * (HLjm1 - self.platos(NL).aliment.H);
                    end
                    Dij(NL-1) = Dij(NL-1) + self.perfil_q(NL-1);

                    if ~isempty(self.platos(iterz).salidaV)
                        Dij(NL-1) = Dij(NL-1) +  self.platos(NL).salidaV*(HVjm1 - HLjm1);
                    end
                    if iterz == 3
                        Dij(NL-1) = Dij(NL-1) - alfa2 * self.perfil_v(2);
                    end
                    Bij(NL-1) = HVj - HLjm1;
                end
                self.perfil_v(3:end) = tridiagThomas(Bij, Aij, Cij, Dij);
                self.perfil_l(2:end) = -self.perfil_v(1);
                if ~isempty(self.platos(1).salidaL)
                    self.perfil_l(2:end) = self.perfil_l(2:end) - self.platos(1).salidaL;
                end
                for itert = 2:self.etapas
                    if itert < self.etapas
                        self.perfil_l(itert) = self.perfil_l(itert) + self.perfil_v(itert + 1);
                    end
                    for iters = 2:itert
                        if ~isempty(self.platos(iters).aliment)
                            self.perfil_l(itert) = self.perfil_l(itert) + self.platos(iters).aliment.molF;
                        end
                        if ~isempty(self.platos(iters).salidaV)
                            if iters ~= 1
                                self.perfil_l(itert) = self.perfil_l(itert) - self.platos(iters).salidaV;
                            end
                        end
                        if ~isempty(self.platos(iters).salidaL) && itert ~= self.etapas
                            self.perfil_l(itert) = self.perfil_l(itert) - self.platos(iters).salidaL;
                        end
                    end
                end
                valI = 0;
                for i = 2:length(Tbackup)
                    valI = valI + ((self.perfil_t(i)-Tbackup(i))/self.perfil_t(i))^2;
                end
                for i = 2:length(Tbackup)
                    valI = valI + ((self.perfil_v(i)-Vbackup(i))/self.perfil_v(i))^2;
                end
                for iteru = 1:self.etapas
                    self.platos(iteru).setV(self.perfil_v(iteru));
                    self.platos(iteru).setL(self.perfil_l(iteru));
                    self.platos(iteru).setyi(self.ysubj(iteru,:));
                    self.platos(iteru).setxi(self.xsubj(iteru,:));
                    self.platos(iteru).K = self.perfil_k(:, iteru);
                end
                self.bsubi = self.xsubj(end,:).*(self.perfil_l(end));
                if ~isempty(self.salidas) && isa(self.salidas, 'cell')
                    if self.salidas{1} == 1 && self.salidas{2} == 1
                        self.dsubi = (self.dflujo - self.salidas{3}).* self.ysubj(1,:) + (self.salidas{3}).*self.xsubj(1,:);
                    else
                        self.dsubi = (self.dflujo).* self.ysubj(1,:);
                    end
                else
                    self.dsubi = (self.dflujo).* self.ysubj(1,:);
                end
                display(valI);
                tamano = size(self.ysubj);
                self.perfil_vi = zeros(tamano(1), tamano(2));
                self.perfil_li = self.perfil_vi;
                for iter1t = 1:self.etapas
                    self.perfil_vi(iter1t,:) = self.perfil_v(iter1t).*(self.ysubj(iter1t,:));
                    self.perfil_li(iter1t,:) = self.perfil_l(iter1t).*(self.xsubj(iter1t,:));
                end
                for iiter0 = 1:self.etapas
                    self.platos(iiter0).v_i = self.perfil_vi(iiter0,:);
                    self.platos(iiter0).l_i = self.perfil_li(iiter0,:);
                    self.platos(iiter0).y_i = self.ysubj(iiter0, :);
                    self.platos(iiter0).x_i = self.xsubj(iiter0, :);
                    self.platos(iiter0).V = self.perfil_v(iiter0);
                    self.platos(iiter0).L = self.perfil_l(iiter0);
                    self.platos(iiter0).T = self.perfil_t(iiter0);
                    self.platos(iiter0).K = self.perfil_k(:,iiter0);
                end
                if ( valI < respaldovalI)
                    respaldovalI = (valI);
                    respaldoV = self.perfil_v;
                    respaldoL = self.perfil_l;
                    respaldovi = self.perfil_vi;
                    respaldoli = self.perfil_li;
                    respaldot = self.perfil_t;
                    respaldok =  self.perfil_k;
                    respaldoplatos = self.platos;
                    respaldoxi = self.xsubj;
                    respaldoyi = self.ysubj;
                    respaldoqc= self.qc;
                    respaldoqb = self.qb;
                    self.respaldoiter = iter;
                end

                self.actualiter = iter;
                self.actualvalI = valI;
                self.ysubj(self.ysubj<0) = 1e-9;
                self.xsubj(self.xsubj<0) = 1e-9;
                self.perfil_vi(self.perfil_vi<0) = 1e-6;
                self.perfil_li(self.perfil_li<0) = 1e-6;
                self.respaldoV(self.actualiter + 1, :) = self.perfil_v;
                self.respaldoL(self.actualiter + 1, :) = self.perfil_l;
                self.respaldovi(self.actualiter.*self.etapas+1:self.actualiter.*self.etapas + self.etapas, :) = self.perfil_vi;
                self.respaldoli(self.actualiter.*self.etapas+1:self.actualiter.*self.etapas + self.etapas, :) = self.perfil_li;
                self.respaldot(self.actualiter+1, :) = self.perfil_t;
                self.respaldoqc(self.actualiter+1) = self.qc;
                self.respaldoqb(self.actualiter+1) = self.qb;

                if valI < self.tol
                    self.laststep = logical(1);
                end
            end
            ajustefinal = (self.bsubi + self.dsubi) - self.Fi;
            if any(ajustefinal ~= 0)
                for i = 1:self.num_sust
                    ajustesubi = self.bsubi(i) + self.dsubi(i);
                    self.bsubi(i) = self.Fi(i).*self.bsubi(i)./ajustesubi;
                    self.dsubi(i) = self.Fi(i).*self.dsubi(i)./ajustesubi;
                end
            end
            catch
                self.perfil_v = respaldoV;
                self.perfil_l = respaldoL;
                self.perfil_vi = respaldovi;
                self.perfil_li = respaldoli;
                self.perfil_t = respaldot;
                self.perfil_k = respaldok;
                self.qc = respaldoqc;
                self.qb = respaldoqb;
                self.platos = respaldoplatos;
                self.xsubj = respaldoxi;
                self.ysubj = respaldoyi;
                display(self);
                fprintf('Fallo la convergencia en la iteracion = %i \n', iter);
                fprintf('El mejor resultado obtenido fue en la iteracion = %i \n', self.respaldoiter)
                fprintf('El mejor error obtenido fue de = %f \n', respaldovalI)
            end
        end
        function avanza1paso(self)
            if isempty(self.perfil_k) || isempty(self.perfil_l) || isempty(self.perfil_v) || isempty(self.perfil_t)
                self.balanmasa();
            end
            if isempty(self.tol)
                self.tol = 1e-8.*self.etapas;
            end
            if isempty(self.respaldovalI)
                valI = 1e7;
            end
            if ~isempty(self.actualiter)
                iter = self.actualiter;
            else
                iter = 0;
            end
            try
                respaldoV = self.perfil_v;
                respaldoL = self.perfil_l;
                respaldovi = self.perfil_vi;
                respaldoli = self.perfil_li;
                respaldot = self.perfil_t;
                respaldovalI = valI;
                respaldok = self.perfil_k;
                respaldoqc = self.qc;
                respaldoqb = self.qb;
                respaldoplatos = self.platos;
                respaldoxi = self.xsubj;
                respaldoyi = self.ysubj;
            if isempty(self.respaldoV)
                self.respaldoV = zeros(2, self.etapas);
                self.respaldoV(1, :) = self.perfil_v;
            else
                self.respaldoV((self.actualiter+1)+1, :) = zeros(1, self.etapas);
            end
            if isempty(self.respaldoL)
                self.respaldoL = zeros(2, self.etapas);
                self.respaldoL = self.perfil_l;
            else
                self.respaldoL((self.actualiter+1)+1, :) = zeros(1, self.etapas);
            end
            if isempty(self.respaldovi)
                self.respaldovi = zeros(2.*self.etapas, self.num_sust);
                self.respaldovi(1:self.etapas, :) = self.perfil_vi;
            else
                self.respaldovi((self.actualiter+1).*self.etapas+1:(self.actualiter+1).*self.etapas + self.etapas, :) = zeros(self.etapas, self.num_sust);
            end
            if isempty(self.respaldoli)
                self.respaldoli = zeros(2.*self.etapas, self.num_sust);
                self.respaldoli(1:self.etapas, :) = self.perfil_li;
            else
                self.respaldoli((self.actualiter+1).*self.etapas+1:(self.actualiter+1).*self.etapas + self.etapas, :) = zeros(self.etapas, self.num_sust);
            end
            if isempty(self.respaldot)
                self.respaldot = zeros(2, self.etapas);
                self.respaldot(1, :) = self.perfil_t;
            else
                self.respaldot((self.actualiter+1)+1, :) = zeros(1, self.etapas);
            end
            if isempty(self.respaldoqc)
                self.respaldoqc = zeros(2,1);
            else
                self.respaldoqc((self.actualiter+1)+1,1) = 0;
            end
            if isempty(self.respaldoqc)
                self.respaldoqb = zeros(2,1);
            else
                self.respaldoqb((self.actualiter+1)+1) = 0;
            end
                iter = iter + 1;
                %%%%%% Voy componente a componente resolviendo Thomas
                tamano = size(self.entradas);
                for ierto =1: self.num_sust
                    Dij = zeros(1, self.etapas);
                    Cij = zeros(1, self.etapas - 1);
                    Bij = zeros(1, self.etapas);
                    Aij = zeros(1, self.etapas - 1);
                    Cij(1,:) = Cij(1, :) + self.perfil_v(2:end).*(self.perfil_k(ierto,2:end));
                    Bij(:) = self.perfil_v(1);
                    Aij(:) = -self.perfil_v(1);
                    for k = 1:self.etapas
                        for j = 1:3:tamano(2)
                            if self.entradas{j} == k
                                Dij(k) = -self.platos(k).aliment.molF.*self.concfeed((j+2)/3, ierto);
                            end
                        end
                        if k < self.etapas
                            Bij(k) = Bij(k) - (self.platos(k+1).V);
                        end
                        if k > 1
                            Aij(k-1) = Aij(k-1) + (self.platos(k).V);
                        end
                        for h = 1:k
                            if ~isempty(self.platos(h).aliment)
                                Bij(k) = Bij(k) - self.platos(h).aliment.molF;
                            end
                        end
                        if k > 1
                            for h = 1:k-1
                                if ~isempty(self.platos(h).aliment)
                                    Aij(k-1) = Aij(k-1) + self.platos(h).aliment.molF;
                                end
                            end
                        end
                        for h = 1:k
                            if ~isempty(self.platos(h).salidaV)
                                if h ~= 1
                                    Bij(k) = Bij(k) + self.platos(h).salidaV;
                                end
                            end
                            if ~isempty(self.platos(h).salidaL)
                                Bij(k) = Bij(k) + self.platos(h).salidaL;
                            end
                        end
                        if k > 1
                            for h = 1:k-1
                                if ~isempty(self.platos(h).salidaL)
                                    if h~= self.etapas
                                        Aij(k-1) = Aij(k-1) - self.platos(h).salidaL;
                                    end
                                end
                                if ~isempty(self.platos(h).salidaV)
                                    if h~= 1
                                        Aij(k-1) = Aij(k-1) - self.platos(h).salidaV;
                                    end
                                end
                            end
                        end
                        if ~isempty(self.platos(k).salidaV)
                            if k > 1
                                Bij(k) = Bij(k) - (self.perfil_v(k) + self.platos(k).salidaV)*(self.perfil_k(ierto,k)) ;
                            else
                                Bij(k) = Bij(k) - (self.perfil_v(k))*(self.perfil_k(ierto,k)) ;
                            end
                        else
                            Bij(k) = Bij(k) - (self.perfil_v(k))*(self.perfil_k(ierto,k)) ;
                        end
                        if ~isempty(self.platos(k).salidaL)
                            Bij(k) = Bij(k) - self.platos(k).salidaL;
                        end
                    end
                    self.xsubj(:,ierto) = tridiagThomas(Bij, Aij, Cij, Dij);
                    self.xsubj(:,ierto) = abs(self.xsubj(:,ierto));
                end
                for i = 1:self.etapas
                    self.xsubj(i,:) = self.xsubj(i,:)./(sum(self.xsubj(i,:))); %Normalizo
                end
                Tbackup = self.perfil_t;
                Vbackup = self.perfil_v;
                for i = 1:self.etapas
                    mezcla = Mezcla(self.sust, self.xsubj(i, :));
                    mezcla = mezcla.definek(self.alimentaciones{2}.mezcla.kij);
                    [self.perfil_t(i), self.ysubj(i,:), self.perfil_k(:,i), flag] = self.MEdE.BubbleT(self.perfil_p(i), mezcla);
                end
                self.bsubi = self.xsubj(end,:).*(self.perfil_l(end));
                if ~isempty(self.salidas) && isa(self.salidas, 'cell')
                    if self.salidas{1} == 1 && self.salidas{2} == 1
                        self.dsubi = (self.dflujo - self.salidas{3}).* self.ysubj(1,:) + (self.salidas{3}).*self.xsubj(1,:);
                    else
                        self.dsubi = (self.dflujo).* self.ysubj(1,:);
                    end
                else
                    self.dsubi = (self.dflujo).* self.ysubj(1,:);
                end
                %V2 is obtained from (10-6) Seader & Henley & Roper:
                %Separation Process Principles 3rd Edition
                %V2 se obtiene de (10-6) Seader & Henley & Roper

                %Condenser Duty
                %Calor del condensador
                Fhf = 0;
                if self.alimentaciones{1} == 1 % Primer plato
                    Fhf = Fhf + self.alimentaciones{2}.molF*self.alimentaciones{2}.H;
                end
                Uh = 0; %Salidas Laterales Liquidas 1
                Wh = 0; %Salidas Laterales Vapor 1
                Lh = 0; %Liquido que va al plato 2
                Vh = 0; %Pudiera haber vapor fuga del plato 1
                Vh2 = 0; %Vapor que entra al plato 1 del plato 2
                cpi_gi = cell.empty(0,self.num_sust);
                Href = zeros(1, self.num_sust);


                for i = 1:self.num_sust
                    Href(i) = self.sust(i).href;
                    cpi_gi{i} = self.sust(i).cp_gi{1};
                end
                mezclaL1 = Mezcla(self.sust, self.xsubj(1, :), self.alimentaciones{2}.mezcla.kij);
                HL1gi = 0;
                HV1gi = 0;
                HV2gi = 0;
                deltaHL1 = zeros(1, self.num_sust);
                for i = 1:self.num_sust
                    deltaHL1(i) = integral(cpi_gi{i}, 273.15, self.perfil_t(1));
                    HL1gi = HL1gi  + deltaHL1(i) * self.xsubj(1,i) + Href(i) * self.xsubj(1,i);
                end
                HdepL1 = self.MEdE.entalpia(self.platos(1).T,  self.perfil_p(1), mezclaL1, 'liq');
                Hdepref = self.MEdE.entalpia(273.15, 101.325, mezclaL1, 'liq');
                HL1 = HL1gi - HdepL1 + Hdepref;
                %HL1 = self.MEdE.entalpia(TL1, PL1, mezcla, 'liq');

                if ~isempty(self.platos(1).salidaL)
                    Uh = Uh + self.platos(1).salidaL;
                end
                Lh = Lh + self.perfil_l(1);
                Vh = Vh + self.perfil_v(1);
                mezclaV1 = Mezcla(self.sust, self.ysubj(1, :), self.alimentaciones{2}.mezcla.kij);
                for i = 1:self.num_sust
                    HV1gi = HV1gi  + deltaHL1(i) * self.ysubj(1, i) + Href(i) * self.ysubj(1, i);
                end
                HdepV1 = self.MEdE.entalpia(self.platos(1).T,  self.perfil_p(1), mezclaV1, 'vap');
                Hdepref = self.MEdE.entalpia(273.15, 101.325, mezclaV1, 'liq');
                HV1 = HV1gi - HdepV1 + Hdepref;
                if ~isempty(self.platos(1).salidaV)
                    Wh = Wh + self.platos(1).salidaV;
                end

                concV2 = self.ysubj(2, :);
                flujoV2 = self.perfil_v(2);
                TV2 = self.perfil_t(2);
                mezclaV2 = Mezcla(self.sust, concV2, self.alimentaciones{2}.mezcla.kij);
                deltaHV2 = zeros(1 , self.num_sust);
                for i = 1:self.num_sust
                    deltaHV2(i) = integral(cpi_gi{i}, 273.15, TV2);
                    HV2gi = HV2gi  + deltaHV2(i) *  self.ysubj(2, i) + Href(i) *  self.ysubj(2, i);
                end
                HdepV2 = self.MEdE.entalpia(self.platos(2).T, self.perfil_p(2), mezclaV2, 'vap' );
                Hdepref = self.MEdE.entalpia(273.15, 101.325, mezclaV2, 'liq');
                HV2 = HV2gi - HdepV2 + Hdepref;
                Vh2 = Vh2 + flujoV2;
                self.qc = Vh2*HV2 + Fhf - (Lh + Uh)*HL1 - (Vh + (Wh - Vh))*HV1;
                Fh = 0;
                Uh = 0;
                Vh = 0;
                Vh1 = 0;
                LhEND = 0;
                Wh = 0;
                Lh = 0;
                for iterx = 1:self.etapas
                    if iterx == 1 && self.alimentaciones{1} == 1
                        Fh = Fh + self.alimentaciones{2}.molF * self.alimentaciones{2}.H;
                    elseif iterx == self.alimentaciones{1}
                        Fh = Fh + self.alimentaciones{2}.molF * self.alimentaciones{2}.H;
                    elseif length(self.alimentaciones)>2
                        if iterx == self.alimentaciones{3}
                            Fh = Fh + self.alimentaciones{4}.molF * self.alimentaciones{4}.H;
                        end
                    elseif length(self.alimentaciones) > 4
                        if iterx == self.alimentaciones{5}
                            Fh = Fh + self.alimentaciones{6}.molF * self.alimentaciones{6}.H;
                        end
                    elseif length(self.alimentaciones ) > 6
                        if iterx == self.alimentaciones{7}
                            Fh = Fh + self.alimentaciones{8}.molF * self.alimentaciones{8}.H;
                        end
                    elseif length(self.alimentaciones) > 8
                        if iterx == self.alimentaciones{9}
                            Fh = Fh + self.alimentaciones{10}.molF * self.alimentaciones{10}.H;
                        end
                    elseif length(self.alimentaciones ) > 10
                        if iterx == self.alimentaciones{11}
                            Fh = Fh + self.alimentaciones{12}.molF * self.alimentaciones{12}.H;
                        end
                    elseif length(self.alimentaciones) > 12
                        if iterx == self.alimentaciones{13}
                            Fh = Fh + self.alimentaciones{14}.molF * self.alimentaciones{14}.H;
                        end
                    end
                    deltaHgi = zeros(1 , self.num_sust);
                    HLgi = 0;
                    HVgi = 0;

                    if ~isempty(self.platos(iterx).salidaV)
                        concVk = self.ysubj(iterx,:);
                        mezclaVk = Mezcla(self.sust, concVk, self.alimentaciones{2}.mezcla.kij);
                        PVk = self.perfil_p(iterx);
                        TVk = self.perfil_t(iterx);

                        for i = 1:self.num_sust
                            deltaHgi(i) = integral(cpi_gi{i}, 273.15, TVk);
                            HVgi = HVgi  + deltaHgi(i) * concVk(i) + Href(i) * concVk(i);
                        end
                        Hdep_Vk_ref = self.MEdE.entalpia(273.15, 101.325, mezclaVk, 'liq');
                        Hdep_Vk = self.MEdE.entalpia(TVk, PVk, mezclaVk, 'vap');
                        HVk = HVgi - Hdep_Vk + Hdep_Vk_ref;
                        if iterx == 1
                            Vh1 = Vh1 +  self.perfil_v(1) * HVk;  % que no se sume
                                %2 veces el destilado vapor
                        else
                            Wh = Wh + self.platos(iterx).salidaV * HVk;
                        end
                    end
                    if ~isempty(self.platos(iterx).salidaL)
                        concLk = self.xsubj(iterx,:);
                        mezclaLk = Mezcla(self.sust, concLk,self.alimentaciones{2}.mezcla.kij);
                        PLk = self.perfil_p(iterx);
                        TLk = self.perfil_t(iterx);

                        for i = 1:self.num_sust
                            if deltaHgi(i) == 0
                                deltaHgi(i) = integral(cpi_gi{i}, 273.15, TLk);
                            end
                            HLgi = HLgi  + deltaHgi(i) * self.xsubj(iterx,i) + Href(i) * self.xsubj(iterx,i) ;
                        end
                        Hdep_Lk = self.MEdE.entalpia(TLk, PLk, mezclaLk, 'liq');
                        Hdep_Lk_ref = self.MEdE.entalpia(273.15, 101.325, mezclaLk, 'liq');
                        HLk = HLgi - Hdep_Lk + Hdep_Lk_ref;
                        if iterx == self.etapas
                            LhEND = LhEND  + self.platos(iterx).salidaL * HLk;
                            %que no se sume dos veces el fondo l�quido
                        else
                            Uh = Uh + self.platos(iterx).salidaL * HLk;
                        end
                    end
                end
                self.qb = Fh - Uh  - Wh - Vh1 - LhEND - sum(self.perfil_q) - self.qc;

                %Matriz didiagonal de platos 3 a N de donde se despeja Vj

                Dij = zeros(1, self.etapas - 2);
                Cij = zeros(1, self.etapas - 3);
                Bij = zeros(1, self.etapas - 2);
                Aij = zeros(1, self.etapas - 3);
                alfa2 = HL1 - HV2;
                Dij = Dij - self.perfil_v(1);
                deltaHgi = zeros(1, self.num_sust);

                for iterz = 3:self.etapas
                    NL = iterz - 1;
                    PL = self.perfil_p(NL);
                    TL = self.perfil_t(NL);
                    mezclaNL = Mezcla(self.sust, self.xsubj(NL, :));
                    mezclaNL.definek(self.alimentaciones{2}.mezcla.kij);
                    HLgijm1 = 0;
                    for i = 1:self.num_sust
                        deltaHgi(i) = integral(cpi_gi{i}, 273.15, TL);
                        HLgijm1 = HLgijm1  + deltaHgi(i) * self.xsubj(NL, i) + Href(i);
                    end
                    HLdep_jm1 = self.MEdE.entalpia(TL, PL, mezclaNL, 'liq');
                    Hdep_Ljm1_ref = self.MEdE.entalpia(273.15, 101.325, mezclaNL, 'liq');
                    HLjm1 = HLgijm1 - HLdep_jm1 + Hdep_Ljm1_ref;
                    NV = iterz;
                    PV = self.perfil_p(NV);
                    TV = self.perfil_t(NV);

                    mezclaNV = Mezcla(self.sust, self.ysubj(NV, :));
                    mezclaNV.definek(self.alimentaciones{2}.mezcla.kij);
                    HVgi = 0;
                    for i = 1:self.num_sust
                        deltaHgi(i) = integral(cpi_gi{i}, 273.15, TV);
                        HVgi = HVgi  + deltaHgi(i) * self.ysubj(NV, i) + Href(i);
                    end
                    HVdep_j = self.MEdE.entalpia(TV, PV, mezclaNV, 'vap');
                    Hdep_Vj_ref = self.MEdE.entalpia(273.15, 101.325, mezclaNV, 'liq');

                    HVj =  HVgi - HVdep_j + Hdep_Vj_ref;

                    mezclaNVm1 = Mezcla(self.sust, self.ysubj(NL, :),self.alimentaciones{2}.mezcla.kij);
                    HVgi = 0;
                    for i = 1:self.num_sust
                        deltaHgi(i) = integral(cpi_gi{i}, 273.15, TL);
                        HVgi = HVgi  + deltaHgi(i) * self.ysubj(NL, i) + Href(i) * self.ysubj(NL, i);
                    end
                    HVdep_jm1 = self.MEdE.entalpia(TL, PL, mezclaNVm1, 'vap');
                    Hdep_Vjm1_ref = self.MEdE.entalpia(273.15, 101.325, mezclaNVm1, 'liq');
                    HVjm1 = HVgi - HVdep_jm1 + Hdep_Vjm1_ref;

                    if iterz < self.etapas
                        Aij(NL-1) = (HLjm1 - HVj);     %alpha j
                    end
                    for iteraa = 1 : NL-1
                        if ~isempty(self.platos(iteraa).aliment)
                            Dij(NL-1) = Dij(NL-1) +  self.platos(iteraa).aliment.molF;
                        end
                        if ~isempty(self.platos(iteraa).salidaV)
                            Dij(NL-1) = Dij(NL-1) -  self.platos(iteraa).salidaV;
                        end
                        if ~isempty(self.platos(iteraa).salidaL)
                            Dij(NL-1) = Dij(NL-1) -  self.platos(iteraa).salidaL;
                        end
                    end
                    %corriente = Corriente(mezclaNL,  TV, 'T', 0, 'x', self.perfil_l(NV),'m' , self.MEdE);

                    %HLj = self.MEdE.entalpia(TV, PV, mezcla, 'liq');
                    Dij(NL-1) = Dij(NL-1)*(HLjm1 - HL1);
                    if ~isempty(self.platos(NL).aliment)
                        Dij(NL-1) = Dij(NL-1) + self.platos(NL).aliment.molF * (HLjm1 - self.platos(NL).aliment.H);
                    end
                    Dij(NL-1) = Dij(NL-1) + self.perfil_q(NL-1);

                    if ~isempty(self.platos(iterz).salidaV)
                        Dij(NL-1) = Dij(NL-1) +  self.platos(NL).salidaV*(HVjm1 - HLjm1);
                    end
                    if iterz == 3
                        Dij(NL-1) = Dij(NL-1) - alfa2 * self.perfil_v(2);
                    end
                    Bij(NL-1) = HVj - HLjm1;
                end
                self.perfil_v(3:end) = tridiagThomas(Bij, Aij, Cij, Dij);
                self.perfil_l(2:end) = -self.perfil_v(1);
                if ~isempty(self.platos(1).salidaL)
                    self.perfil_l(2:end) = self.perfil_l(2:end) - self.platos(1).salidaL;
                end
                for itert = 2:self.etapas
                    if itert < self.etapas
                        self.perfil_l(itert) = self.perfil_l(itert) + self.perfil_v(itert + 1);
                    end
                    for iters = 2:itert
                        if ~isempty(self.platos(iters).aliment)
                            self.perfil_l(itert) = self.perfil_l(itert) + self.platos(iters).aliment.molF;
                        end
                        if ~isempty(self.platos(iters).salidaV)
                            if iters ~= 1
                                self.perfil_l(itert) = self.perfil_l(itert) - self.platos(iters).salidaV;
                            end
                        end
                        if ~isempty(self.platos(iters).salidaL) && itert ~= self.etapas
                            self.perfil_l(itert) = self.perfil_l(itert) - self.platos(iters).salidaL;
                        end
                    end
                end
                valI = 0;
                for i = 2:length(Tbackup)
                    valI = valI + ((self.perfil_t(i)-Tbackup(i))/self.perfil_t(i))^2;
                end
                for i = 2:length(Tbackup)
                    valI = valI + ((self.perfil_v(i)-Vbackup(i))/self.perfil_v(i))^2;
                end
                for iteru = 1:self.etapas
                    self.platos(iteru).setV(self.perfil_v(iteru));
                    self.platos(iteru).setL(self.perfil_l(iteru));
                    self.platos(iteru).setyi(self.ysubj(iteru,:));
                    self.platos(iteru).setxi(self.xsubj(iteru,:));
                    self.platos(iteru).K = self.perfil_k(:, iteru);
                end
                self.bsubi = self.xsubj(end,:).*(self.perfil_l(end));
                if ~isempty(self.salidas) && isa(self.salidas, 'cell')
                    if self.salidas{1} == 1 && self.salidas{2} == 1
                        self.dsubi = (self.dflujo - self.salidas{3}).* self.ysubj(1,:) + (self.salidas{3}).*self.xsubj(1,:);
                    else
                        self.dsubi = (self.dflujo).* self.ysubj(1,:);
                    end
                else
                    self.dsubi = (self.dflujo).* self.ysubj(1,:);
                end
                display(valI);
                tamano = size(self.ysubj);
                self.perfil_vi = zeros(tamano(1), tamano(2));
                self.perfil_li = self.perfil_vi;
                for iter1t = 1:self.etapas
                    self.perfil_vi(iter1t,:) = self.perfil_v(iter1t).*(self.ysubj(iter1t,:));
                    self.perfil_li(iter1t,:) = self.perfil_l(iter1t).*(self.xsubj(iter1t,:));
                end
                for iiter0 = 1:self.etapas
                    self.platos(iiter0).v_i = self.perfil_vi(iiter0,:);
                    self.platos(iiter0).l_i = self.perfil_li(iiter0,:);
                    self.platos(iiter0).y_i = self.ysubj(iiter0, :);
                    self.platos(iiter0).x_i = self.xsubj(iiter0, :);
                    self.platos(iiter0).V = self.perfil_v(iiter0);
                    self.platos(iiter0).L = self.perfil_l(iiter0);
                    self.platos(iiter0).T = self.perfil_t(iiter0);
                    self.platos(iiter0).K = self.perfil_k(:,iiter0);
                end
                if ( valI < respaldovalI)
                    respaldovalI = (valI);
                    respaldoV = self.perfil_v;
                    respaldoL = self.perfil_l;
                    respaldovi = self.perfil_vi;
                    respaldoli = self.perfil_li;
                    respaldot = self.perfil_t;
                    respaldok =  self.perfil_k;
                    respaldoplatos = self.platos;
                    respaldoxi = self.xsubj;
                    respaldoyi = self.ysubj;
                    respaldoqc= self.qc;
                    respaldoqb = self.qb;
                    self.respaldoiter = iter;
                end
                self.actualiter = iter;
                self.actualvalI = valI;
                self.respaldoV(self.actualiter + 1, :) = self.perfil_v;
                self.respaldoL(self.actualiter + 1, :) = self.perfil_l;
                self.respaldovi(self.actualiter.*self.etapas+1:self.actualiter.*self.etapas + self.etapas, :) = self.perfil_vi;
                self.respaldoli(self.actualiter.*self.etapas+1:self.actualiter.*self.etapas + self.etapas, :) = self.perfil_li;
                self.respaldot(self.actualiter+1, :) = self.perfil_t;
                self.respaldoqc(self.actualiter+1) = self.qc;
                self.respaldoqb(self.actualiter+1) = self.qb;

                if valI < self.tol
                    self.laststep = logical(1);
                end
                ajustefinal = (self.bsubi + self.dsubi) - self.Fi;
                if any(ajustefinal ~= 0)
                    for i = 1:self.num_sust
                        ajustesubi = self.bsubi(i) + self.dsubi(i);
                        self.bsubi(i) = self.Fi(i).*self.bsubi(i)./ajustesubi;
                        self.dsubi(i) = self.Fi(i).*self.dsubi(i)./ajustesubi;
                    end
                end
            catch
                self.perfil_v = respaldoV;
                self.perfil_l = respaldoL;
                self.perfil_vi = respaldovi;
                self.perfil_li = respaldoli;
                self.perfil_t = respaldot;
                self.perfil_k = respaldok;
                self.qc = respaldoqc;
                self.qb = respaldoqb;
                self.platos = respaldoplatos;
                self.xsubj = respaldoxi;
                self.ysubj = respaldoyi;
                display(self);
                fprintf('Fallo la convergencia en la iteracion = %i \n', iter);
                fprintf('El mejor resultado obtenido fue en la iteracion = %i \n', self.respaldoiter)
                fprintf('El mejor error obtenido fue de = %f \n', respaldovalI)
            end

        end

        function generar_etapas(self, perfil_v, perfil_k)
            self.platos = Plato.empty(0, self.etapas);
            if nargin > 1 && ~isempty(perfil_v)
                self.perfil_v = perfil_v;
            end
            if nargin > 2 && ~isempty(perfil_k)
                self.perfil_k = perfil_k;
            end
            etapaentrada = zeros(1, floor(length(self.entradas)/3));
            for i = 1:length(etapaentrada)
                etapaentrada(i) = self.alimentaciones{i*2-1};
            end
            aliment = Corriente.empty(0,length(self.alimentaciones)/2);
            for i = 2:2:length(self.alimentaciones)
                aliment(floor(i/2+1)) = self.alimentaciones{i};
            end
            etapasalida = zeros(1, floor(length(self.salidas)/3));
            for i = 1:3:length(self.salidas)
                etapasalida(i) = self.salidas{i};
            end
            salida = zeros(1, floor(length(self.salidas))/3);
            for i = 3:3:length(self.salidas)
                salida(floor(i/3)) =  self.salidas{i};
            end
            for itei = 1:self.etapas
                plato_aliment = find(etapaentrada == itei, 1, 'first');
                plato_salida = find(etapasalida == itei, 1, 'first');
                if itei == 1
                    if isempty(plato_salida)
                        if self.dflujo > 1e-3
                            plato_salida = 1;
                        end
                    end
                end

                if isempty(plato_aliment) && isempty(plato_salida)
                    self.platos(1, itei) = Plato(itei, self.perfil_t(itei), self.perfil_p(itei), [], [], self.xsubj(itei,:), self.ysubj(itei,:), self.perfil_v(itei), [], self.perfil_k(:, itei), [], [], [], [], []);
                elseif ~isempty(plato_aliment) && isempty(plato_salida)
                    self.platos(1, itei) = Plato(itei, self.perfil_t(itei), self.perfil_p(itei), [], [], self.xsubj(itei,:), self.ysubj(itei,:), self.perfil_v(itei), [], self.perfil_k(:, itei),  self.alimentaciones{plato_aliment*2}, [], [], [], []);
                elseif ~isempty(plato_salida) && isempty(plato_aliment)
                    if itei == 1
                        salidasy = self.dflujo;
                        self.platos(1, itei) = Plato(itei, self.perfil_t(itei), self.perfil_p(itei), [], [], self.xsubj(itei,:), self.ysubj(itei,:), self.perfil_v(itei),[],self.perfil_k(:, itei),  [], salidasy,[], [], []);
                    elseif abs(self.salidas{plato_salida*3-1}-1) <= 1e-4 %Si es liquido
                        salidasx = self.salidas{plato_salida*3};
                        salidasy = self.dflujo;
                        self.platos(1, itei) = Plato(itei, self.perfil_t(itei), self.perfil_p(itei), [], [], self.xsubj(itei,:), self.ysubj(itei,:), self.perfil_v(itei), [], self.perfil_k(:, itei), [], salidasy, salidasx, [], []);
                    elseif abs(self.salidas(plato_salida*3-1)) <= 1e-4 %Si es vapor
                        salidasy = self.salidas{plato_salida*3};
                        self.platos(1, itei) = Plato(itei, self.perfil_t(itei), self.perfil_p(itei), [], [], self.xsubj(itei,:), self.ysubj(itei,:), self.perfil_v(itei),[],self.perfil_k(:, itei),  [], salidasy,[], [], []);
                    end
                end
            end
            if self.dflujo > 1e-3
                self.perfil_v(1) = self.dflujo;
            end
            if self.platos(1).v_i == zeros(1, self.num_sust)
                self.platos(1).v_i = self.platos(1).l_i;
            end
        end
        function recalculo_etapas(self)
            etapaentrada = zeros(1, floor(length(self.entradas)/2));
            for i = 1:3:length(self.entradas)
                etapaentrada(floor(i/2+1)) = self.entradas{i};
            end
            aliment = Corriente.empty(0,length(self.alimentaciones)/2);
            for i = 2:2:length(self.alimentaciones)
                aliment(floor(i/2+1)) = self.alimentaciones{i};
            end
            etapasalida = zeros(1, floor(length(self.salidas)/3));
            for i = 1:3:length(self.salidas)
                etapasalida(i) = self.salidas{i};
            end
            salida = zeros(1, floor(length(self.salidas))/3);
            for i = 3:3:length(self.salidas)
                salida(floor(i/3)) =  self.salidas{i};
            end
            self.platos = Plato.empty(0,self.etapas);
            for itei = 1:self.etapas
                plato_aliment = find(etapaentrada == itei, 1, 'first');
                plato_salida = find(etapasalida == itei, 1, 'first');
                if isempty(plato_aliment) && isempty(plato_salida)
                    self.platos(1, itei) = Plato(itei, self.perfil_t(itei), self.perfil_p(itei), [], [], self.xsubj(itei,:), self.ysubj(itei,:), self.perfil_v(itei), self.perfil_l(itei), self.perfil_k(:, itei), [], [], [], [], []);
                    self.platos(1, itei).T = self.perfil_t(itei);
                    self.platos(1,itei).P = self.perfil_p(itei);
                    self.platos(1, itei) = self.platos(1, itei).setxi(self.xsubj(itei,:));
                    self.platos(1, itei) = self.platos(1, itei).setyi(self.ysubj(itei,:));
                end
                if ~isempty(plato_aliment) && isempty(plato_salida)
                    self.platos(1, itei) = Plato(itei, self.perfil_t(itei), self.perfil_p(itei), [], [], self.xsubj(itei,:), self.ysubj(itei,:), self.perfil_v(itei), self.perfil_l(itei), self.perfil_k(:, itei),  self.alimentaciones{plato_aliment*2}, [], [], [], []);
                end
                if ~isempty(plato_salida) && isempty(plato_aliment)
                    if abs(self.salidas{plato_salida*3-1}-1) <= 1e-4 %Si es l�quido
                        salidasx = self.salidas{plato_salida*3};
                        salidasy = self.dflujo - self.salidas{plato_salida*3};
                        self.platos(1, itei) = Plato(itei, self.perfil_t(itei), self.perfil_p(itei), [], [], self.xsubj(:, itei), self.ysubj(:,itei), self.perfil_v(itei), self.perfil_l(itei), self.perfil_k(:, itei), [],[] , salidasx,[], []);
                    elseif abs(self.salidas(plato_salida*3-1)) <= 1e-4  %Si es vapor
                        salidasy = self.salidas{plato_salida*3};
                        self.platos(1, itei) = Plato(itei, self.perfil_t(itei), self.perfil_p(itei), [], [], self.xsubj(:, itei), self.ysubj(:,itei), self.perfil_v(itei), self.perfil_l(itei),self.perfil_k(:, itei),  [], salidasy,[], [], []);
                    end
                end
            end
            if self.platos(1).v_i == zeros(1, self.num_sust)
                self.platos(1).v_i = self.platos(1).l_i;
            end
        end
    end
end
