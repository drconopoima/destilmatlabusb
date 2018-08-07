classdef NaphtaliSandholm < handle
    %NaphtaliSandholm Método riguroso de Naphtali Sandholm para resolución de
    %   columnas de destilación complejas.
    %   
    %   Referencias: 
    %       - Seader, Henley, Roper. "Separation Process Principles". 3rd
    %       Edition. 
    %   Luis Jesús Díaz Manzo
    %
    %
    %
%Separaci�n de una mezcla binaria por m�todos rigurosos
%
%Una mezcla binaria de Butano, Pentano, 11% molar y con un par�metro de interacci�n binaria de 0.0011
%
%mezcla1 = Mezcla ([Sustancia('Butane'), Sustancia('Pentane')], [0.11,0.89], [1.1e-3]);
%
%A 100 PSIA, en l�quido saturado
%
%corriente1 = Corriente(mezcla1, 0, 'x', 689.475728, 'P', 100, 'm', RMVdW(PREdE));
%   
%   
%El m�todo de NEWTON RAPHSON:
%
%El m�todo de Newton Raphson, que toma en cuenta los balances de energ�a. 
%
%En este m�todo, el usuario podr�a introducir tambi�n un conjunto de variables similar al del m�todo de Bubble Point
%lo cual iterar�a y algunos casos converge, sin embargo, es extremadamente sensible a los cambios, es en ocasiones 
%necesario realizar ajustes bajando el damping factor al rango de 0.2-0.6, aunque es posible que esto sea consecuencia
%de que mi versi�n de Matlab sea de 32 bits.
%
%Una variante del m�todo es que en este se permite la introducci�n de estimados iniciales para cualquiera de los flujos
%de destilado 'di' , residuo 'bi' , perfiles de vapor 'v' o l�quido 'l' en la columna, perfil de coeficientes de reparto K 'k',
%perfiles de vapor o de l�quido por componentes 'vi' o 'li'.
%
%Con el m�nimo n�mero de par�metros la torre se puede pasar a resolver de la siguiente forma:

%Torre3 = NaphtaliSandholm(8, {4, corriente1}, Corriente.empty(0,1), Corriente.empty(0,1), {1,1,40}, 689.475728, 0, [],  [], [], [], [], 0.6, 60 , 'B', 3, 'R');
% NaphtaliSandholm(N�m.Platos, Corriente Alimentaci�n,  Corriente destilado y residuo, Salida etapa 1, l�quida q=1, 40 moles, Presi�n, adiab�tica q =0, ,,,,0.6 de damping factor, Especificaciones: 60 flujo molar residuo  y 3 coeficiente de reflujo (L/D))
% 
% Los resultados son:
% 
% Comienza el m�todo con converger los balances de masa los que estima a ser: Destilado moles: [11, 29] 11moles tope Butano, 21 Pentano, Moles fondo: [0, 60] todo el fondo es Pentano
% 
% Valores convergidos 
% 
%Flujos molares de destilado
%Torre2.dsubi = 10.4725, 29.5275
%Flujos molares de residuo de fondo
%Torre2.bsubi = 0.5275,  59.4725
%
%Perfil de temperaturas:
%
%Torre3.perfil_t = 356.9981, 371.6742, 374.7421, 376.0672, 377.4402, 378.3983, 379.0472, 379.4780
%
%Perfil de vapor:
%
%Torre3.perfil_v = 0.5172 , 160 , 169.3314, 169.9981, 170.8733,  171.6058, 172.1624, 172.5605

%Perfil de l�quido:

%Torre3.perfil_l = 120, 129.3314, 129.9981, 230.8733, 231.6058, 232.1624, 232.5605, 60
%
%Perfiles por componentes individuales: Torre2.perfil_vi y Torre2.perfil_li

%La debilidad del m�todo es que no converge los balances de energ�a, �nicamente cumple los balances de masa.
%

%Se puede hacer la columna tambi�n utilizando destilado vapor, con ambos m�todos, se escribir�a igual aunque reemplazando el t�rmino {1,1,40} de las salidas.
%
%Torre3 = NaphtaliSandholm(8, {4, corriente1}, Corriente.empty(0,1), Corriente.empty(0,1), {}, 689.475728, 0, 3,  [], [], [], [], 0.8, 40 , 'D', 3, 'R'); (alternativamente tambi�n doy equivalentemente flujo molar de D destilado sin problema, resuelve internamente los balances)
%                                                                                          {} = Sin salida = Salida vapor        el damping factor puede elevarse en este caso a 0.8 y converge
%
%Ese caso convergen a estos perfiles de temperaturas: Torre3.perfil_t = [371.6378, 374.665, 375.953, 376.480, 377.6947, 378.5475, 379.1302, 379.5209] 
%                                                                                          
%                                                                                          
%Para dar los estimados iniciales y facilitar la convergencia basta con escribirlos al final de la redacci�n
%
%Torre3 = NaphtaliSandholm(8, {4, corriente1}, Corriente.empty(0,1), Corriente.empty(0,1), {}, 689.475728, 0, 3,  [], [], [], [], 0.7, 40 , 'D', 3, 'R', [0.5,59.5], 'bi', [365.2696, 371.7052, 374.7180, 375.9715, 377.3161, 378.2996, 378.9864, 379.4471], 'tj');
%                                       Despu�s del damping factor van 2 especificaciones, despu�s pueden ir arbitrariamente todos los estimados iniciales que se usen, por ejemplo en este caso bi de [0.5, 59.5] o perfilesTj de 365.2696, 371.7052, 374.7180, 375.9715, 377.3161, 378.2996, 378.9864, 379.4471
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
        qc
        qb
        dflujo
        bflujo
        reflujo
        iteraciones
        num_sust
        alimentaciones
        cdest
        cfond
        salidas = cell.empty(0,1)
        entradas
        tol
        damping
        gdlibertad
        comp
        sust
        MEdE
        dsubi
        bsubi
        ysubj
        xsubj
        varX
        concfeed
        matriz_tridiag
        lk
        hk
        recover
        Fi
        UWi
        platos
        convergence
        platos_entradas
        Aij
        Bij
        Cij
        Fk
        F
        UmW = 0;
        respaldoiter = 0;
        laststep = logical(0);
        diferencia
        respaldoV
        respaldoL
        respaldovi
        respaldoli
        respaldot
        respaldoqc
        respaldoqb
        actualiter
        actualvalI
        respaldovalI
    end
    methods
        function self = NaphtaliSandholm(platos, alim, dest, bot, salid, Pperfil, Qperfil, num_iter, dobflow, d_o_b, perfil_vi, perfil_li, damping, cond1, whoscond1, cond2, whoscond2, varargin)
            if nargin > 12 && ~isempty(damping)
                if isa(damping, 'double') && length(damping) == 1
                    self.damping = damping;
                else
                    self.damping = 1;
                end
            end
            if nargin > 0 && ~isempty(platos) 
                self.etapas = platos;
            end     
            if nargin > 1 && ~isempty(alim) && isa(alim, 'cell') && isa(alim{2}, 'Corriente')
                self.sust = alim{2}.mezcla.comp;
            end
            self.comp = cell.empty(0,0);
            
            if nargin > 1 && ~isempty(alim) && isa(alim, 'cell') && isa(alim{2}, 'Corriente')
                self.definealim(platos, alim)
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
            if nargin > 5 && isempty(Pperfil) && ~isempty(platos)
            %Si no se especifica presion, se usará el de la primera alimentación en toda la columna
                self.perfil_p = ones(1, platos).*alim{2}.P;
            elseif nargin > 5 && isa(Pperfil, 'cell') && ~isempty(platos)
            %Si se provee una celda, se supone que cada 2do elemento contiene presión y cada elem. 1ero
            %contiene un número de plato, se tiene que proveer necesariamente presión del condensador y el rehervidor
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
            elseif nargin > 5  && ~isempty(platos) && length(Pperfil) == platos
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
            if nargin > 7 && ~isempty(num_iter)
                self.iteraciones = num_iter;
            else
                self.iteraciones = 500;
            end
            if nargin > 8 && ~isempty(dobflow) && strcmpi(d_o_b, 'D')
                self.dflujo = dobflow;
                tamano = length(self.alimentaciones);
                self.F = 0;
                for ooi = 2:2:tamano
                    self.F = self.F + self.alimentaciones{ooi}.molF;
                end
                self.UmW = 0;
                if ~isempty(salid) && isa(salid, 'cell')
                    tamano = length(salid);
                    for ooiiu = 1:3:tamano
                        if salid{ooiiu} ~= 1
                            self.UmW = self.UmW + salid{ooiiu+2};
                        end
                    end
                end
                self.bflujo = self.F - self.UmW -self.dflujo;
            elseif nargin > 8 && ~isempty(dobflow) && strcmpi(d_o_b, 'B')
                self.bflujo = dobflow;
                tamano = length(self.alimentaciones);
                self.F = 0;
                for ooi = 2:2:tamano
                    self.F = self.F + self.alimentaciones{ooi}.molF;
                end
                self.UmW = 0;
                if ~isempty(salid) && isa(salid, 'cell')
                    tamano = length(salid);
                    for ooiiu = 1:3:tamano
                        if salid{ooiiu} ~= 1
                            self.UmW = self.UmW + salid{ooiiu+2};
                            
                        end
                    end
                end
                self.dflujo = self.F - self.UmW -self.bflujo;
            end
            if nargin > 10 && ~isempty(perfil_vi) && ~isempty(platos)
                self.perfil_vi = perfil_vi;
                if isempty(salid) || (isa(salid, 'cell') && salid{1} ~= 1)
                    self.dsubi = sum(perfil_vi(1,:));
                elseif (isa(salid, 'cell') && salid{1} == 1)
                    if ~isempty(perfil_li)
                        self.dsubi = perfil_vi(1,:) + perfil_li(1,:);
                    end
                end
                for i = 1:platos
                    self.ysubj(i,:) = self.perfil_vi(i,:)/sum(self.perfil_vi(i,:));
                    self.perfil_v(i) = sum(self.perfil_vi(i,:)); 
                end
            end
            if nargin > 11 && ~isempty(perfil_li)
                self.perfil_li = perfil_li;
                self.bsubi = perfil_li(end,:);
                
                for i = 1:platos
                    self.xsubj(i,:) = self.perfil_li(i,:)/sum(self.perfil_li(i,:));
                    self.perfil_l(i) = sum(self.perfil_li(i,:)); 
                end
            end
            if nargin > 1 && ~isempty(self.etapas) && ~isempty(self.alimentaciones)
                flujofj = 0;
                for i = 2:2:length(self.alimentaciones)
                    flujofj = flujofj + (self.alimentaciones{i}.molF);
                end
                self.convergence = self.etapas*(2*length(self.sust) + 1)*(flujofj)^2*1e-6;
            end
            if nargin > 17
                if ~isempty(varargin)
                    tamano = length(varargin);
                    for i=2:2:tamano
                        if strcmpi(varargin{i}, 't')
                            self.perfil_t = varargin{i-1};
                        elseif strcmpi(varargin{i}, 'k')
                            self.perfil_k = varargin{i-1};
                        elseif strcmpi(varargin{i}, 'di')
                            self.dsubi = varargin{i-1};
                            if isempty(self.salidas) ||  ~isa(self.salidas,'cell')
                                self.dflujo = sum(self.dsubi);
                            elseif self.salidas{1} ~= 1
                                self.dflujo = sum(self.dsubi);
                            end
                        elseif strcmpi(varargin{i}, 'bi')
                            self.bsubi = varargin{i-1};
                        elseif strcmpi(varargin{i}, 'vi')
                            self.perfil_vi = varargin{i-1};
                            for oou = 1:self.etapas
                                self.ysubj(oou,:) = self.perfil_vi(oou,:)/sum(self.perfil_vi(oou,:));
                                self.perfil_v(oou) = sum(self.perfil_vi(oou,:)); 
                            end
                        elseif strcmpi(varargin{i}, 'v')
                            self.perfil_v = varargin{i-1};
                        elseif strcmpi(varargin{i}, 'l' )
                            self.perfil_l = varargin{i-1};
                        elseif strcmpi(varargin{i}, 'li')
                            self.perfil_li = varargin{i-1};
                            for oou = 1:self.etapas
                                self.xsubj(oou,:) = self.perfil_li(oou,:)/sum(self.perfil_li(oou,:));
                                self.perfil_l(oou) = sum(self.perfil_li(oou,:)); 
                            end
                        elseif strcmpi(varargin{i}, 'tj')
                            self.perfil_t = varargin{i-1};
                        end
                    end
                end
            end
            if nargin > 16 && ~isempty(cond1) && ~isempty(cond2) && isa(whoscond2, 'char') && isa(whoscond1, 'char')
                if strcmpi(whoscond1, 'qc') && strcmpi(whoscond2, 'qb')
                    self.qc = cond1;
                    
                    self.qb = cond2;
                end
                if strcmpi(whoscond1, 'R') && strcmpi(whoscond2, 'B')
                    self.bflujo = cond2;
                    if isempty(self.UmW) || ~isa(self.UmW, 'double')
                        self.dflujo = self.F - self.UmW - self.bflujo;
                    else
                         self.dflujo = self.F - self.bflujo - self.UmW;
                    end
                    self.reflujo = cond1;
                end
                if strcmpi(whoscond1, 'B') && strcmpi(whoscond2, 'R')
                    self.bflujo = cond1;
                    if isempty(self.UmW) || ~isa(self.UmW, 'double')
                        self.dflujo = self.F - self.bflujo;
                    else
                         self.dflujo = self.F - self.bflujo - self.UmW;
                    end
                    self.reflujo = cond2;
                end
                if strcmpi(whoscond1, 'R') && strcmpi(whoscond2, 'D')
                    self.dflujo = cond2;
                    if isempty(self.UmW) || ~isa(self.UmW, 'double')
                        self.bflujo = self.F - self.dflujo;
                    else
                        self.bflujo = self.F - self.UmW - self.dflujo;
                    end
                    self.reflujo = cond1;
                end
                if strcmpi(whoscond1, 'D') && strcmpi(whoscond2, 'R')
                    self.dflujo = cond2;
                    if isempty(self.UmW) || ~isa(self.UmW, 'double')
                        self.bflujo = self.F - self.dflujo;
                    else
                        self.bflujo = self.F - self.UmW - self.dflujo;
                    end
                    self.reflujo = cond1;
                end
                if strcmpi(whoscond1, 'bi') && strcmpi(whoscond2, 'di')
                    hk = cond1{1};
                    bidato = cond1{2};
                    lk = cond2{1};
                    didato = cond2{2};
                    if isempty(self.reflujo) || self.reflujo == 0
                        self.reflujo = 3;
                    end
                    if ~isempty(self.salidas) &&  isa(self.salidas, 'cell')
                        saliddd = self.salidas;
                        self.salidas = double.empty(0,3);
                    end
                    if isempty(self.dflujo) || self.dflujo == 0
                        self.dflujo = 0;
                    
                        for ittyhr = 1:hk
                            if ittyhr == lk
                                self.dflujo = self.dflujo + didato;
                            elseif ittyhr == hk
                                self.dflujo = self.dflujo + self.Fi(ittyhr) - bidato;
                            else
                                self.dflujo = self.dflujo + self.Fi(ittyhr);
                            end
                        end
                    end
                    self.balanmasa();
                    if isempty(self.salidas) && ~isa(self.salidas, 'cell') 
                        self.salidas = saliddd;
                    end
                    self.lk = lk;
                    self.hk = hk;
                    self.bsubi(self.hk) = bidato;
                    self.dsubi(self.lk) = didato;
                    self.bsubi(self.lk) = self.Fi(self.lk) - didato;
                    self.dsubi(self.hk) = self.Fi(self.hk) - bidato;
                    if ~isempty(self.salidas) && isa(self.salidas, 'cell')
                        self.salidas{3} = sum(self.dflujo);
                    end
                end
                if strcmpi(whoscond1, 'di') && strcmpi(whoscond2, 'bi')
                    hk = cond2{1};
                    bidato = cond2{2};
                    lk = cond1{1};
                    didato = cond1{2};
                    if isempty(self.reflujo) || self.reflujo == 0
                        self.reflujo = 3;
                    end
                    if isempty(self.salidas) && isa(self.salidas, 'double')
                        saliddd = self.salidas;
                        self.salidas = double.empty(0,3);
                    end
                    if isempty(self.dflujo) || self.dflujo == 0
                        self.dflujo = 0;
                        for ittyhr = 1:lk
                            if ittyhr == lk
                                self.dflujo = self.dflujo + didato;
                            elseif ittyhr == hk
                                self.dflujo = self.dflujo + self.Fi(ittyhr) - bidato;
                            else
                                self.dflujo = self.dflujo + self.Fi(ittyhr);
                            end
                        end
                    end
                    self.balanmasa();
                    self.lk = lk;
                    self.hk = hk;
                    if ~isempty(self.salidas) && ~isa(self.salidas, 'cell') 
                        self.salidas = saliddd;
                    end
                    self.bsubi(self.hk) = bidato;
                    self.dsubi(self.lk) = didato;
                    self.bsubi(self.lk) = self.Fi(self.lk) - didato;
                    self.dsubi(self.hk) = self.Fi(self.hk) - bidato;
                    if isempty(self.salidas) && isa(self.salidas, 'double')
                        self.salidas{3} = sum(self.dflujo);
                    end
                end
                if strcmpi(whoscond1, 'qc') 
                    self.qc = cond1;
                end
                if strcmpi(whoscond2, 'qb') 
                    self.qb = cond2;
                end
            elseif nargin > 16 && isempty(cond1) && ~isempty(cond2) && isa(whoscond2, 'char') 
                if strcmpi(whoscond2, 'qc') && strcmpi(whoscond2, 'qb')
                    self.qc = cond1;
                    self.qb = cond2;
                elseif strcmpi(whoscond2, 'B') && strcmpi(whoscond1, 'R')  
                    self.bflujo = cond2;
                    if isempty(self.UmW) || ~isa(self.UmW, 'double')
                        self.dflujo = self.F - self.UmW - self.bflujo;
                    else
                         self.dflujo = self.F - self.bflujo - self.UmW;
                    end
                    self.reflujo = cond1;
                elseif strcmpi(whoscond2, 'R') && strcmpi(whoscond1, 'B') 
                    self.bflujo = cond2;
                    if isempty(self.UmW) || ~isa(self.UmW, 'double')
                        self.dflujo = self.F - self.UmW - self.bflujo;
                    else
                         self.dflujo = self.F - self.bflujo - self.UmW;
                    end
                    self.reflujo = cond1;
                elseif strcmpi(whoscond1, 'D') && strcmpi(whoscond2, 'R')
                    self.dflujo = cond2;
                    if isempty(self.UmW) || ~isa(self.UmW, 'double')
                        self.bflujo = self.F - self.dflujo;
                    else
                        self.bflujo = self.F - self.UmW - self.dflujo;
                    end
                    self.reflujo = cond1;
                elseif strcmpi(whoscond2, 'D') && strcmpi(whoscond1, 'R')
                    self.dflujo = cond2;
                    if isempty(self.UmW) || ~isa(self.UmW, 'double')
                        self.bflujo = self.F - self.dflujo;
                    else
                        self.bflujo = self.F - self.UmW - self.dflujo;
                    end
                    self.reflujo = cond1;
                end
                
            elseif nargin > 15 && ~isempty(cond1) && isa(whoscond1, 'char') 
                if strcmpi(whoscond1, 'qc')
                    self.qc = cond1;
                end
                if strcmpi(whoscond1, 'qb') 
                    self.qb = cond2;
                end
            end
        end
        function definealim(self, platos, alim)
            if nargin > 1 && ~isempty(alim) && isa(alim, 'cell') && isa(alim{2}, 'Corriente')
                long = length(alim);
                self.entradas = cell.empty(0,0);
                tamano = length(alim);
                self.F = 0;
                for ooi = 2:2:tamano
                    self.F = self.F + alim{ooi}.molF;
                end
                for i = 2:2:long
                    corr = alim{i};
                    susta = corr.comp;
                    %Extraigo las clases sustancia no repetidas a la mezcla de la torre
                    for j = 1:length(susta)
                        if ~any(strcmp(self.comp, susta(j).id))
                            self.comp{length(self.comp) + 1} = susta(j).id;
                            
                        end
                    end
                    %El primer elemento es la etapa en donde entra la alimentacion
                    self.entradas{length(self.entradas) + 1} = alim{i-1};
                    %El segundo elemento es la condición térmica de alimentación 
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
                    %Extraigo las concentraciones de cada una de las sustancias en la mezcla
                    for j = 1:length(susta)
                        indice = find(strcmp(susta(j).id, self.comp),1, 'first');
                        self.concfeed(floor(i/2), indice) = corr.conc(j);
                        %Las contengo en una matriz "concfeed" de tantas filas como 
                        %alimentaciones y de tantas columnas como componentes tienen 
                        %las alimentaciones ordenadas en volatilidades
                    end
                end
                self.xsubj = zeros(platos, length(self.comp));
                self.ysubj = zeros(platos, length(self.comp));
                self.alimentaciones = alim;
                %El usuario debe usar un mismo modelo termo en todas las alim, así que se extrae de la principal
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
        end
        function newtonraphson(self)
            % Se construye la matriz tridiagonal de derivadas
            %  [ B(1)  C(1)                                  ] [  x(1)  ]   [  F(1)  ]
            %  [ A(2)  B(2)  C(2)                            ] [  x(2)  ]   [  F(2)  ]
            %  [       A(3)  B(3)  C(3)                      ] [        ]   [        ]
            %  [            ...   ...   ...                  ] [  ...   ] = [  ...   ]
            %  [                    ...    ...    ...        ] [        ]   [        ]
            %  [                        A(n-1) B(n-1) C(n-1) ] [ x(n-1) ]   [ F(n-1) ]
            %  [                                 A(n)  B(n)  ] [  x(n)  ]   [  F(n)  ]
            % en ella, las matrices B(j) son la contribuci�n al plato j
            % de las variables de estado de los platos j-1, j y j+1
            % Las matrices A(j-1) son la dependencia con el plato j-1 del factor
            % de cambio de las variables de estado del plato j
            % Las matrices C(j+1) son la dependencia con el plato j+1 del factor
            % de cambio de las variables de estado del plato j
            %
            % B(j), A(j-1), C(j+1)  son matrices de la forma
            %                 1                          2*n+1
            % dFj/dX        X = v1   � X = lij          X = Tj
            % 1     F = Hj [dHj/dv1  � dHj/dlij  �  � dHj/dTj  ]                             ]
            % 2     F = M1j[dM1j/dv1 � dM1j/dlij �  � dM1j/dTj ]
            % �           [   �        dFj/dX     �      �    ] 
            % 2*n+1        [   �          �       �       �   ]
%             if isempty(self.perfil_k) || isempty(self.perfil_v) || isempty(self.perfil_t)
%                 self = self.balanmasa();
%                 self.generar_etapas();
%             end
            if isempty(self.perfil_k) || isempty(self.perfil_vi) || isempty(self.perfil_t) || isempty(self.perfil_li)
                self = self.balanmasa();
            else
                self.generar_etapas();
            end
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
            self.perfil_v(1) = self.dflujo;
            if ~isempty(self.salidas) && isa(self.salidas, 'cell') 
                if self.salidas{1} == 1 && self.salidas{2} == 1
                    self.dflujo = self.dflujo + self.salidas{3};
                elseif self.salidas{1} == 1 && self.salidas{2} == 0
                    self.dflujo = self.dflujo + self.salidas{3};
                end
            end
            if isempty(self.tol)
                self.tol = self.convergence;
            end
            if isempty(self.damping) || ~isa(self.damping, 'double') 
                self.damping = 1;
            elseif isa(self.damping, 'double') && (self.damping < 0  || self.damping > 2) 
                self.damping = 1;
            end
            nuevas_varX = zeros((2.*self.num_sust + 1)*self.etapas, 1);
            self.varX = nuevas_varX;
            Vbackup = zeros(1, self.etapas);
            Tbackup = zeros(1, self.etapas);
            Lbackup = zeros(1, self.etapas);
            for iitter0 = 1 : self.etapas
                self.varX(1+(iitter0-1)*(2*(self.num_sust)+1):1+(iitter0-1)*(2*(self.num_sust)+1) + self.num_sust - 1) = self.perfil_vi(iitter0,:);
                self.varX(1+self.num_sust+(iitter0-1)*(2*(self.num_sust)+1):1+(iitter0-1)*(2*(self.num_sust)+1) + 2*self.num_sust-1) = self.perfil_li(iitter0,:);
                self.varX(1+(iitter0-1)*(2*(self.num_sust)+1)+2*self.num_sust) = self.perfil_t(iitter0);
            end
            
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
            valI = 1e308;
            iter = 0;
            while valI > self.tol && iter < self.iteraciones
                valI = 0;
                iter = iter + 1;
                self.perfil_k
                if iter == 1
                    for i = 1:self.etapas
                        self.platos(i).L = self.perfil_l(i);
                        self.platos(i).V = self.perfil_v(i);
                        self.platos(i).y_i = self.ysubj(i,:);
                        self.platos(i).x_i = self.xsubj(i,:);
                        self.platos(i).v_i = self.perfil_vi(i,:);
                        self.platos(i).l_i = self.perfil_li(i,:);
                    end
                end
                
                self.Bij = zeros(2*self.num_sust +1, (self.etapas)*(2*self.num_sust + 1));
                self.Aij = zeros(2*self.num_sust +1, (self.etapas-1)*(2*self.num_sust + 1));
                self.Cij = zeros(2*self.num_sust +1, (self.etapas-1)*(2*self.num_sust + 1));
                self.Fk = zeros((2*self.num_sust+1).*self.etapas, 1);
                tamano = size(self.entradas);

                for etapa = 1:self.etapas 
                    if etapa < self.etapas
                        for dMidli = 1:self.num_sust
                            self.Cij(dMidli+1, (etapa-1)*(2*self.num_sust + 1)+dMidli) = -1;
                        end
                    end
                    if etapa > 1
                        for dMidli = 1:self.num_sust
                            self.Aij(dMidli+1, self.num_sust + (etapa-2)*(2*self.num_sust + 1)+dMidli) = -1;
                        end
                    end
                end
                self.platos_entradas = zeros(1, length(self.entradas)/3);
                alimen = Corriente.empty(0, length(self.entradas)/3 );
                for itere=1:3:length(self.entradas)
                    self.platos_entradas((itere-1)/3+1) = self.entradas{itere};
                end
                for itere = 1:2:length(self.alimentaciones)
                    alimen((itere-1)/2+1) = self.alimentaciones{itere + 1};
                end

                for num_plato = 1:self.etapas
                    for l = 1:self.num_sust
                        for m = 3:3:tamano(2)
                            if num_plato == self.entradas{m-2}
                                self.Fk(1+l+(num_plato - 1)*(2.*self.num_sust + 1)) = self.Fk(1+l+(num_plato - 1)*(2.*self.num_sust + 1)) + self.entradas{m}.*self.concfeed(m/3, l);
                            end
                        end
                    end
                    if num_plato ~= 1
                        Xcj = self.platos(num_plato-1).x_i;
                        TLm1 = self.platos(num_plato-1).T;
                        PLm1 = self.platos(num_plato-1).P;
                        mezclaXm1 = Mezcla(self.sust, Xcj, self.alimentaciones{2}.mezcla.kij);
                        HgiLm1 = zeros(1, self.num_sust);
                    else 
                        HgiLm1 = 0;
                    end
                    if num_plato ~= self.etapas
                        Ycj = self.platos(num_plato+1).y_i;
                        TVp1 = self.platos(num_plato+1).T;
                        PVp1 = self.platos(num_plato+1).P;
                        mezclaYp1 = Mezcla(self.sust, Ycj, self.alimentaciones{2}.mezcla.kij);
                        HgiVp1 = zeros(1, self.num_sust);
                    else
                        HgiVp1 = 0;
                    end
                    T = self.varX(1+(num_plato-1)*(2*(self.num_sust)+1)+2*self.num_sust);   % Temperatura Alternativamente self.platos(num_plato).T
                    P = self.platos(num_plato).P;
                    mezclaX = Mezcla(self.sust, self.platos(num_plato).x_i, self.alimentaciones{2}.mezcla.kij);
                    mezclaY = Mezcla(self.sust, self.platos(num_plato).y_i, self.alimentaciones{2}.mezcla.kij);
                    HgiL = 0;
                    HgiV = 0;
                    Href = zeros(1, self.num_sust);
                    for i = 1:self.num_sust
                        Href(i) = self.sust(i).href;
                        try 
                            cp = self.sust(i).cp_gi{1};
                        catch ME
                            error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                        end
                        if num_plato ~= 1
                            deltaHLm1 = integral(@(t) cp(t), 273.15, TLm1);
                            HgiLm1(i) = deltaHLm1*mezclaXm1.conc(i) + Href(i)*mezclaXm1.conc(i);
                        end
                        if num_plato ~= self.etapas
                            deltaHVp1 = integral(@(t) cp(t), 273.15, TVp1);
                            HgiVp1(i) = deltaHVp1*mezclaYp1.conc(i) + Href(i)*mezclaYp1.conc(i);
                        end
                        deltaHL = integral(@(t) cp(t), 273.15, T);
                        HgiL =  HgiL + deltaHL*mezclaX.conc(i) + Href(i)*mezclaX.conc(i);
                        deltaHV = deltaHL;
                        HgiV = HgiV + deltaHV*mezclaY.conc(i) + Href(i)*mezclaY.conc(i);
                    end
                    if num_plato ~= 1
                        HgiLm1 = sum(HgiLm1);
                    else
                        HgiLm1 = 0;
                    end
                    if num_plato ~= self.etapas
                        HgiVp1 = sum(HgiVp1);
                    end
                    if num_plato ~= 1
                        Hdep_bubX = self.MEdE.entalpia(TLm1, PLm1, mezclaXm1, 'liq');
                        Hdep_refX = self.MEdE.entalpia(273.15, 101.325, mezclaXm1, 'liq');
                        HLjm1 = HgiLm1 - Hdep_bubX + Hdep_refX; % HLjm1-1
                    else
                        HLjm1 = 0;
                    end
                    if num_plato ~= self.etapas
                        Hdep_dewY = self.MEdE.entalpia(TVp1, PVp1, mezclaYp1, 'vap');
                        Hdep_refY = self.MEdE.entalpia(273.15, 101.325, mezclaYp1, 'liq');
                        HVjp1 = HgiVp1 - Hdep_dewY + Hdep_refY;  % HVjp1+1
                    else 
                        HVjp1 = 0;
                    end
                    
                    Hdep_bubX = self.MEdE.entalpia(T, P, mezclaX, 'liq');
                    Hdep_refX = self.MEdE.entalpia(273.15, 101.325, mezclaX, 'liq');
                    HLj = HgiL - Hdep_bubX + Hdep_refX; % HLj 
                    if num_plato ~= 1
                        Hdep_dewY = self.MEdE.entalpia(T, P, mezclaY, 'vap');
                        Hdep_refY = self.MEdE.entalpia(273.15, 101.325, mezclaY, 'liq');
                        HVj = HgiV - Hdep_dewY + Hdep_refY;  % HVj
                    elseif self.platos(num_plato).V < 1e-5
                        HVj = 0;
                    else
                        Hdep_dewY = self.MEdE.entalpia(T, P, mezclaY, 'vap');
                        Hdep_refY = self.MEdE.entalpia(273.15, 101.325, mezclaY, 'liq');
                        HVj = HgiV - Hdep_dewY + Hdep_refY;  % HVj
                    end
                    %+HLjm1-1 se necesita en la dHj-1/dLj-1 y +HVjp1+1 se necesita en dHj/dVj+1
                    if ~isempty(self.platos(num_plato).salidaV)
                        if num_plato ~=1
                            if self.platos(num_plato).salidaV > 1e-4
                                Sj = self.platos(num_plato).salidaV/self.platos(num_plato).V;
                            else
                                Sj = 0;
                            end
                        else
                            Sj = 0;
                        end
                    else
                        Sj = 0;
                    end
                    if ~isempty(self.platos(num_plato).salidaL);
                        if num_plato ~= self.etapas
                            if self.platos(num_plato).salidaL > 1e-4
                                sj = self.platos(num_plato).salidaL/self.platos(num_plato).L;
                            else
                                sj = 0;
                            end
                        else
                            sj = 0;
                        end
                    else
                        sj = 0;
                    end
                    mezclaY.conc = self.platos(num_plato).y_i;
                    mezclaX = mezclaY;
                    mezclaX.conc = self.platos(num_plato).x_i;
                    %por diferenciacion numerica de 3 puntos centrada de Ej = 0 = Kij*lij*(sum(vj))/(sum(lj))
                    % siendo la funci�n de la temperatura solo dependiente de Kij
                    % Kijp1(Tj) = Kij(T0+dif) = Kij(T0 + 1E-7)
                    delta = 1e-8;
                    fGijp1 = self.MEdE.fugF(T+delta, P, mezclaY, 'vap');
                    fLijp1 = self.MEdE.fugF(T+delta, P, mezclaX, 'liq');
                    fGijm1 = self.MEdE.fugF(T-delta, P, mezclaY, 'vap');
                    fLijm1 = self.MEdE.fugF(T-delta, P, mezclaX, 'liq');
                %                 fGij = MEdE.fugF(Tj, P, mezcla, 'vap');
                %                 fLij = MEdE.fugF(Tj, P, mezcla, 'liq');
                    HgiLm1T = 0;
                    HgiLp1T = 0;

                    HgiVm1T = 0;
                    HgiVp1T = 0;

                    for i = 1:self.num_sust
                        try 
                            cp = self.sust(i).cp_gi{1};
                        catch ME
                            error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                        end
                        deltaHLm1T = integral(@(t) cp(t), 273.15, T - delta);
                        deltaHLp1T = integral(@(t) cp(t), 273.15, T + delta);
                        HgiLm1T = HgiLm1T + deltaHLm1T*mezclaX.conc(i);
                        HgiLp1T = HgiLp1T + deltaHLp1T*mezclaX.conc(i);
                        deltaHVm1T = deltaHLm1T;
                        deltaHVp1T = deltaHLp1T;
                        HgiVm1T = HgiVm1T +  deltaHVm1T*mezclaY.conc(i);
                        HgiVp1T = HgiVp1T + deltaHVp1T*mezclaY.conc(i);
                    end

                    Hdep_bubXm1 = self.MEdE.entalpia(T - delta, P, mezclaX, 'liq');
                    Hdep_bubXp1 = self.MEdE.entalpia(T + delta, P, mezclaX, 'liq');
                        Hdep_dewYm1 = self.MEdE.entalpia(T - delta, P, mezclaY, 'vap');
                        Hdep_dewYp1 = self.MEdE.entalpia(T + delta, P, mezclaY, 'vap');
                        dHVjdT = (HgiVp1T - HgiVm1T)/(2*delta) - (Hdep_dewYp1 - Hdep_dewYm1)/(2*delta);
                    
                    dHLjdT = (HgiLp1T - HgiLm1T)/(2*delta) - (Hdep_bubXp1 - Hdep_bubXm1)/(2*delta);

                    Kijp1 = fLijp1 ./ fGijp1;
                    Kijm1 = fLijm1 ./ fGijm1;
                    
                    fGij = self.MEdE.fugF(T, P, mezclaY, 'vap');
                    fLij = self.MEdE.fugF(T, P, mezclaX, 'liq');
                    self.perfil_k(:,num_plato) = fLij ./ fGij;
                    HgiLm1m1T = 0;
                    HgiLm1p1T = 0;       % Hgasideal para calculo L etapa j-1 -delta
                    HgiVp1m1T = 0;
                    HgiVp1p1T = 0;

                    if num_plato > 1
                        for i = 1:self.num_sust
                            try 
                                cp = self.sust(i).cp_gi{1};
                            catch ME
                                error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                            end
                            
                            deltaHLm1m1T = integral(@(t) cp(t), 273.15, TLm1 - delta);
                            deltaHLm1p1T = integral(@(t) cp(t), 273.15, TLm1 + delta);
                            HgiLm1m1T = HgiLm1m1T + deltaHLm1m1T*mezclaXm1.conc(i);
                            HgiLm1p1T = HgiLm1p1T + deltaHLm1p1T*mezclaXm1.conc(i); 
                        end
                    else 
                        HgiLm1p1T = 0;
                        HgiLm1m1T = 0;
                    end
                    if num_plato > 1
                        Hdep_bubXm1 = self.MEdE.entalpia(TLm1 - delta, P, mezclaXm1, 'liq');
                        Hdep_bubXp1 = self.MEdE.entalpia(TLm1 + delta, P, mezclaXm1, 'liq');
                    else
                        Hdep_bubXp1 = 0;
                        Hdep_bubXm1 = 0;
                    end
                    dHLm1dT = (HgiLm1p1T - HgiLm1m1T)/(2*delta) - (Hdep_bubXp1 - Hdep_bubXm1)/(2*delta);   %Diferencial de HL etapa j-1 

                    if num_plato < self.etapas
                        for i = 1:self.num_sust
                            try 
                                cp = self.sust(i).cp_gi{1};
                            catch ME
                                error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                            end
                            deltaHVp1m1T = integral(@(t) cp(t), 273.15, TVp1 - delta);
                            deltaHVp1p1T = integral(@(t) cp(t), 273.15, TVp1 + delta);
                            HgiVp1m1T = HgiVp1m1T + deltaHVp1m1T*mezclaYp1.conc(i);
                            HgiVp1p1T = HgiVp1p1T + deltaHVp1p1T*mezclaYp1.conc(i);
                        end
                    else
                        HgiVp1m1T = 0;
                        HgiVp1p1T = 0;
                    end
                    if num_plato < self.etapas
                        Hdep_dewYp1 = self.MEdE.entalpia(TVp1 + delta, P, mezclaYp1, 'vap');
                        Hdep_dewYm1 = self.MEdE.entalpia(TVp1 - delta, P, mezclaYp1, 'vap');
                    else
                        Hdep_dewYp1 = 0;
                        Hdep_dewYm1 = 0;
                    end
                    dHVp1dT = (HgiVp1p1T - HgiVp1m1T)/(2*delta) - (Hdep_dewYp1 - Hdep_dewYm1)/(2*delta);   %Diferencial de HV etapa j+1
                    
%                     if num_plato == 1
%                         dQjdTjm1 = 0;
%                     else
%                         dQjdTjm1 = (dHLm1dT)*sum((self.platos(num_plato - 1).l_i));
%                     end
%                     if num_plato == self.etapas
%                         dQjdTjp1 = 0;
%                     else
%                         dQjdTjp1 = (dHVp1dT)*sum((self.platos(num_plato + 1).v_i));
%                     end
                    dQjdTjp1 = 0;
                    dQjdTjm1 = 0;
                    if num_plato > 1
                        self.Aij(1, (num_plato - 2)*((2*self.num_sust + 1))+2.*self.num_sust+1) =  -sum(self.platos(num_plato-1).l_i(:))*(dHLm1dT) - dQjdTjm1;
                    end
                    if num_plato < self.etapas
                        self.Cij(1, (num_plato - 1)*((2*self.num_sust + 1))+2.*self.num_sust+1) =  -sum(self.platos(num_plato+1).v_i(:))*(dHVp1dT) - dQjdTjp1;
                    end
%                     if num_plato == 1
%                         dQjdTj = - (dHVjdT)*(1+Sj).*(sum((self.platos(num_plato).v_i))) - (dHLjdT)*(1+sj).*(sum((self.platos(num_plato).l_i)));  
%                         % [(HVjp1p1 - HVjp1m1)/2delta]*Vjp1 - 
%                         % [(HVjp1 - HVjm1)/2delta]*Vj*(1+Sj) -
%                         % [(HLjp1 - HLjm1)/2delta]*Lj*(1+sj)
%                     elseif num_plato == self.etapas
%                         dQjdTj = - (dHVjdT)*(1+Sj).*(sum((self.platos(num_plato).v_i))) - (dHLjdT)*(1+sj).*(sum((self.platos(num_plato).l_i)));  
%                     else
%                         dQjdTj = - (dHVjdT)*(1+Sj).*(sum((self.platos(num_plato).v_i))) - (dHLjdT)*(1+sj).*(sum((self.platos(num_plato).l_i)));  
%                     end
                    dQjdTj = 0;
                    if num_plato ~= 1 
                        self.Bij(1, (num_plato - 1)*((2*self.num_sust + 1))+2*self.num_sust+1) = dHLjdT*(1+sj)*sum(self.platos(num_plato).l_i(:)) + dHVjdT.*(1+Sj)*sum(self.platos(num_plato).v_i(:)) + dQjdTj;
                    else
                        if sum(self.platos(num_plato).v_i(:)) > 1e-5
                            self.Bij(1, 1:self.num_sust - 1) = dHLjdT*(1+sj)*sum(self.platos(num_plato).l_i(:)) + dHVjdT.*(1+Sj)*sum(self.platos(num_plato).v_i(:)) + dQjdTj;
                        else
                            self.Bij(1,1:self.num_sust - 1) = dHLjdT*(1+sj)*sum(self.platos(num_plato).l_i(:)) + dHVjdT.*(1+sj)*sum(self.platos(num_plato).l_i(:)) + dQjdTj;
                        end
                    end
                    for dFidvi = 1:self.num_sust
                        %Funcion Fj siendo Ej
                        dEjdTj = ((Kijp1(dFidvi) - Kijm1(dFidvi))./(2*delta))*(self.platos(num_plato).l_i(dFidvi)*(sum(self.platos(num_plato).v_i(:)))/(sum(self.platos(num_plato).l_i(:))));
                        self.Bij(self.num_sust + 1 + dFidvi, (num_plato-1)*(2*self.num_sust + 1) + 2*self.num_sust + 1) = dEjdTj;
                        %Funcion Fj siendo Hj
                        self.Bij(1, (num_plato-1)*(2*self.num_sust + 1) + dFidvi) = HVj*(1+Sj);
                        %Funcion Fj siendo Mij
                        self.Bij(dFidvi + 1, (num_plato-1)*(2*self.num_sust + 1) + dFidvi) = 1+Sj;

                        %fprintf(1, '%f * %f * 1/(%f) - 1 = %f \n', self.platos(num_plato).K(num_plato, dFidvi),self.platos(num_plato).l_i(num_plato, dFidvi),sum(self.platos(num_plato).l_i(num_plato, :)), self.platos(num_plato).K(num_plato, dFidvi)*self.platos(num_plato).l_i(num_plato, dFidvi)*(1)/(sum(self.platos(num_plato).l_i(num_plato, :))) - 1)
                        self.Bij(self.num_sust + 1 + dFidvi, (num_plato-1)*(2*self.num_sust + 1) + dFidvi) = self.platos(num_plato).K(dFidvi)*self.platos(num_plato).l_i(dFidvi)*(1)/(sum(self.platos(num_plato).l_i(:))) - 1;
                        if num_plato < self.etapas
                            self.Cij(1, (num_plato-1)*(2*self.num_sust + 1)+dFidvi) = -HVjp1;
                        end
                    end
                    for dFidvi = self.num_sust:2
                        for dFidvj = 1:1:self.num_sust
                            self.Bij(self.num_sust + dFidvi, (num_plato-1)*(2*self.num_sust + 1) + dFidvi) =  self.platos(num_plato).K(dFidvj)*self.platos(num_plato).l_i(dFidvj)*(1)/(sum(self.platos(num_plato).l_i(:)));
                        end
                    end
                    for dFidvi = 1:self.num_sust
                        for dFidvj = 1:1:self.num_sust
                            if dFidvj ~=  dFidvi 
                                self.Bij(self.num_sust + 1 + dFidvj, (num_plato-1)*(2*self.num_sust + 1) + dFidvi) = self.platos(num_plato).K(dFidvj)*self.platos(num_plato).l_i(dFidvj)*(1)/(sum(self.platos(num_plato).l_i(:)));
                            end
                        end
                    end
                    for dFidli = 1:self.num_sust
                        %Funcion Fj siendo Hj
                        self.Bij(1, self.num_sust + (num_plato - 1)*(2*self.num_sust + 1)+dFidli) = HLj*(1+sj);
                        %Funcion Fj siendo Mij
                        self.Bij(dFidli + 1, (num_plato-1)*(2*self.num_sust + 1) + dFidli + self.num_sust) = 1+sj;
                        if num_plato > 1
                            self.Aij(1, self.num_sust + (num_plato-2)*(2*self.num_sust + 1)+dFidli) = -HLjm1;
                        end
                        %fprintf(1, '%f\n', self.platos(num_plato).K(num_plato, dFidli)*((sum(self.platos(num_plato).l_i(num_plato, :))-self.platos(num_plato).l_i(num_plato, dFidli))*(sum(self.platos(num_plato).v_i(num_plato, :))))/(sum(self.platos(num_plato).l_i(num_plato, :))^2))
                        self.Bij(self.num_sust + 1 + dFidli, self.num_sust + (num_plato - 1)*(2*self.num_sust + 1)+dFidli) = self.platos(num_plato).K(dFidli)*((sum(self.platos(num_plato).l_i(:))-self.platos(num_plato).l_i(dFidli))*(sum(self.platos(num_plato).v_i(:))))/(sum(self.platos(num_plato).l_i(:))^2);
                    end
                    for dFidli = 1:self.num_sust
                        for dFidlj = self.num_sust:-1:1
                            %Funcion Fj siendo Ej
                            if dFidlj ~= dFidli
                                self.Bij(1 + dFidlj + self.num_sust, self.num_sust + (num_plato-1)*(2*self.num_sust + 1) + dFidli) = -self.platos(num_plato).K(dFidlj)*((sum(self.platos(num_plato).v_i(:)))*(self.platos(num_plato).l_i(dFidlj)))/(sum(self.platos(num_plato).l_i(:))^2);
                            end
                        end
                    end
                    for i = 1:self.num_sust
                        self.Fk(1+i+(num_plato - 1)*(2.*self.num_sust + 1)) = self.Fk(1+i+(num_plato - 1)*(2.*self.num_sust + 1)) - self.platos(num_plato).l_i(i)*(1+sj) - self.platos(num_plato).v_i(i)*(1+Sj);
                        if num_plato > 1
                            self.Fk(1+i+(num_plato - 1)*(2.*self.num_sust + 1)) =  self.Fk(1+i+(num_plato - 1)*(2.*self.num_sust + 1)) + self.platos(num_plato-1).l_i(i);
                        end
                        if num_plato < self.etapas
                            self.Fk(1+i+(num_plato - 1)*(2.*self.num_sust + 1)) =  self.Fk(1+i+(num_plato - 1)*(2.*self.num_sust + 1)) + self.platos(num_plato+1).v_i(i);
                        end      
                    end
                    if num_plato == self.etapas
                        self.Fk(1+(num_plato - 1)*(2.*self.num_sust + 1)) = self.Fk(1+(num_plato - 1)*(2.*self.num_sust + 1)) + HLjm1*(sum(self.platos(num_plato-1).l_i(:)));
                    end
                    if num_plato == 1
                        self.Fk(1+(num_plato - 1)*(2.*self.num_sust + 1)) = self.Fk(1+(num_plato - 1)*(2.*self.num_sust + 1)) + HVjp1*(sum(self.platos(num_plato+1).v_i(:)));
                    end
                    for i = 1:self.num_sust
                        self.Fk(1+self.num_sust+i +(num_plato - 1)*(2.*self.num_sust + 1)) = self.Fk(1+self.num_sust+i +(num_plato - 1)*(2.*self.num_sust + 1)) - self.platos(num_plato).K(i)*self.platos(num_plato).l_i(i)*((sum(self.platos(num_plato).v_i(:)))./(sum(self.platos(num_plato).l_i(:))))+self.platos(num_plato).v_i(i);     
                    end
                    indicee = 1;
                    if  any(num_plato == self.platos_entradas(:))    %Entalpia de alimentaci�n HF
                        if length(self.platos_entradas) == 1
                            indicee = 1;
                        else
                            indicee = find((self.platos_entradas(:) == num_plato) ~= 0, 1, 'first');
                        end
                        HF = alimen(indicee).H;
                    else
                        HF = 0;
                    end
                    if num_plato == 1
                        self.Fk(1+(num_plato - 1)*(2.*self.num_sust + 1)) = self.Fk(1+(num_plato - 1)*(2.*self.num_sust + 1))  - HLj*(1+sj)*(sum(self.platos(num_plato).l_i(:))) - HVj*(1+Sj)*(sum(self.platos(num_plato).v_i(:))) + HF*(sum(alimen(indicee).conc(:)))*alimen(indicee).molF - self.qc;
                    elseif num_plato == self.etapas
                        self.Fk(1+(num_plato - 1)*(2.*self.num_sust + 1)) = self.Fk(1+(num_plato - 1)*(2.*self.num_sust + 1))  - HLj*(1+sj)*(sum(self.platos(num_plato).l_i(:))) - HVj*(1+Sj)*(sum(self.platos(num_plato).v_i(:))) + HF*(sum(alimen(indicee).conc(:)))*alimen(indicee).molF - self.qb;
                    else
                        self.Fk(1+(num_plato - 1)*(2.*self.num_sust + 1)) = self.Fk(1+(num_plato - 1)*(2.*self.num_sust + 1))  - HLj*(1+sj)*(sum(self.platos(num_plato).l_i(:))) - HVj*(1+Sj)*(sum(self.platos(num_plato).v_i(:))) + HF*(sum(alimen(indicee).conc(:)))*alimen(indicee).molF + self.perfil_q(num_plato-1) + HLjm1 * sum(self.platos(num_plato-1).l_i(:)) + HVjp1 * sum(self.platos(num_plato + 1).v_i(:));
                    end
                end
                
                [deltaXvar, self.matriz_tridiag ] = tridiagThomas(self.Bij, self.Aij, self.Cij, self.Fk);
                condition = 0;
                for num_plato = 1:self.etapas
                    nuevas_varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust) = self.varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust) + self.damping.*deltaXvar(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust);   
                    if condition == 1 || any(nuevas_varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust-1) < 0)
                        condition = 1;
                        nuevas_varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust - 1) = self.varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust-1).*exp((self.damping.*deltaXvar(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust-1))./(self.varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust-1)));
                    end
                    if num_plato == self.etapas && condition == 1
                        condition = 0;
                    end
                    Tbackup(num_plato) = self.perfil_t(num_plato);
                    Vbackup(num_plato) = self.perfil_v(num_plato);
                    Lbackup(num_plato) = self.perfil_l(num_plato);
                    self.perfil_t(num_plato) = nuevas_varX(1+(num_plato-1)*(2*(self.num_sust)+1)+2*self.num_sust);
                    if sum(nuevas_varX(1+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + self.num_sust-1)) > self.tol
                        self.perfil_v(num_plato) = sum(nuevas_varX(1+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + self.num_sust-1));
                    end
                    self.perfil_l(num_plato) = sum(nuevas_varX(1+self.num_sust+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + 2*self.num_sust-1));

                    valI = valI + ((self.perfil_t(num_plato) - Tbackup(num_plato))/self.perfil_t(num_plato))^2;
                    if self.perfil_v(num_plato) > 0.00001 * self.tol
                        valI = valI + ((self.perfil_v(num_plato) - Vbackup(num_plato))/self.perfil_v(num_plato))^2;
                    end
                    if self.perfil_l(num_plato) > 0.00001 * self.tol
                        valI = valI + ((self.perfil_l(num_plato) - Lbackup(num_plato))/self.perfil_l(num_plato))^2;
                    end
                    self.varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust) = nuevas_varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust);

                    self.perfil_vi(num_plato, :) = self.varX(1+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + self.num_sust-1)';
                    self.perfil_li(num_plato, :) = self.varX(1+self.num_sust+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + 2*self.num_sust-1)';
                    if num_plato == 1
                        if ~isempty(self.salidas) && isa(self.salidas, 'cell')
                            if self.salidas{1} == 1 && self.salidas{2} == 1 && self.perfil_v(num_plato)/(self.perfil_v(num_plato)+self.perfil_l(num_plato)) < 1e-4
                                self.ysubj(num_plato,:) = self.varX(1+(num_plato)*(2*(self.num_sust)+1):1+(num_plato)*(2*(self.num_sust)+1) + self.num_sust-1)'./sum(self.varX(1+(num_plato)*(2*(self.num_sust)+1):1+(num_plato)*(2*(self.num_sust)+1) + self.num_sust-1));
                            else
                                self.ysubj(num_plato,:) = self.varX(1:self.num_sust)'./sum(self.varX(1:self.num_sust));
                            end
                        else
                            self.ysubj(num_plato,:) = self.varX(1:self.num_sust)'./sum(self.varX(1:self.num_sust));
                        end
                    else
                        self.ysubj(num_plato,:) = self.varX(1+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + self.num_sust-1)'./sum(self.varX(1+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + self.num_sust-1));
                    end
                    self.xsubj(num_plato,:) = self.varX(1+self.num_sust+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + 2*self.num_sust-1)'./sum(self.varX(1+self.num_sust + (num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + 2*self.num_sust-1));                    
                    self.ysubj(num_plato,:) = self.ysubj(num_plato,:);
                    self.xsubj(num_plato,:) = self.xsubj(num_plato,:);
                    self.perfil_t(num_plato) = self.varX(1+(num_plato-1)*(2*(self.num_sust)+1)+2*self.num_sust);
                    
                    fprintf(1,  '%f ', Tbackup(num_plato));
                end
                %fprintf(1, 'Tbackup = %f ', Tbackup(1));
                %for i = 2:self.etapas - 1
                %    fprintf(1, ' %f ', Tbackup(i));
                %end
                %fprintf(1, ' %f \n', Tbackup(end));
                for iiiiter = 1:self.etapas
                    if iiiiter == self.etapas
                        fprintf(1, ' %f \n', self.perfil_t(iiiiter));
                        break
                    end
                    if iiiiter == 1
                        fprintf(1, '\n %f ', self.perfil_t(iiiiter));
                    else
                        fprintf(1, ' %f ', self.perfil_t(iiiiter));
                    end
                    
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
                
                
                for i = 1:self.etapas
                    self.platos(i).L = self.perfil_l(i);
                    self.platos(i).V = self.perfil_v(i);
                    self.platos(i).y_i = self.ysubj(i,:);
                    self.platos(i).x_i = self.xsubj(i,:);
                    self.platos(i).v_i = self.perfil_vi(i,:);
                    self.platos(i).l_i = self.perfil_li(i,:);
                end
                self.bflujo = sum(self.perfil_li(end,:));
                if isempty(self.salidas) 
                    self.dflujo = sum(self.perfil_vi(1,:)); 
                else
                    if isa(self.salidas, 'cell')
                        if self.salidas{1} == 1 && self.salidas{2} == 1
                            self.dflujo =  sum(self.perfil_vi(1,:)) + self.F - self.UmW - sum(self.perfil_li(end,:));
                        end
                    else
                        self.dflujo = sum(self.perfil_vi(1,:)); 
                    end
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
                self.respaldovalI(self.actualiter) = valI;
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
        end
        function newtonraphson2(self)  %Especificaciones R y B
            % Se construye la matriz tridiagonal de derivadas
            %  [ B(1)  C(1)                                  ] [  x(1)  ]   [  F(1)  ]
            %  [ A(2)  B(2)  C(2)                            ] [  x(2)  ]   [  F(2)  ]
            %  [       A(3)  B(3)  C(3)                      ] [        ]   [        ]
            %  [            ...   ...   ...                  ] [  ...   ] = [  ...   ]
            %  [                    ...    ...    ...        ] [        ]   [        ]
            %  [                        A(n-1) B(n-1) C(n-1) ] [ x(n-1) ]   [ F(n-1) ]
            %  [                                 A(n)  B(n)  ] [  x(n)  ]   [  F(n)  ]
            % en ella, las matrices B(j) son la contribuci�n al plato j
            % de las variables de estado de los platos j-1, j y j+1
            % Las matrices A(j-1) son la dependencia con el plato j-1 del factor
            % de cambio de las variables de estado del plato j
            % Las matrices C(j+1) son la dependencia con el plato j+1 del factor
            % de cambio de las variables de estado del plato j
            %
            %La ecuaci�n de reemplazo del balance de energ�a en la etapa 1
            %es la OD: (Henley & Seader, 2008)
            %
            % OD = sum(li,1) - R*sum(vi,1) = 0
            %
            % La ecuaci�n de reemplazo del balance de energ�a en la etapa N
            % es la OW
            %
            % OW = sum(li,N) - B = 0
            %
            % B(j), A(j-1), C(j+1)  son matrices de la forma
            %                 1                          2*n+1
            % dFj/dX        X = v1   � X = lij          X = Tj
            % 1     F = ODj [dODj/dv1  � dODj/dlij  �  � dODj/dTj]
            % 2     F = M1j [dM1j/dv1 � dM1j/dlij �  � dM1j/dTj ]
            % �             [   �        dFj/dX     �      �    ] 
            % 2*n+1         [   �          �       �       �   ]
%             if isempty(self.perfil_k) || isempty(self.perfil_v) || isempty(self.perfil_t)
%                 self = self.balanmasa();
%                 self.generar_etapas();
%             end
            if isempty(self.perfil_k) || isempty(self.perfil_vi) || isempty(self.perfil_t) || isempty(self.perfil_li)
                self = self.balanmasa();
            else
                self.generar_etapas();
            end
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
            self.perfil_v(1) = self.dflujo;
            if ~isempty(self.salidas) && isa(self.salidas, 'cell') 
                if self.salidas{1} == 1 && self.salidas{2} == 1
                    self.dflujo = self.dflujo + self.salidas{3};
                elseif self.salidas{1} == 1 && self.salidas{2} == 0
                    self.dflujo = self.dflujo + self.salidas{3};
                end
            end
            if isempty(self.tol)
                self.tol = self.convergence;
            end
            if isempty(self.damping) || ~isa(self.damping, 'double') 
                self.damping = 1;
            elseif isa(self.damping, 'double') && (self.damping < 0  || self.damping > 2) 
                self.damping = 1;
            end
            nuevas_varX = zeros((2.*self.num_sust + 1)*self.etapas, 1);
            self.varX = nuevas_varX;
            Vbackup = zeros(1, self.etapas);
            Tbackup = zeros(1, self.etapas);
            Lbackup = zeros(1, self.etapas);
            for iitter0 = 1 : self.etapas
                self.varX(1+(iitter0-1)*(2*(self.num_sust)+1):1+(iitter0-1)*(2*(self.num_sust)+1) + self.num_sust - 1) = self.perfil_vi(iitter0,:);
                self.varX(1+self.num_sust+(iitter0-1)*(2*(self.num_sust)+1):1+(iitter0-1)*(2*(self.num_sust)+1) + 2*self.num_sust-1) = self.perfil_li(iitter0,:);
                self.varX(1+(iitter0-1)*(2*(self.num_sust)+1)+2*self.num_sust) = self.perfil_t(iitter0);
            end
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
            if isempty(self.respaldovalI)
                valI = 1e308;
            else
                valI = self.respaldovalI;
            end
            
            iter = 0;
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
                respaldobi = self.bsubi;
                respaldodi = self.dsubi;
            
                while valI > self.tol && iter < self.iteraciones
                    valI = 0;
                    iter = iter + 1;
                    respaldandoiter = iter;
                    if iter == 1
                        for i = 1:self.etapas
                            self.platos(i).L = self.perfil_l(i);
                            self.platos(i).V = self.perfil_v(i);
                            self.platos(i).y_i = self.ysubj(i,:);
                            self.platos(i).x_i = self.xsubj(i,:);
                            self.platos(i).v_i = self.perfil_vi(i,:);
                            self.platos(i).l_i = self.perfil_li(i,:);
                        end
                    end
                    self.Bij = zeros(2*self.num_sust +1, (self.etapas)*(2*self.num_sust + 1));
                    self.Aij = zeros(2*self.num_sust +1, (self.etapas-1)*(2*self.num_sust + 1));
                    self.Cij = zeros(2*self.num_sust +1, (self.etapas-1)*(2*self.num_sust + 1));
                    self.Fk = zeros((2*self.num_sust+1).*self.etapas, 1);
                    tamano = size(self.entradas);
                    HVj = zeros(1, self.etapas);
                    HLj = zeros(1, self.etapas);
                    sj = zeros(1, self.etapas);
                    Sj = zeros(1, self.etapas);
                    for etapa = 1:self.etapas 
                        if etapa < self.etapas
                            for dMidli = 1:self.num_sust
                                self.Cij(dMidli+1, (etapa-1)*(2*self.num_sust + 1)+dMidli) = -1;
                            end
                        end
                        if etapa > 1
                            for dMidli = 1:self.num_sust
                                self.Aij(dMidli+1, self.num_sust + (etapa-2)*(2*self.num_sust + 1)+dMidli) = -1;
                            end
                        end
                    end
                    self.platos_entradas = zeros(1, length(self.entradas)/3);
                    alimen = Corriente.empty(0, length(self.entradas)/3 );
                    for itere=1:3:length(self.entradas)
                        self.platos_entradas((itere-1)/3+1) = self.entradas{itere};
                    end
                    for itere = 1:2:length(self.alimentaciones)
                        alimen((itere-1)/2+1) = self.alimentaciones{itere + 1};
                    end

                    for num_plato = 1:self.etapas
                        if ~isempty(self.platos(num_plato).salidaL);
                                if num_plato ~= self.etapas
                                    if self.platos(num_plato).salidaL > 1e-4
                                        sj(num_plato) = self.platos(num_plato).salidaL/self.platos(num_plato).L;
                                    else
                                        sj(num_plato) = 0;
                                    end
                                else
                                    sj(num_plato) = 0;
                                end
                        else
                            sj(num_plato) = 0;
                        end
                        if num_plato == 1
                            if sj(num_plato) ~= 0 
                                self.Bij(1,  1:self.num_sust) = -self.reflujo-1;
                                self.Bij(1 , 1+self.num_sust + (num_plato-1)*(2*self.num_sust +1):  1+self.num_sust + (num_plato-1)*(2*self.num_sust +1) + self.num_sust - 1) = sj(num_plato)*(-self.reflujo -1);
                                self.Cij(1, 1:self.num_sust) = 1;
                            else
                                self.Bij(1,  1:self.num_sust) = -self.reflujo;
                                self.Bij(1 , 1+self.num_sust + (num_plato-1)*(2*self.num_sust +1):  1+self.num_sust + (num_plato-1)*(2*self.num_sust +1) + self.num_sust - 1) = 1;
                            end
                        end
                        % La derivada de la funci�n de reemplazo OW es siempre +1
                        % para todos los Li y -reflujo (osea -R) para todos los vi
                        if num_plato == self.etapas
                            self.Bij(1, 1+self.num_sust+(num_plato-1)*(2*self.num_sust + 1):1+self.num_sust+(num_plato-1)*(2*self.num_sust + 1) + self.num_sust-1) = 1;
                        end
                        for l = 1:self.num_sust
                            for m = 3:3:tamano(2)
                                if num_plato == self.entradas{m-2}
                                    self.Fk(1+l+(num_plato - 1)*(2.*self.num_sust + 1)) = self.Fk(1+l+(num_plato - 1)*(2.*self.num_sust + 1)) + self.entradas{m}.*self.concfeed(m/3, l);
                                end
                            end
                        end
                        if num_plato == 1
                            T = self.varX(1+(num_plato-1)*(2*(self.num_sust)+1)+2*self.num_sust);   % Temperatura Alternativamente self.platos(num_plato).T
                            P = self.platos(num_plato).P;
                            mezclaY = Mezcla(self.sust, self.platos(num_plato).y_i, self.alimentaciones{2}.mezcla.kij);
                            HgiV = 0;
                            Href = zeros(1, self.num_sust);
                            for i = 1:self.num_sust
                                Href(i) = self.sust(i).href;
                                try 
                                    cp = self.sust(i).cp_gi{1};
                                catch ME
                                    error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                                end
                                deltaHL = integral(@(t) cp(t), 273.15, T);
                                deltaHV = deltaHL;
                                HgiV = HgiV + deltaHV*mezclaY.conc(i) + Href(i)*mezclaY.conc(i);
                            end
                            Hdep_dewY = self.MEdE.entalpia(T, P, mezclaY, 'vap');
                            Hdep_refY = self.MEdE.entalpia(273.15, 101.325, mezclaY, 'liq');
                            HVj(num_plato) = HgiV - Hdep_dewY + Hdep_refY;  % HVj
                        end
                        mezclaX = Mezcla(self.sust, self.platos(num_plato).x_i, self.alimentaciones{2}.mezcla.kij);
                            HgiL = 0;
                            Href = zeros(1, self.num_sust);

                            for i = 1:self.num_sust
                                Href(i) = self.sust(i).href;
                                try 
                                    cp = self.sust(i).cp_gi{1};
                                catch ME
                                    error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                                end
                                deltaHL = integral(@(t) cp(t), 273.15, T);
                                HgiL =  HgiL + deltaHL*mezclaX.conc(i) + Href(i)*mezclaX.conc(i);
                            end
                            Hdep_bubX = self.MEdE.entalpia(T, P, mezclaX, 'liq');
                            Hdep_refX = self.MEdE.entalpia(273.15, 101.325, mezclaX, 'liq');
                            HLj(num_plato) = HgiL - Hdep_bubX + Hdep_refX; % HLj 
                        if num_plato < self.etapas
                            TVp1 = self.varX(1+(num_plato)*(2*(self.num_sust)+1)+2*self.num_sust);   % Temperatura Alternativamente self.platos(num_plato).T
                            PVp1 = self.platos(num_plato+1).P;
                            Ycj = self.ysubj(num_plato + 1,:);
                            mezclaYp1 = Mezcla(self.sust, Ycj, self.alimentaciones{2}.mezcla.kij);
                            if num_plato > 1
                                TLm1 =  self.varX(1+(num_plato-2)*(2*(self.num_sust)+1)+2*self.num_sust);
                                PLm1 = self.platos(num_plato-1).P;
                                Xcjm1 = self.xsubj(num_plato - 1,:);
                                mezclaXm1 = Mezcla(self.sust, Xcjm1, self.alimentaciones{2}.mezcla.kij);
                            end
                            HgiVp1 =0;
                            HgiL = 0;
                            HgiV = 0;
                            Href = zeros(1, self.num_sust);
                            for i = 1:self.num_sust
                                Href(i) = self.sust(i).href;
                                try 
                                    cp = self.sust(i).cp_gi{1};
                                catch ME
                                    error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                                end
                                deltaHV = integral(@(t) cp(t), 273.15, TVp1);
                                HgiVp1 = HgiVp1 + deltaHV*mezclaYp1.conc(i) + Href(i)*mezclaYp1.conc(i);
                            end
                            Hdep_dewY = self.MEdE.entalpia(TVp1, PVp1, mezclaYp1, 'vap');
                            Hdep_refY = self.MEdE.entalpia(273.15, 101.325, mezclaYp1, 'liq');
                            HVj(num_plato+1) = HgiVp1 - Hdep_dewY + Hdep_refY;  % HVj
                        end
                        if ~isempty(self.platos(num_plato).salidaV)
                            if num_plato ~=1
                                if self.platos(num_plato).salidaV > 1e-4
                                    Sj(num_plato) = self.platos(num_plato).salidaV/self.platos(num_plato).V;
                                else
                                    Sj(num_plato) = 0;
                                end
                            else
                                Sj(num_plato) = 0;
                            end
                        else
                            Sj(num_plato) = 0;
                        end
                        HgiLm1T = 0;
                        HgiLp1T = 0;

                        HgiVm1T = 0;
                        HgiVp1T = 0;
                        delta = 1e-8;


                        mezclaY.conc = self.platos(num_plato).y_i;
                        HgiVm1T = 0;
                        HgiVp1T = 0;
                        for i = 1:self.num_sust
                            try 
                                cp = self.sust(i).cp_gi{1};
                            catch ME
                                error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                            end
                            if num_plato > 1                            
                                deltaHLm1T = integral(@(t) cp(t), 273.15, T - delta);
                                deltaHLp1T = integral(@(t) cp(t), 273.15, T + delta);
                                HgiLm1T = HgiLm1T + deltaHLm1T*mezclaX.conc(i) + Href(i) * mezclaX.conc(i);
                                HgiLp1T = HgiLp1T + deltaHLp1T*mezclaX.conc(i) + Href(i) * mezclaX.conc(i);
                            else 
                                HgiLm1T = 0;
                                HgiLp1T = 0;
                            end
                            deltaHVm1T = integral(@(t) cp(t), 273.15, T - delta);
                            deltaHVp1T = integral(@(t) cp(t), 273.15, T + delta);

                            HgiVm1T = HgiVm1T + deltaHVm1T*mezclaY.conc(i) + Href(i) * mezclaY.conc(i);
                            HgiVp1T = HgiVp1T + deltaHVp1T*mezclaY.conc(i) + Href(i) * mezclaY.conc(i);
                        end
                        if num_plato > 1
                            Hdep_bubXm1 = self.MEdE.entalpia(T - delta, P, mezclaX, 'liq');
                            Hdep_bubXp1 = self.MEdE.entalpia(T + delta, P, mezclaX, 'liq');
                        else 
                            Hdep_bubXm1 =  0;
                            Hdep_bubXp1 = 0;
                        end
                            Hdep_dewYm1 = self.MEdE.entalpia(T - delta, P, mezclaY, 'vap');
                            Hdep_dewYp1 = self.MEdE.entalpia(T + delta, P, mezclaY, 'vap');
                            dHVjdT = (HgiVp1T - HgiVm1T)/(2*delta) -(Hdep_dewYp1 - Hdep_dewYm1)/(2*delta);

                        dHLjdT = (HgiLp1T - HgiLm1T)/(2*delta) - (Hdep_bubXp1 - Hdep_bubXm1)/(2*delta);
                        HgiLm1m1T = 0;
                        HgiLm1p1T = 0;       % Hgasideal para calculo L etapa j-1 -delta
                        HgiVp1m1T = 0;
                        HgiVp1p1T = 0;
                        HgiLp1T = 0;
                        HgiLm1T = 0;
                        for i = 1:self.num_sust
                            try 
                                cp = self.sust(i).cp_gi{1};
                            catch ME
                                error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                            end
                            if num_plato > 1
                                deltaHVm1T = integral(@(t) cp(t), 273.15, TLm1 - delta);
                                deltaHVp1T = integral(@(t) cp(t), 273.15, TLm1 + delta);

                                HgiLm1T = HgiLm1T +  deltaHVm1T*mezclaXm1.conc(i)+ Href(i) * mezclaXm1.conc(i);
                                HgiLp1T = HgiLp1T + deltaHVp1T*mezclaXm1.conc(i)+ Href(i) * mezclaXm1.conc(i);
                            else
                                HgiLm1T = 0;
                                HgiLp1T = 0;
                            end
                        end

                        if num_plato > 1
                            Hdep_bubXm1 = self.MEdE.entalpia(TLm1 - delta, PLm1, mezclaXm1, 'liq');
                            Hdep_bubXp1 = self.MEdE.entalpia(TLm1 + delta, PLm1, mezclaXm1, 'liq');
                        else
                            Hdep_bubXp1 = 0;
                            Hdep_bubXm1 = 0;
                        end
                        dHLm1dT = (HgiLp1T - HgiLm1T)/(2*delta) - (Hdep_bubXp1 - Hdep_bubXm1)/(2*delta);   %Diferencial de HL etapa j-1 
                        HgiVp1m1T = 0;
                        HgiVp1p1T = 0;
                        if num_plato < self.etapas
                            for i = 1:self.num_sust
                                try 
                                    cp = self.sust(i).cp_gi{1};
                                catch ME
                                    error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                                end
                                deltaHVp1m1T = integral(@(t) cp(t), 273.15, TVp1 - delta);
                                deltaHVp1p1T = integral(@(t) cp(t), 273.15, TVp1 + delta);
                                HgiVp1m1T = HgiVp1m1T + deltaHVp1m1T*mezclaYp1.conc(i);
                                HgiVp1p1T = HgiVp1p1T + deltaHVp1p1T*mezclaYp1.conc(i);
                            end
                        else
                            HgiVp1m1T = 0;
                            HgiVp1p1T = 0;
                        end
                        if num_plato < self.etapas
                            Hdep_dewYp1 = self.MEdE.entalpia(TVp1 + delta, PVp1, mezclaYp1, 'vap');
                            Hdep_dewYm1 = self.MEdE.entalpia(TVp1 - delta, PVp1, mezclaYp1, 'vap');
                        else
                            Hdep_dewYp1 = 0;
                            Hdep_dewYm1 = 0;
                        end
                        dHVp1dT = (HgiVp1p1T - HgiVp1m1T)/(2*delta) - (Hdep_dewYp1 - Hdep_dewYm1)/(2*delta);   %Diferencial de HV etapa j+1

                        dQjdTjp1 = 0;
                        dQjdTjm1 = 0;
                        if num_plato > 1 && num_plato < self.etapas
                            self.Aij(1, (num_plato - 2)*((2*self.num_sust + 1))+2.*self.num_sust+1) =  -sum(self.platos(num_plato-1).l_i(:))*(dHLm1dT) - dQjdTjm1;
                        end
                        if num_plato < self.etapas && num_plato > 1
                            self.Cij(1, (num_plato - 1)*((2*self.num_sust + 1))+2.*self.num_sust+1) =  -sum(self.platos(num_plato+1).v_i(:))*(dHVp1dT) - dQjdTjp1;
                        end
                        for dFidli = 1:self.num_sust
                            %Funcion Fj siendo Hj
                            if num_plato < self.etapas && num_plato > 1
                                self.Bij(1, self.num_sust + (num_plato - 1)*(2*self.num_sust + 1)+dFidli) = HLj(num_plato)*(1+sj(num_plato));
                            end
                            if num_plato > 1 && num_plato < self.etapas
                                self.Aij(1, self.num_sust + (num_plato-2)*(2*self.num_sust + 1)+dFidli) = -HLj(num_plato - 1);
                            end
                        end

                        if num_plato ~= 1  && num_plato ~= self.etapas
                            self.Bij(1, (num_plato - 1)*((2*self.num_sust + 1))+2*self.num_sust+1) = dHLjdT*(1+sj(num_plato))*sum(self.platos(num_plato).l_i(:)) + dHVjdT.*(1+Sj(num_plato))*sum(self.platos(num_plato).v_i(:));
                        end
                        for dFidvi = 1:self.num_sust
                            %Funcion Fj siendo Hj
                            if num_plato ~= 1 && num_plato ~= self.etapas
                                self.Bij(1, (num_plato-1)*(2*self.num_sust + 1) + dFidvi) = HVj(num_plato)*(1+Sj(num_plato));
                            end
                            if num_plato < self.etapas && num_plato > self.etapas
                                self.Cij(1, (num_plato-1)*(2*self.num_sust + 1)+dFidvi) = -HVj(num_plato + 1);
                            end

                        end
                        mezclaX = mezclaY;
                        mezclaX.conc = self.platos(num_plato).x_i;
                        %por diferenciacion numerica de 3 puntos centrada de Ej = 0 = Kij*lij*(sum(vj))/(sum(lj))
                        % siendo la funci�n de la temperatura solo dependiente de Kij
                        % Kijp1(Tj) = Kij(T0+dif) = Kij(T0 + 1E-7)
                        fGijp1 = self.MEdE.fugF(T+delta, P, mezclaY, 'vap');
                        fLijp1 = self.MEdE.fugF(T+delta, P, mezclaX, 'liq');
                        fGijm1 = self.MEdE.fugF(T-delta, P, mezclaY, 'vap');
                        fLijm1 = self.MEdE.fugF(T-delta, P, mezclaX, 'liq');
                    %                 fGij = MEdE.fugF(Tj, P, mezcla, 'vap');
                    %                 fLij = MEdE.fugF(Tj, P, mezcla, 'liq');
                        Kijp1 = fLijp1 ./ fGijp1;
                        Kijm1 = fLijm1 ./ fGijm1;
                        fGij = self.MEdE.fugF(T, P, mezclaY, 'vap');
                        fLij = self.MEdE.fugF(T, P, mezclaX, 'liq');
                        self.perfil_k(:,num_plato) = fLij ./ fGij;
                        for dFidvi = 1:self.num_sust
                            %Funcion Fj siendo Ej
                            dEjdTj = ((Kijp1(dFidvi) - Kijm1(dFidvi))./(2*delta))*(self.platos(num_plato).l_i(dFidvi)*(sum(self.platos(num_plato).v_i(:)))/(sum(self.platos(num_plato).l_i(:))));
                            self.Bij(self.num_sust + 1 + dFidvi, (num_plato-1)*(2*self.num_sust + 1) + 2*self.num_sust + 1) = dEjdTj;

                            self.Bij(dFidvi + 1, (num_plato-1)*(2*self.num_sust + 1) + dFidvi) = 1+Sj(num_plato);

                            self.Bij(self.num_sust + 1 + dFidvi, (num_plato-1)*(2*self.num_sust + 1) + dFidvi) = self.platos(num_plato).K(dFidvi)*self.platos(num_plato).l_i(dFidvi)*(1)/(sum(self.platos(num_plato).l_i(:))) - 1;
                        end
                        for dFidvi = self.num_sust:2
                            for dFidvj = 1:1:self.num_sust
                                self.Bij(self.num_sust + dFidvi, (num_plato-1)*(2*self.num_sust + 1) + dFidvi) =  self.platos(num_plato).K(dFidvj)*self.platos(num_plato).l_i(dFidvj)*(1)/(sum(self.platos(num_plato).l_i(:)));
                            end
                        end
                        for dFidvi = 1:self.num_sust
                            for dFidvj = 1:1:self.num_sust
                                if dFidvj ~=  dFidvi 
                                    self.Bij(self.num_sust + 1 + dFidvj, (num_plato-1)*(2*self.num_sust + 1) + dFidvi) = self.platos(num_plato).K(dFidvj)*self.platos(num_plato).l_i(dFidvj)*(1)/(sum(self.platos(num_plato).l_i(:)));
                                end
                            end
                        end
                        for dFidli = 1:self.num_sust
                            %Funcion Fj siendo Mij
                            self.Bij(dFidli + 1, (num_plato-1)*(2*self.num_sust + 1) + dFidli + self.num_sust) = 1+sj(num_plato);
                            %fprintf(1, '%f\n', self.platos(num_plato).K(num_plato, dFidli)*((sum(self.platos(num_plato).l_i(num_plato, :))-self.platos(num_plato).l_i(num_plato, dFidli))*(sum(self.platos(num_plato).v_i(num_plato, :))))/(sum(self.platos(num_plato).l_i(num_plato, :))^2))
                            self.Bij(self.num_sust + 1 + dFidli, self.num_sust + (num_plato - 1)*(2*self.num_sust + 1)+dFidli) = self.platos(num_plato).K(dFidli)*((sum(self.platos(num_plato).l_i(:))-self.platos(num_plato).l_i(dFidli))*(sum(self.platos(num_plato).v_i(:))))/(sum(self.platos(num_plato).l_i(:))^2);
                        end
                        for dFidli = 1:self.num_sust
                            for dFidlj = self.num_sust:-1:1
                                %Funcion Fj siendo Ej
                                if dFidlj ~= dFidli
                                    self.Bij(1 + dFidlj + self.num_sust, self.num_sust + (num_plato-1)*(2*self.num_sust + 1) + dFidli) = -self.platos(num_plato).K(dFidlj)*((sum(self.platos(num_plato).v_i(:)))*(self.platos(num_plato).l_i(dFidlj)))/(sum(self.platos(num_plato).l_i(:))^2);
                                end
                            end
                        end
                        for i = 1:self.num_sust
                            self.Fk(1+i+(num_plato - 1)*(2.*self.num_sust + 1)) = self.Fk(1+i+(num_plato - 1)*(2.*self.num_sust + 1)) - self.platos(num_plato).l_i(i)*(1+sj(num_plato)) - self.platos(num_plato).v_i(i)*(1+Sj(num_plato));
                            if num_plato > 1
                                self.Fk(1+i+(num_plato - 1)*(2.*self.num_sust + 1)) =  self.Fk(1+i+(num_plato - 1)*(2.*self.num_sust + 1)) + self.platos(num_plato-1).l_i(i);
                            end
                            if num_plato < self.etapas
                                self.Fk(1+i+(num_plato - 1)*(2.*self.num_sust + 1)) =  self.Fk(1+i+(num_plato - 1)*(2.*self.num_sust + 1)) + self.platos(num_plato+1).v_i(i);
                            end      
                        end

                        for i = 1:self.num_sust
                            self.Fk(1+self.num_sust+i +(num_plato - 1)*(2.*self.num_sust + 1)) = self.Fk(1+self.num_sust+i +(num_plato - 1)*(2.*self.num_sust + 1)) - self.platos(num_plato).K(i)*self.platos(num_plato).l_i(i)*((sum(self.platos(num_plato).v_i(:)))./(sum(self.platos(num_plato).l_i(:))))+self.platos(num_plato).v_i(i);     
                        end
                        indicee = 1;
                        if  any(num_plato == self.platos_entradas(:))    %Entalpia de alimentaci�n HF
                            if length(self.platos_entradas) == 1
                                indicee = 1;
                            else
                                indicee = find((self.platos_entradas(:) == num_plato) ~= 0, 1, 'first');
                            end
                            HF(num_plato) = alimen(indicee).H;
                        else
                            HF(num_plato) = 0;
                        end
                        if num_plato == 1
                            if sj(num_plato) == 0
                                self.Fk(1) = -self.platos(1).L  + self.reflujo * self.platos(1).V;
                            else
                                self.Fk(1) = -sum(self.varX(2*self.num_sust + 1+1:2*self.num_sust + 1 + self.num_sust))  + self.reflujo * (sum(self.varX(1:self.num_sust)) + sj(num_plato)*(sum(self.varX(self.num_sust+1:2*self.num_sust)))) + sj(num_plato)*(sum(self.varX(self.num_sust+1:2*self.num_sust))) + sum(self.varX(1:self.num_sust));
                            end
                        elseif num_plato == self.etapas
                            self.Fk(1+(num_plato - 1)*(2.*self.num_sust + 1)) = -(sum(self.varX(1 + self.num_sust + (num_plato - 1)*(2*self.num_sust + 1):1+self.num_sust +(num_plato - 1)*(2*self.num_sust + 1)+self.num_sust-1)) - self.bflujo);
                        else 
                            self.Fk(1+(num_plato - 1)*(2.*self.num_sust + 1)) =  - HLj(num_plato)*(1+sj(num_plato))*(sum(self.platos(num_plato).l_i(:))) - HVj(num_plato)*(1+Sj(num_plato))*(sum(self.platos(num_plato).v_i(:))) + HF(num_plato)*alimen(indicee).molF + self.perfil_q(num_plato-1) + HLj(num_plato - 1) * sum(self.platos(num_plato-1).l_i(:)) + HVj(num_plato + 1) * sum(self.platos(num_plato + 1).v_i(:));
                        end
                        if num_plato > 1 && num_plato < self.etapas
                            self.Cij(1,1+(num_plato-1)*(2*self.num_sust + 1):1+(num_plato-1)*(2*self.num_sust + 1)+self.num_sust -1) = -HVj(num_plato + 1);
                        end
                        T = TVp1;
                        P = PVp1;
                        mezclaY = mezclaYp1;
                    end

                    [deltaXvar, self.matriz_tridiag ] = tridiagThomas(self.Bij, self.Aij, self.Cij, self.Fk);
                    condition = 0;
                    for num_plato = 1:self.etapas
                        nuevas_varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust) = self.varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust) + self.damping.*deltaXvar(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust);   
                        if any(nuevas_varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust-1) < 0)
                            nuevas_varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust - 1) = self.varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust-1).*exp((self.damping.*deltaXvar(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust-1))./(self.varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust-1)));
                        end
                        
                        if num_plato == self.etapas && condition == 1
                            condition = 0;
                        end
                        Tbackup(num_plato) = self.perfil_t(num_plato);
                        Vbackup(num_plato) = self.perfil_v(num_plato);
                        Lbackup(num_plato) = self.perfil_l(num_plato);
                        self.perfil_t(num_plato) = nuevas_varX(1+(num_plato-1)*(2*(self.num_sust)+1)+2*self.num_sust);
                        if sum(nuevas_varX(1+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + self.num_sust-1)) > self.tol
                            self.perfil_v(num_plato) = sum(nuevas_varX(1+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + self.num_sust-1));
                        end
                        self.perfil_l(num_plato) = sum(nuevas_varX(1+self.num_sust+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + 2*self.num_sust-1));
                        
                        if (sum(self.Fk.^2)) < respaldovalI
                            respaldovalI = (sum(self.Fk.^2));
                            respaldobi = self.bsubi;
                            respaldodi = self.dsubi;
                            respaldoV = self.perfil_v;
                            respaldoL = self.perfil_l;
                            respaldovi = self.perfil_vi;
                            respaldoli = self.perfil_li;
                            respaldot = self.perfil_t;
                            respaldok =  self.perfil_k;
                            respaldoplatos = self.platos;
                            respaldoxi =  self.xsubj;
                            respaldoyi = self.ysubj;
                            if isempty(self.salidas) || ~isa(self.salidas, 'cell')
                                self.qc = -self.perfil_v(1)*(HVj(1)) - HLj(1)*self.perfil_l(1) + self.perfil_v(2)*HVj(2);
                            else
                                if self.salidas{1} == 1 && self.salidas{2} == 1 
                                    self.qc = -self.perfil_v(1)*(HVj(1)) - HLj(1)*self.perfil_l(1) - sj(1)*self.perfil_l(1)*HLj(1) + self.perfil_v(2)*HVj(2);
                                elseif self.salidas{1} == 1 && self.salidas{2} == 0
                                    self.qc = -self.perfil_v(1)*(HVj(1)) - HLj(1)*self.perfil_l(1) - sj(1)*self.perfil_l(1)*HLj(1) + self.perfil_v(2)*HVj(2);
                                else
                                    self.qc = -self.perfil_v(1)*(HVj(1)) - HLj(1)*self.perfil_l(1) + self.perfil_v(2)*HVj(2);
                                end
                            end
                            Fh = 0;
                            Vh1 = 0;
                            Uh = 0; %Salidas Laterales Liquidas 1
                            Wh = 0; %Salidas Laterales Vapor 1
                            Lh = 0; %Liquido que va al plato 2
                            Vh = 0; %Pudiera haber vapor fuga del plato 1
                            Vh2 = 0; %Vapor que entra al plato 1 del plato 2
                            LhEND = 0;
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
                                if ~isempty(self.platos(iterx).salidaV)
                                    if iterx == 1
                                        Vh1 = Vh1 +  self.perfil_v(1) * HVj(iterx);  % que no se sume
                                            %2 veces el destilado vapor
                                    else    
                                        Wh = Wh + self.platos(iterx).salidaV * HVj(iterx);
                                    end
                                end
                                if ~isempty(self.platos(iterx).salidaL)
                                    if iterx == self.etapas
                                        LhEND = LhEND  + self.platos(iterx).salidaL * HLj(num_plato);
                                        %que no se sume dos veces el fondo l�quido
                                    else                            
                                        Uh = Uh + self.platos(iterx).salidaL * HLj(num_plato);
                                    end
                                end
                            end
                            self.qb = Fh - Uh  - Wh - Vh1 - LhEND - sum(self.perfil_q) - self.qc;
                            respaldoqc= self.qc;
                            respaldoqb = self.qb;
                            self.respaldoiter = respaldandoiter;
                        end
                        
                        valI = sum(self.Fk.^2);
    %                     valI = valI + ((self.perfil_t(num_plato) - Tbackup(num_plato))/self.perfil_t(num_plato))^2;
    %                     if self.perfil_v(num_plato) > 0.00001 .* self.tol
    %                         valI = valI + ((self.perfil_v(num_plato) - Vbackup(num_plato))/self.perfil_v(num_plato))^2;
    %                     end
    %                     if self.perfil_l(num_plato) > 0.00001.* self.tol
    %                         valI = valI + ((self.perfil_l(num_plato) - Lbackup(num_plato))/self.perfil_l(num_plato))^2;
    %                     end
                        self.varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust) = nuevas_varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust);

                        self.perfil_vi(num_plato, :) = self.varX(1+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + self.num_sust-1)';
                        self.perfil_li(num_plato, :) = self.varX(1+self.num_sust+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + 2*self.num_sust-1)';
                        if num_plato == 1
                            if ~isempty(self.salidas) && isa(self.salidas, 'cell') 
                                if self.salidas{1} == 1 && self.salidas{2} == 1 && self.perfil_v(num_plato)/(self.perfil_v(num_plato)+self.perfil_l(num_plato)) < 1e-4
                                    self.ysubj(num_plato,:) = self.varX(1+(num_plato)*(2*(self.num_sust)+1):1+(num_plato)*(2*(self.num_sust)+1) + self.num_sust-1)'./sum(self.varX(1+(num_plato)*(2*(self.num_sust)+1):1+(num_plato)*(2*(self.num_sust)+1) + self.num_sust-1));
                                else
                                    self.ysubj(num_plato,:) = self.varX(1:self.num_sust)'./sum(self.varX(1:self.num_sust));
                                end
                            else
                                self.ysubj(num_plato,:) = self.varX(1:self.num_sust)'./sum(self.varX(1:self.num_sust));
                            end
                        else
                            self.ysubj(num_plato,:) = self.varX(1+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + self.num_sust-1)'./sum(self.varX(1+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + self.num_sust-1));
                        end
                        self.xsubj(num_plato,:) = self.varX(1+self.num_sust+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + 2*self.num_sust-1)'./sum(self.varX(1+self.num_sust + (num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + 2*self.num_sust-1));                    
                        self.ysubj(num_plato,:) = self.ysubj(num_plato,:);
                        self.xsubj(num_plato,:) = self.xsubj(num_plato,:);
                        self.perfil_t(num_plato) = self.varX(1+(num_plato-1)*(2*(self.num_sust)+1)+2*self.num_sust);

                    end

                    for iiiiter = 1:self.etapas
                        if iiiiter == self.etapas
                            fprintf(1, ' %f \n', self.perfil_t(iiiiter));
                            break
                        end
                        if iiiiter == 1
                            fprintf(1, '\n %f ', self.perfil_t(iiiiter));
                        else
                            fprintf(1, ' %f ', self.perfil_t(iiiiter));
                        end

                    end

                    for iteru = 1:self.etapas

                        self.platos(iteru).setV(self.perfil_v(iteru));
                        self.platos(iteru).setL(self.perfil_l(iteru));
                        self.platos(iteru).setyi(self.ysubj(iteru,:));                    
                        self.platos(iteru).setxi(self.xsubj(iteru,:));
                        self.platos(iteru).K = self.perfil_k(:, iteru);
                    end
                    self.bsubi = self.xsubj(end,:).*(self.perfil_l(end));
    %                 self.perfil_t(1) = Backupt1;

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


                    for i = 1:self.etapas
                        self.platos(i).L = self.perfil_l(i);
                        self.platos(i).V = self.perfil_v(i);
                        self.platos(i).y_i = self.ysubj(i,:);
                        self.platos(i).x_i = self.xsubj(i,:);
                        self.platos(i).v_i = self.perfil_vi(i,:);
                        self.platos(i).l_i = self.perfil_li(i,:);
                    end
                    self.bflujo = sum(self.perfil_li(end,:));
                    if isempty(self.salidas) 
                        self.dflujo = sum(self.perfil_vi(1,:)); 
                    else
                        if isa(self.salidas, 'cell')
                            if self.salidas{1} == 1 && self.salidas{2} == 1
                                self.dflujo =  sum(self.perfil_vi(1,:)) + self.F - self.UmW - sum(self.perfil_li(end,:));
                            end
                        else
                            self.dflujo = sum(self.perfil_vi(1,:)); 
                        end
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
                    self.respaldovalI(self.actualiter) = valI;

                    if valI < self.tol
                        self.laststep = logical(1);
                    end
                    if self.perfil_t(1) < 0 
                        self.perfil_t(1) = Tbackup(1);
                        self.varX(2.*self.num_sust + 1) = Tbackup(1);
                    end
                end

                if isempty(self.salidas) || ~isa(self.salidas, 'cell')
                    mezclaX = Mezcla(self.sust, self.xsubj(1,:), self.alimentaciones{2}.mezcla.kij);
                    [~, yy ] = self.MEdE.BubbleT(mezclaX, self.perfil_p(1), self.alimenciones{2}.mezcla.kij);
                    [Tb, ~ ] = self.MEdE.BubbleT(mezclaX, self.perfil_p(1), self.alimenciones{2}.mezcla.kij, yy);
                    mezclaY = mezclaX;
                    mezclaY.conc = self.ysubj(1,:);
                    for ittwee = 1:self.num_sust
                        try 
                            cp = self.sust(i).cp_gi{1};
                        catch ME
                            error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                        end
                        deltaHV1 = integral(@(t) cp(t), 273.15, Tb);
                        HgiV1 = deltaHV1*mezclaY.conc(i) + Href(ittwee) * mezclaY.conc(i);                
                    end
                    HVj(1) = HgiV1 - self.MEdE.entalpia( Tb, self.perfil_p(1), mezclaY, 'vap') + self.MEdE.entalpia( 273.15, 101.325, mezclaY, 'liq');
                    self.qc = -self.perfil_v(1)*(HVj(1)) - HLj(1)*self.perfil_l(1) + self.perfil_v(2)*HVj(2);
                else
                    if self.salidas{1} == 1 && self.salidas{2} == 1 
                        self.qc = -self.perfil_v(1)*(HVj(1)) - HLj(1)*self.perfil_l(1) - sj(1)*self.perfil_l(1)*HLj(1) + self.perfil_v(2)*HVj(2);
                    elseif self.salidas{1} == 1 && self.salidas{2} == 0
                        mezclaX = Mezcla(self.sust, self.xsubj(1,:), self.alimentaciones{2}.mezcla.kij);
                        [~, yy ] = self.MEdE.BubbleT(mezclaX, self.perfil_p(1), self.alimenciones{2}.mezcla.kij);
                        [Tb, ~ ] = self.MEdE.BubbleT(mezclaX, self.perfil_p(1), self.alimenciones{2}.mezcla.kij, yy);
                        mezclaY = mezclaX;
                        mezclaY.conc = self.ysubj(1,:);
                        for ittwee = 1:self.num_sust
                            try 
                                cp = self.sust(i).cp_gi{1};
                            catch ME
                                error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                            end
                            deltaHV1 = integral(@(t) cp(t), 273.15, Tb);
                            HgiV1 = deltaHV1*mezclaY.conc(i) + Href(ittwee) * mezclaY.conc(i);                
                        end
                        HVj(1) = HgiV1 - self.MEdE.entalpia( Tb, self.perfil_p(1), mezclaY, 'vap') + self.MEdE.entalpia( 273.15, 101.325, mezclaY, 'liq');
                        self.qc = -self.perfil_v(1)*(HVj(1)) - HLj(1)*self.perfil_l(1) - sj(1)*self.perfil_l(1)*HLj(1) + self.perfil_v(2)*HVj(2);
                    else
                         mezclaX = Mezcla(self.sust, self.xsubj(1,:), self.alimentaciones{2}.mezcla.kij);
                        [~, yy ] = self.MEdE.BubbleT(mezclaX, self.perfil_p(1), self.alimenciones{2}.mezcla.kij);
                        [Tb, ~ ] = self.MEdE.BubbleT(mezclaX, self.perfil_p(1), self.alimenciones{2}.mezcla.kij, yy);
                        mezclaY = mezclaX;
                        mezclaY.conc = self.ysubj(1,:);
                        for ittwee = 1:self.num_sust
                            try 
                                cp = self.sust(i).cp_gi{1};
                            catch ME
                                error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                            end
                            deltaHV1 = integral(@(t) cp(t), 273.15, Tb);
                            HgiV1 = deltaHV1*mezclaY.conc(i) + Href(ittwee) * mezclaY.conc(i);                
                        end
                        HVj(1) = HgiV1 - self.MEdE.entalpia( Tb, self.perfil_p(1), mezclaY, 'vap') + self.MEdE.entalpia( 273.15, 101.325, mezclaY, 'liq');
                        self.qc = -self.perfil_v(1)*(HVj(1)) - HLj(1)*self.perfil_l(1) + self.perfil_v(2)*HVj(2);
                    end

                end
                Fh = 0;
                Vh1 = 0;
                Uh = 0; %Salidas Laterales Liquidas 1
                Wh = 0; %Salidas Laterales Vapor 1
                Lh = 0; %Liquido que va al plato 2
                Vh = 0; %Pudiera haber vapor fuga del plato 1
                Vh2 = 0; %Vapor que entra al plato 1 del plato 2
                LhEND = 0;
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
                    if ~isempty(self.platos(iterx).salidaV)
                        if iterx == 1
                            Vh1 = Vh1 +  self.perfil_v(1) * HVj(iterx);  % que no se sume
                                %2 veces el destilado vapor
                        else    
                            Wh = Wh + self.platos(iterx).salidaV * HVj(iterx);
                        end
                    end
                    if ~isempty(self.platos(iterx).salidaL)
                        if iterx == self.etapas
                            LhEND = LhEND  + self.platos(iterx).salidaL * HLj(num_plato);
                            %que no se sume dos veces el fondo l�quido
                        else                            
                            Uh = Uh + self.platos(iterx).salidaL * HLj(num_plato);
                        end
                    end
                end
                self.qb = Fh - Uh  - Wh - Vh1 - LhEND - sum(self.perfil_q) - self.qc;
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
                self.bsubi = respaldobi;
                self.dsubi = respaldodi;
                self.dflujo = sum(respaldodi);
                self.bflujo = sum(respaldobi);
                fprintf('Fallo la convergencia en la iteracion = %i \n', iter);
                fprintf('El mejor resultado obtenido fue en la iteracion = %i \n', self.respaldoiter)
                display(self);
                fprintf('El mejor error obtenido fue de = %f \n', respaldovalI)
            end
        end
        function newtonraphson3(self)  %Especificaciones R y B
            % Se construye la matriz tridiagonal de derivadas
            %  [ B(1)  C(1)                                  ] [  x(1)  ]   [  F(1)  ]
            %  [ A(2)  B(2)  C(2)                            ] [  x(2)  ]   [  F(2)  ]
            %  [       A(3)  B(3)  C(3)                      ] [        ]   [        ]
            %  [            ...   ...   ...                  ] [  ...   ] = [  ...   ]
            %  [                    ...    ...    ...        ] [        ]   [        ]
            %  [                        A(n-1) B(n-1) C(n-1) ] [ x(n-1) ]   [ F(n-1) ]
            %  [                                 A(n)  B(n)  ] [  x(n)  ]   [  F(n)  ]
            % en ella, las matrices B(j) son la contribuci�n al plato j
            % de las variables de estado de los platos j-1, j y j+1
            % Las matrices A(j-1) son la dependencia con el plato j-1 del factor
            % de cambio de las variables de estado del plato j
            % Las matrices C(j+1) son la dependencia con el plato j+1 del factor
            % de cambio de las variables de estado del plato j
            %
            %La ecuaci�n de reemplazo del balance de energ�a en la etapa 1
            %es la OD: (Henley & Seader, 2008)
            %
            % OD = vi,1 - di,1 = 0
            %
            % La ecuaci�n de reemplazo del balance de energ�a en la etapa N
            % es la OW
            %
            % li,N - wi = 0
            %
            % B(j), A(j-1), C(j+1)  son matrices de la forma
            %                 1                          2*n+1
            % dFj/dX        X = v1   � X = lij          X = Tj
            % 1     F = ODj [dODj/dv1  � dODj/dlij  �  � dODj/dTj]
            % 2     F = M1j [dM1j/dv1 � dM1j/dlij �  � dM1j/dTj ]
            % �             [   �        dFj/dX     �      �    ] 
            % 2*n+1         [   �          �       �       �   ]
%             if isempty(self.perfil_k) || isempty(self.perfil_v) || isempty(self.perfil_t)
%                 self = self.balanmasa();
%                 self.generar_etapas();
%             end
            if isempty(self.perfil_k) || isempty(self.perfil_vi) || isempty(self.perfil_t) || isempty(self.perfil_li)
                self = self.balanmasa();
            else
                self.generar_etapas();
            end
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
            self.perfil_v(1) = self.dflujo;
            if ~isempty(self.salidas) && isa(self.salidas, 'cell') 
                if self.salidas{1} == 1 && self.salidas{2} == 1
                    self.dflujo = self.dflujo + self.salidas{3};
                elseif self.salidas{1} == 1 && self.salidas{2} == 0
                    self.dflujo = self.dflujo + self.salidas{3};
                end
            end
            if isempty(self.tol)
                self.tol = self.convergence;
            end
            if isempty(self.damping) || ~isa(self.damping, 'double') 
                self.damping = 1;
            elseif isa(self.damping, 'double') && (self.damping < 0  || self.damping > 2) 
                self.damping = 1;
            end
            nuevas_varX = zeros((2.*self.num_sust + 1)*self.etapas, 1);
            self.varX = nuevas_varX;
            Vbackup = zeros(1, self.etapas);
            Tbackup = zeros(1, self.etapas);
            Lbackup = zeros(1, self.etapas);
            for iitter0 = 1 : self.etapas
                self.varX(1+(iitter0-1)*(2*(self.num_sust)+1):1+(iitter0-1)*(2*(self.num_sust)+1) + self.num_sust - 1) = self.perfil_vi(iitter0,:);
                self.varX(1+self.num_sust+(iitter0-1)*(2*(self.num_sust)+1):1+(iitter0-1)*(2*(self.num_sust)+1) + 2*self.num_sust-1) = self.perfil_li(iitter0,:);
                self.varX(1+(iitter0-1)*(2*(self.num_sust)+1)+2*self.num_sust) = self.perfil_t(iitter0);
            end
            if isempty(self.respaldovalI)
                valI = 1e308;
            else
                valI = self.respaldovalI;
            end
            iter = 0;
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
                respaldodi = self.dsubi;
                respaldobi = self.bsubi;
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
                    valI = 0;
                    iter = iter + 1;
                    if iter == 1
                        for i = 1:self.etapas
                            self.platos(i).L = self.perfil_l(i);
                            self.platos(i).V = self.perfil_v(i);
                            self.platos(i).y_i = self.ysubj(i,:);
                            self.platos(i).x_i = self.xsubj(i,:);
                            self.platos(i).v_i = self.perfil_vi(i,:);
                            self.platos(i).l_i = self.perfil_li(i,:);
                        end
                    end
                    self.Bij = zeros(2*self.num_sust +1, (self.etapas)*(2*self.num_sust + 1));
                    self.Aij = zeros(2*self.num_sust +1, (self.etapas-1)*(2*self.num_sust + 1));
                    self.Cij = zeros(2*self.num_sust +1, (self.etapas-1)*(2*self.num_sust + 1));
                    self.Fk = zeros((2*self.num_sust+1).*self.etapas, 1);
                    tamano = size(self.entradas);
                    HVj = zeros(1, self.etapas);
                    HLj = zeros(1, self.etapas);
                    sj = zeros(1, self.etapas);
                    Sj = zeros(1, self.etapas);
                    for etapa = 1:self.etapas 
                        if etapa < self.etapas
                            for dMidli = 1:self.num_sust
                                self.Cij(dMidli+1, (etapa-1)*(2*self.num_sust + 1)+dMidli) = -1;
                            end
                        end
                        if etapa > 1
                            for dMidli = 1:self.num_sust
                                self.Aij(dMidli+1, self.num_sust + (etapa-2)*(2*self.num_sust + 1)+dMidli) = -1;
                            end
                        end
                    end
                    self.platos_entradas = zeros(1, length(self.entradas)/3);
                    alimen = Corriente.empty(0, length(self.entradas)/3 );
                    for itere=1:3:length(self.entradas)
                        self.platos_entradas((itere-1)/3+1) = self.entradas{itere};
                    end
                    for itere = 1:2:length(self.alimentaciones)
                        alimen((itere-1)/2+1) = self.alimentaciones{itere + 1};
                    end

                    for num_plato = 1:self.etapas
                        if ~isempty(self.platos(num_plato).salidaL);
                            if num_plato ~= self.etapas
                                if self.platos(num_plato).salidaL > 1e-4
                                    sj(num_plato) = self.platos(num_plato).salidaL/self.platos(num_plato).L;
                                else
                                    sj(num_plato) = 0;
                                end
                            else
                                sj(num_plato) = 0;
                            end
                        else
                            sj(num_plato) = 0;
                        end
                        % La derivada de la funci�n de reemplazo OD

                        if num_plato == 1
                            self.Bij(1, (self.lk)) = 1;
                        end
                        % La derivada de la funci�n de reemplazo OW es siempre +1
                        % para todos los Li y -reflujo (osea -R) para todos los vi
                        if num_plato == self.etapas
                            self.Bij(1, (num_plato-1)*(2*self.num_sust+1) + self.num_sust + self.hk) = 1;
                        end
                        for l = 1:self.num_sust
                            for m = 3:3:tamano(2)
                                if num_plato == self.entradas{m-2}
                                    self.Fk(1+l+(num_plato - 1)*(2.*self.num_sust + 1)) = self.Fk(1+l+(num_plato - 1)*(2.*self.num_sust + 1)) + self.entradas{m}.*self.concfeed(m/3, l);
                                end
                            end
                        end
                        if num_plato == 1
                            T = self.varX(1+(num_plato-1)*(2*(self.num_sust)+1)+2*self.num_sust);   % Temperatura Alternativamente self.platos(num_plato).T
                            P = self.platos(num_plato).P;
                            mezclaY = Mezcla(self.sust, self.platos(num_plato).y_i, self.alimentaciones{2}.mezcla.kij);
                            HgiV = 0;
                            Href = zeros(1, self.num_sust);
                            for i = 1:self.num_sust
                                Href(i) = self.sust(i).href;
                                try 
                                    cp = self.sust(i).cp_gi{1};
                                catch ME
                                    error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                                end
                                deltaHL = integral(@(t) cp(t), 273.15, T);
                                deltaHV = deltaHL;
                                HgiV = HgiV + deltaHV*mezclaY.conc(i) + Href(i)*mezclaY.conc(i);
                            end
                            Hdep_dewY = self.MEdE.entalpia(T, P, mezclaY, 'vap');
                            Hdep_refY = self.MEdE.entalpia(273.15, 101.325, mezclaY, 'liq');
                            HVj(num_plato) = HgiV - Hdep_dewY + Hdep_refY;  % HVj
                        end
                        mezclaX = Mezcla(self.sust, self.platos(num_plato).x_i, self.alimentaciones{2}.mezcla.kij);
                            HgiL = 0;
                            Href = zeros(1, self.num_sust);
                        for i = 1:self.num_sust
                                Href(i) = self.sust(i).href;
                                try 
                                    cp = self.sust(i).cp_gi{1};
                                catch ME
                                    error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                                end
                                deltaHL = integral(@(t) cp(t), 273.15, T);
                                HgiL =  HgiL + deltaHL*mezclaX.conc(i) + Href(i)*mezclaX.conc(i);
                        end
                            Hdep_bubX = self.MEdE.entalpia(T, P, mezclaX, 'liq');
                            Hdep_refX = self.MEdE.entalpia(273.15, 101.325, mezclaX, 'liq');
                            HLj(num_plato) = HgiL - Hdep_bubX + Hdep_refX; % HLj 
                        if num_plato < self.etapas
                            TVp1 = self.varX(1+(num_plato)*(2*(self.num_sust)+1)+2*self.num_sust);   % Temperatura Alternativamente self.platos(num_plato).T
                            PVp1 = self.platos(num_plato+1).P;
                            Ycj = self.ysubj(num_plato + 1,:);
                            mezclaYp1 = Mezcla(self.sust, Ycj, self.alimentaciones{2}.mezcla.kij);
                            if num_plato > 1
                                TLm1 =  self.varX(1+(num_plato-2)*(2*(self.num_sust)+1)+2*self.num_sust);
                                PLm1 = self.platos(num_plato-1).P;
                                Xcjm1 = self.xsubj(num_plato - 1,:);
                                mezclaXm1 = Mezcla(self.sust, Xcjm1, self.alimentaciones{2}.mezcla.kij);
                            end
                            HgiVp1 =0;
                            HgiL = 0;
                            HgiV = 0;
                            Href = zeros(1, self.num_sust);
                            for i = 1:self.num_sust
                                Href(i) = self.sust(i).href;
                                try 
                                    cp = self.sust(i).cp_gi{1};
                                catch ME
                                    error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                                end
                                deltaHV = integral(@(t) cp(t), 273.15, TVp1);
                                HgiVp1 = HgiVp1 + deltaHV*mezclaYp1.conc(i) + Href(i)*mezclaYp1.conc(i);
                            end
                            Hdep_dewY = self.MEdE.entalpia(TVp1, PVp1, mezclaYp1, 'vap');
                            Hdep_refY = self.MEdE.entalpia(273.15, 101.325, mezclaYp1, 'liq');
                            HVj(num_plato+1) = HgiVp1 - Hdep_dewY + Hdep_refY;  % HVj
                        end
                        if ~isempty(self.platos(num_plato).salidaV)
                            if num_plato ~=1
                                if self.platos(num_plato).salidaV > 1e-4
                                    Sj(num_plato) = self.platos(num_plato).salidaV/self.platos(num_plato).V;
                                else
                                    Sj(num_plato) = 0;
                                end
                            else
                                Sj(num_plato) = 0;
                            end
                        else
                            Sj(num_plato) = 0;
                        end
                     if num_plato == 1
                            T = self.varX(1+(num_plato-1)*(2*(self.num_sust)+1)+2*self.num_sust);   % Temperatura Alternativamente self.platos(num_plato).T
                            P = self.platos(num_plato).P;
                            mezclaY = Mezcla(self.sust, self.platos(num_plato).y_i, self.alimentaciones{2}.mezcla.kij);
                            HgiV = 0;
                            Href = zeros(1, self.num_sust);
                            for i = 1:self.num_sust
                                Href(i) = self.sust(i).href;
                                try 
                                    cp = self.sust(i).cp_gi{1};
                                catch ME
                                    error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                                end
                                deltaHL = integral(@(t) cp(t), 273.15, T);
                                deltaHV = deltaHL;
                                HgiV = HgiV + deltaHV*mezclaY.conc(i) + Href(i)*mezclaY.conc(i);
                            end
                            Hdep_dewY = self.MEdE.entalpia(T, P, mezclaY, 'vap');
                            Hdep_refY = self.MEdE.entalpia(273.15, 101.325, mezclaY, 'liq');
                            HVj(num_plato) = HgiV - Hdep_dewY + Hdep_refY;  % HVj
                        end
                        mezclaX = Mezcla(self.sust, self.platos(num_plato).x_i, self.alimentaciones{2}.mezcla.kij);
                            HgiL = 0;
                            Href = zeros(1, self.num_sust);
                        for i = 1:self.num_sust
                                Href(i) = self.sust(i).href;
                                try 
                                    cp = self.sust(i).cp_gi{1};
                                catch ME
                                    error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                                end
                                deltaHL = integral(@(t) cp(t), 273.15, T);
                                HgiL =  HgiL + deltaHL*mezclaX.conc(i) + Href(i)*mezclaX.conc(i);
                        end
                            Hdep_bubX = self.MEdE.entalpia(T, P, mezclaX, 'liq');
                            Hdep_refX = self.MEdE.entalpia(273.15, 101.325, mezclaX, 'liq');
                            HLj(num_plato) = HgiL - Hdep_bubX + Hdep_refX; % HLj 
                        if num_plato < self.etapas
                            TVp1 = self.varX(1+(num_plato)*(2*(self.num_sust)+1)+2*self.num_sust);   % Temperatura Alternativamente self.platos(num_plato).T
                            PVp1 = self.platos(num_plato+1).P;
                            Ycj = self.ysubj(num_plato + 1,:);
                            mezclaYp1 = Mezcla(self.sust, Ycj, self.alimentaciones{2}.mezcla.kij);
                            if num_plato > 1
                                TLm1 =  self.varX(1+(num_plato-2)*(2*(self.num_sust)+1)+2*self.num_sust);
                                PLm1 = self.platos(num_plato-1).P;
                                Xcjm1 = self.xsubj(num_plato - 1,:);
                                mezclaXm1 = Mezcla(self.sust, Xcjm1, self.alimentaciones{2}.mezcla.kij);
                            end
                            HgiVp1 =0;
                            HgiL = 0;
                            HgiV = 0;
                            Href = zeros(1, self.num_sust);
                            for i = 1:self.num_sust
                                Href(i) = self.sust(i).href;
                                try 
                                    cp = self.sust(i).cp_gi{1};
                                catch ME
                                    error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                                end
                                deltaHV = integral(@(t) cp(t), 273.15, TVp1);
                                HgiVp1 = HgiVp1 + deltaHV*mezclaYp1.conc(i) + Href(i)*mezclaYp1.conc(i);
                            end
                            Hdep_dewY = self.MEdE.entalpia(TVp1, PVp1, mezclaYp1, 'vap');
                            Hdep_refY = self.MEdE.entalpia(273.15, 101.325, mezclaYp1, 'liq');
                            HVj(num_plato+1) = HgiVp1 - Hdep_dewY + Hdep_refY;  % HVj
                        end
                        if ~isempty(self.platos(num_plato).salidaV)
                            if num_plato ~=1
                                if self.platos(num_plato).salidaV > 1e-4
                                    Sj(num_plato) = self.platos(num_plato).salidaV/self.platos(num_plato).V;
                                else
                                    Sj(num_plato) = 0;
                                end
                            else
                                Sj(num_plato) = 0;
                            end
                        else
                            Sj(num_plato) = 0;
                        end
                        HgiLm1T = 0;
                        HgiLp1T = 0;

                        HgiVm1T = 0;
                        HgiVp1T = 0;
                        delta = 1e-8;


                        mezclaY.conc = self.platos(num_plato).y_i;
                        HgiVm1T = 0;
                        HgiVp1T = 0;
                        for i = 1:self.num_sust
                            try 
                                cp = self.sust(i).cp_gi{1};
                            catch ME
                                error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                            end
                            if num_plato > 1                            
                                deltaHLm1T = integral(@(t) cp(t), 273.15, T - delta);
                                deltaHLp1T = integral(@(t) cp(t), 273.15, T + delta);
                                HgiLm1T = HgiLm1T + deltaHLm1T*mezclaX.conc(i) + Href(i) * mezclaX.conc(i);
                                HgiLp1T = HgiLp1T + deltaHLp1T*mezclaX.conc(i) + Href(i) * mezclaX.conc(i);
                            else 
                                HgiLm1T = 0;
                                HgiLp1T = 0;
                            end
                            deltaHVm1T = integral(@(t) cp(t), 273.15, T - delta);
                            deltaHVp1T = integral(@(t) cp(t), 273.15, T + delta);

                            HgiVm1T = HgiVm1T + deltaHVm1T*mezclaY.conc(i) + Href(i) * mezclaY.conc(i);
                            HgiVp1T = HgiVp1T + deltaHVp1T*mezclaY.conc(i) + Href(i) * mezclaY.conc(i);
                        end
                        if num_plato > 1
                            Hdep_bubXm1 = self.MEdE.entalpia(T - delta, P, mezclaX, 'liq');
                            Hdep_bubXp1 = self.MEdE.entalpia(T + delta, P, mezclaX, 'liq');
                        else 
                            Hdep_bubXm1 =  0;
                            Hdep_bubXp1 = 0;
                        end
                            Hdep_dewYm1 = self.MEdE.entalpia(T - delta, P, mezclaY, 'vap');
                            Hdep_dewYp1 = self.MEdE.entalpia(T + delta, P, mezclaY, 'vap');
                            dHVjdT = (HgiVp1T - HgiVm1T)/(2*delta) -(Hdep_dewYp1 - Hdep_dewYm1)/(2*delta);

                        dHLjdT = (HgiLp1T - HgiLm1T)/(2*delta) - (Hdep_bubXp1 - Hdep_bubXm1)/(2*delta);
                        HgiLm1m1T = 0;
                        HgiLm1p1T = 0;       % Hgasideal para calculo L etapa j-1 -delta
                        HgiVp1m1T = 0;
                        HgiVp1p1T = 0;
                        HgiLp1T = 0;
                        HgiLm1T = 0;
                        for i = 1:self.num_sust
                            try 
                                cp = self.sust(i).cp_gi{1};
                            catch ME
                                error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                            end
                            if num_plato > 1
                                deltaHVm1T = integral(@(t) cp(t), 273.15, TLm1 - delta);
                                deltaHVp1T = integral(@(t) cp(t), 273.15, TLm1 + delta);

                                HgiLm1T = HgiLm1T +  deltaHVm1T*mezclaXm1.conc(i)+ Href(i) * mezclaXm1.conc(i);
                                HgiLp1T = HgiLp1T + deltaHVp1T*mezclaXm1.conc(i)+ Href(i) * mezclaXm1.conc(i);
                            else
                                HgiLm1T = 0;
                                HgiLp1T = 0;
                            end
                        end

                        if num_plato > 1
                            Hdep_bubXm1 = self.MEdE.entalpia(TLm1 - delta, PLm1, mezclaXm1, 'liq');
                            Hdep_bubXp1 = self.MEdE.entalpia(TLm1 + delta, PLm1, mezclaXm1, 'liq');
                        else
                            Hdep_bubXp1 = 0;
                            Hdep_bubXm1 = 0;
                        end
                        dHLm1dT = (HgiLp1T - HgiLm1T)/(2*delta) - (Hdep_bubXp1 - Hdep_bubXm1)/(2*delta);   %Diferencial de HL etapa j-1 
                        HgiVp1m1T = 0;
                        HgiVp1p1T = 0;
                        if num_plato < self.etapas
                            for i = 1:self.num_sust
                                try 
                                    cp = self.sust(i).cp_gi{1};
                                catch ME
                                    error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                                end
                                deltaHVp1m1T = integral(@(t) cp(t), 273.15, TVp1 - delta);
                                deltaHVp1p1T = integral(@(t) cp(t), 273.15, TVp1 + delta);
                                HgiVp1m1T = HgiVp1m1T + deltaHVp1m1T*mezclaYp1.conc(i);
                                HgiVp1p1T = HgiVp1p1T + deltaHVp1p1T*mezclaYp1.conc(i);
                            end
                        else
                            HgiVp1m1T = 0;
                            HgiVp1p1T = 0;
                        end
                        if num_plato < self.etapas
                            Hdep_dewYp1 = self.MEdE.entalpia(TVp1 + delta, PVp1, mezclaYp1, 'vap');
                            Hdep_dewYm1 = self.MEdE.entalpia(TVp1 - delta, PVp1, mezclaYp1, 'vap');
                        else
                            Hdep_dewYp1 = 0;
                            Hdep_dewYm1 = 0;
                        end
                        dHVp1dT = (HgiVp1p1T - HgiVp1m1T)/(2*delta) - (Hdep_dewYp1 - Hdep_dewYm1)/(2*delta);   %Diferencial de HV etapa j+1

                        dQjdTjp1 = 0;
                        dQjdTjm1 = 0;
                        if num_plato > 1 && num_plato < self.etapas
                            self.Aij(1, (num_plato - 2)*((2*self.num_sust + 1))+2.*self.num_sust+1) =  -sum(self.platos(num_plato-1).l_i(:))*(dHLm1dT) - dQjdTjm1;
                        end
                        if num_plato < self.etapas && num_plato > 1
                            self.Cij(1, (num_plato - 1)*((2*self.num_sust + 1))+2.*self.num_sust+1) =  -sum(self.platos(num_plato+1).v_i(:))*(dHVp1dT) - dQjdTjp1;
                        end
                        for dFidli = 1:self.num_sust
                            %Funcion Fj siendo Hj
                            if num_plato < self.etapas && num_plato > 1
                                self.Bij(1, self.num_sust + (num_plato - 1)*(2*self.num_sust + 1)+dFidli) = HLj(num_plato)*(1+sj(num_plato));
                            end
                            if num_plato > 1 && num_plato < self.etapas
                                self.Aij(1, self.num_sust + (num_plato-2)*(2*self.num_sust + 1)+dFidli) = -HLj(num_plato - 1);
                            end
                        end
                        
                        if num_plato ~= 1  && num_plato ~= self.etapas
                            self.Bij(1, (num_plato - 1)*((2*self.num_sust + 1))+2*self.num_sust+1) = dHLjdT*(1+sj(num_plato))*sum(self.platos(num_plato).l_i(:)) + dHVjdT.*(1+Sj(num_plato))*sum(self.platos(num_plato).v_i(:));
                        end
                        for dFidvi = 1:self.num_sust
                            %Funcion Fj siendo Hj
                            if num_plato ~= 1 && num_plato ~= self.etapas
                                self.Bij(1, (num_plato-1)*(2*self.num_sust + 1) + dFidvi) = HVj(num_plato)*(1+Sj(num_plato));
                            end
                            if num_plato < self.etapas && num_plato > self.etapas
                                self.Cij(1, (num_plato-1)*(2*self.num_sust + 1)+dFidvi) = -HVj(num_plato + 1);
                            end

                        end
                        mezclaX = mezclaY;
                        mezclaX.conc = self.platos(num_plato).x_i;
                        %por diferenciacion numerica de 3 puntos centrada de Ej = 0 = Kij*lij*(sum(vj))/(sum(lj))
                        % siendo la funci�n de la temperatura solo dependiente de Kij
                        % Kijp1(Tj) = Kij(T0+dif) = Kij(T0 + 1E-7)
                        fGijp1 = self.MEdE.fugF(T+delta, P, mezclaY, 'vap');
                        fLijp1 = self.MEdE.fugF(T+delta, P, mezclaX, 'liq');
                        fGijm1 = self.MEdE.fugF(T-delta, P, mezclaY, 'vap');
                        fLijm1 = self.MEdE.fugF(T-delta, P, mezclaX, 'liq');
                    %                 fGij = MEdE.fugF(Tj, P, mezcla, 'vap');
                    %                 fLij = MEdE.fugF(Tj, P, mezcla, 'liq');
                        Kijp1 = fLijp1 ./ fGijp1;
                        Kijm1 = fLijm1 ./ fGijm1;
                        fGij = self.MEdE.fugF(T, P, mezclaY, 'vap');
                        fLij = self.MEdE.fugF(T, P, mezclaX, 'liq');
                        self.perfil_k(:,num_plato) = fLij ./ fGij;
                        for dFidvi = 1:self.num_sust
                            %Funcion Fj siendo Ej
                            dEjdTj = ((Kijp1(dFidvi) - Kijm1(dFidvi))./(2*delta))*(self.platos(num_plato).l_i(dFidvi)*(sum(self.platos(num_plato).v_i(:)))/(sum(self.platos(num_plato).l_i(:))));
                            self.Bij(self.num_sust + 1 + dFidvi, (num_plato-1)*(2*self.num_sust + 1) + 2*self.num_sust + 1) = dEjdTj;

                            self.Bij(dFidvi + 1, (num_plato-1)*(2*self.num_sust + 1) + dFidvi) = 1+Sj(num_plato);

                            self.Bij(self.num_sust + 1 + dFidvi, (num_plato-1)*(2*self.num_sust + 1) + dFidvi) = self.platos(num_plato).K(dFidvi)*self.platos(num_plato).l_i(dFidvi)*(1)/(sum(self.platos(num_plato).l_i(:))) - 1;
                        end
                        for dFidvi = self.num_sust:2
                            for dFidvj = 1:1:self.num_sust
                                self.Bij(self.num_sust + dFidvi, (num_plato-1)*(2*self.num_sust + 1) + dFidvi) =  self.platos(num_plato).K(dFidvj)*self.platos(num_plato).l_i(dFidvj)*(1)/(sum(self.platos(num_plato).l_i(:)));
                            end
                        end
                        for dFidvi = 1:self.num_sust
                            for dFidvj = 1:1:self.num_sust
                                if dFidvj ~=  dFidvi 
                                    self.Bij(self.num_sust + 1 + dFidvj, (num_plato-1)*(2*self.num_sust + 1) + dFidvi) = self.platos(num_plato).K(dFidvj)*self.platos(num_plato).l_i(dFidvj)*(1)/(sum(self.platos(num_plato).l_i(:)));
                                end
                            end
                        end
                        for dFidli = 1:self.num_sust
                            %Funcion Fj siendo Mij
                            self.Bij(dFidli + 1, (num_plato-1)*(2*self.num_sust + 1) + dFidli + self.num_sust) = 1+sj(num_plato);
                            %fprintf(1, '%f\n', self.platos(num_plato).K(num_plato, dFidli)*((sum(self.platos(num_plato).l_i(num_plato, :))-self.platos(num_plato).l_i(num_plato, dFidli))*(sum(self.platos(num_plato).v_i(num_plato, :))))/(sum(self.platos(num_plato).l_i(num_plato, :))^2))
                            self.Bij(self.num_sust + 1 + dFidli, self.num_sust + (num_plato - 1)*(2*self.num_sust + 1)+dFidli) = self.platos(num_plato).K(dFidli)*((sum(self.platos(num_plato).l_i(:))-self.platos(num_plato).l_i(dFidli))*(sum(self.platos(num_plato).v_i(:))))/(sum(self.platos(num_plato).l_i(:))^2);
                        end
                        for dFidli = 1:self.num_sust
                            for dFidlj = self.num_sust:-1:1
                                %Funcion Fj siendo Ej
                                if dFidlj ~= dFidli
                                    self.Bij(1 + dFidlj + self.num_sust, self.num_sust + (num_plato-1)*(2*self.num_sust + 1) + dFidli) = -self.platos(num_plato).K(dFidlj)*((sum(self.platos(num_plato).v_i(:)))*(self.platos(num_plato).l_i(dFidlj)))/(sum(self.platos(num_plato).l_i(:))^2);
                                end
                            end
                        end
                        for i = 1:self.num_sust
                            self.Fk(1+i+(num_plato - 1)*(2.*self.num_sust + 1)) = self.Fk(1+i+(num_plato - 1)*(2.*self.num_sust + 1)) - self.platos(num_plato).l_i(i)*(1+sj(num_plato)) - self.platos(num_plato).v_i(i)*(1+Sj(num_plato));
                            if num_plato > 1
                                self.Fk(1+i+(num_plato - 1)*(2.*self.num_sust + 1)) =  self.Fk(1+i+(num_plato - 1)*(2.*self.num_sust + 1)) + self.platos(num_plato-1).l_i(i);
                            end
                            if num_plato < self.etapas
                                self.Fk(1+i+(num_plato - 1)*(2.*self.num_sust + 1)) =  self.Fk(1+i+(num_plato - 1)*(2.*self.num_sust + 1)) + self.platos(num_plato+1).v_i(i);
                            end      
                        end

                        for i = 1:self.num_sust
                            self.Fk(1+self.num_sust+i +(num_plato - 1)*(2.*self.num_sust + 1)) = self.Fk(1+self.num_sust+i +(num_plato - 1)*(2.*self.num_sust + 1)) - self.platos(num_plato).K(i)*self.platos(num_plato).l_i(i)*((sum(self.platos(num_plato).v_i(:)))./(sum(self.platos(num_plato).l_i(:))))+self.platos(num_plato).v_i(i);     
                        end
                        indicee = 1;
                        if  any(num_plato == self.platos_entradas(:))    %Entalpia de alimentaci�n HF
                            if length(self.platos_entradas) == 1
                                indicee = 1;
                            else
                                indicee = find((self.platos_entradas(:) == num_plato) ~= 0, 1, 'first');
                            end
                            HF(num_plato) = alimen(indicee).H;
                        else
                            HF(num_plato) = 0;
                        end
                        if num_plato == 1
                            if isempty(self.salidas) || ~isa(self.salidas, 'cell')
                                self.Fk(1) = -(self.varX(self.lk) - self.dsubi(self.lk));
                            else
                                self.Fk(1) = -(self.varX(self.lk) - self.dsubi(self.lk));
                            end
                        elseif num_plato == self.etapas
                            self.Fk(1+self.num_sust+(num_plato - 1)*(2.*self.num_sust + 1)) = -(self.varX((num_plato-1)*(2*self.num_sust + 1) + self.num_sust + self.hk) - self.bsubi(self.hk));
                        else 
                            self.Fk(1+(num_plato - 1)*(2.*self.num_sust + 1)) =  - HLj(num_plato)*(1+sj(num_plato))*(sum(self.platos(num_plato).l_i(:))) - HVj(num_plato)*(1+Sj(num_plato))*(sum(self.platos(num_plato).v_i(:))) + HF(num_plato)*alimen(indicee).molF + self.perfil_q(num_plato-1) + HLj(num_plato - 1) * sum(self.platos(num_plato-1).l_i(:)) + HVj(num_plato + 1) * sum(self.platos(num_plato + 1).v_i(:));
                        end
                        if num_plato > 1 && num_plato < self.etapas
                            self.Cij(1,1+(num_plato-1)*(2*self.num_sust + 1):1+(num_plato-1)*(2*self.num_sust + 1)+self.num_sust -1) = -HVj(num_plato + 1);
                        end
                        T = TVp1;
                        P = PVp1;
                        mezclaY = mezclaYp1;
                    end
                    
                    [deltaXvar, self.matriz_tridiag ] = tridiagThomas(self.Bij, self.Aij, self.Cij, self.Fk);
                    condition = 0;
                    for num_plato = 1:self.etapas
                        nuevas_varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust) = self.varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust) + self.damping.*deltaXvar(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust);   
                        if any(nuevas_varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust-1) < 0)
                            nuevas_varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust - 1) = self.varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust-1).*exp((self.damping.*deltaXvar(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust-1))./(self.varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust-1)));
                        end
                        if num_plato == self.etapas && condition == 1
                            condition = 0;
                        end
                        Tbackup(num_plato) = self.perfil_t(num_plato);
                        Vbackup(num_plato) = self.perfil_v(num_plato);
                        Lbackup(num_plato) = self.perfil_l(num_plato);
                        self.perfil_t(num_plato) = nuevas_varX(1+(num_plato-1)*(2*(self.num_sust)+1)+2*self.num_sust);
                        if sum(nuevas_varX(1+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + self.num_sust-1)) > self.tol
                            self.perfil_v(num_plato) = sum(nuevas_varX(1+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + self.num_sust-1));
                        end
                        self.perfil_l(num_plato) = sum(nuevas_varX(1+self.num_sust+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + 2*self.num_sust-1));
                        self.varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust) = nuevas_varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust);

                        self.perfil_vi(num_plato, :) = self.varX(1+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + self.num_sust-1)';
                        self.perfil_li(num_plato, :) = self.varX(1+self.num_sust+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + 2*self.num_sust-1)';
                        if num_plato == 1
                            if ~isempty(self.salidas) && isa(self.salidas, 'cell') 
                                if self.salidas{1} == 1 && self.salidas{2} == 1 && self.perfil_v(num_plato)/(self.perfil_v(num_plato)+self.perfil_l(num_plato)) < 1e-4
                                    self.ysubj(num_plato,:) = self.varX(1+(num_plato)*(2*(self.num_sust)+1):1+(num_plato)*(2*(self.num_sust)+1) + self.num_sust-1)'./sum(self.varX(1+(num_plato)*(2*(self.num_sust)+1):1+(num_plato)*(2*(self.num_sust)+1) + self.num_sust-1));
                                else
                                    self.ysubj(num_plato,:) = self.varX(1:self.num_sust)'./sum(self.varX(1:self.num_sust));
                                end
                            else
                                self.ysubj(num_plato,:) = self.varX(1:self.num_sust)'./sum(self.varX(1:self.num_sust));
                            end
                        else
                            self.ysubj(num_plato,:) = self.varX(1+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + self.num_sust-1)'./sum(self.varX(1+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + self.num_sust-1));
                        end
                        self.xsubj(num_plato,:) = self.varX(1+self.num_sust+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + 2*self.num_sust-1)'./sum(self.varX(1+self.num_sust + (num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + 2*self.num_sust-1));                    
                        self.ysubj(num_plato,:) = self.ysubj(num_plato,:);
                        self.xsubj(num_plato,:) = self.xsubj(num_plato,:);
                        self.perfil_t(num_plato) = self.varX(1+(num_plato-1)*(2*(self.num_sust)+1)+2*self.num_sust);
                    end
                    for iiiiter = 1:self.etapas
                        if iiiiter == self.etapas
                            fprintf(1, ' %f \n', self.perfil_t(iiiiter));
                            break
                        end
                        fprintf(1, ' %f ', self.perfil_t(iiiiter));
                    end

                    for iteru = 1:self.etapas
                        self.platos(iteru).setV(self.perfil_v(iteru));
                        self.platos(iteru).setL(self.perfil_l(iteru));
                        self.platos(iteru).setyi(self.ysubj(iteru,:));                    
                        self.platos(iteru).setxi(self.xsubj(iteru,:));
                        self.platos(iteru).K = self.perfil_k(:, iteru);
                    end
                    bsubiclavepesado = self.bsubi(self.hk);
                    self.bsubi = self.xsubj(end,:).*(self.perfil_l(end));
                    self.bsubi(self.hk) = bsubiclavepesado;
                    if ~isempty(self.salidas) && isa(self.salidas, 'cell') 
                        if self.salidas{1} == 1 && self.salidas{2} == 1
                            dsubiclaveliviano = self.dsubi(self.lk);
                            self.dsubi = (self.dflujo - self.salidas{3}).* self.ysubj(1,:) + (self.salidas{3}).*self.xsubj(1,:);
                            self.dsubi(self.lk) = dsubiclaveliviano;
                        else
                            dsubiclaveliviano = self.dsubi(self.lk);
                            self.dsubi = (self.dflujo).* self.ysubj(1,:);     
                            self.dsubi(self.lk) = dsubiclaveliviano;
                        end
                    else
                        dsubiclaveliviano = self.dsubi(self.lk);
                        self.dsubi = (self.dflujo).* self.ysubj(1,:);
                        self.dsubi(self.lk) = dsubiclaveliviano;
                    end


                    for i = 1:self.etapas
                        self.platos(i).L = self.perfil_l(i);
                        self.platos(i).V = self.perfil_v(i);
                        self.platos(i).y_i = self.ysubj(i,:);
                        self.platos(i).x_i = self.xsubj(i,:);
                        self.platos(i).v_i = self.perfil_vi(i,:);
                        self.platos(i).l_i = self.perfil_li(i,:);
                    end
                    self.bflujo = sum(self.perfil_li(end,:));
                    if isempty(self.salidas) 
                        self.dflujo = sum(self.perfil_vi(1,:)); 
                    else
                        if isa(self.salidas, 'cell')
                            if self.salidas{1} == 1 && self.salidas{2} == 1
                                self.dflujo =  sum(self.perfil_vi(1,:)) + self.F - self.UmW - sum(self.perfil_li(end,:));
                            end
                        else
                            self.dflujo = sum(self.perfil_vi(1,:)); 
                        end
                    end
                    self.perfil_v = sum(self.perfil_vi, 2);
                    self.perfil_l = sum(self.perfil_li, 2);
                        if (sum(self.Fk.^2)) < respaldovalI
                            respaldovalI = (sum(self.Fk.^2));
                            respaldobi = self.bsubi;
                            respaldodi = self.dsubi;
                            respaldoV = self.perfil_v;
                            respaldoL = self.perfil_l;
                            respaldovi = self.perfil_vi;
                            respaldoli = self.perfil_li;
                            respaldot = self.perfil_t;
                            respaldok =  self.perfil_k;
                            respaldoplatos = self.platos;
                            respaldoxi = self.xsubj;
                            respaldoyi =  self.ysubj;
                            if isempty(self.salidas) || ~isa(self.salidas, 'cell')
                                self.qc = -self.perfil_v(1)*(HVj(1)) - HLj(1)*self.perfil_l(1) + self.perfil_v(2)*HVj(2);
                            else
                                if self.salidas{1} == 1 && self.salidas{2} == 1 
                                    self.qc = -self.perfil_v(1)*(HVj(1)) - HLj(1)*self.perfil_l(1) - sj(1)*self.perfil_l(1)*HLj(1) + self.perfil_v(2)*HVj(2);
                                elseif self.salidas{1} == 1 && self.salidas{2} == 0
                                    self.qc = -self.perfil_v(1)*(HVj(1)) - HLj(1)*self.perfil_l(1) - sj(1)*self.perfil_l(1)*HLj(1) + self.perfil_v(2)*HVj(2);
                                else
                                    self.qc = -self.perfil_v(1)*(HVj(1)) - HLj(1)*self.perfil_l(1) + self.perfil_v(2)*HVj(2);
                                end
                            end
                            Fh = 0;
                            Vh1 = 0;
                            Uh = 0; %Salidas Laterales Liquidas 1
                            Wh = 0; %Salidas Laterales Vapor 1
                            Lh = 0; %Liquido que va al plato 2
                            Vh = 0; %Pudiera haber vapor fuga del plato 1
                            Vh2 = 0; %Vapor que entra al plato 1 del plato 2
                            LhEND = 0;
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
                                if ~isempty(self.platos(iterx).salidaV)
                                    if iterx == 1
                                        Vh1 = Vh1 +  self.perfil_v(1) * HVj(iterx);  % que no se sume
                                            %2 veces el destilado vapor
                                    else    
                                        Wh = Wh + self.platos(iterx).salidaV * HVj(iterx);
                                    end
                                end
                                if ~isempty(self.platos(iterx).salidaL)
                                    if iterx == self.etapas
                                        LhEND = LhEND  + self.platos(iterx).salidaL * HLj(num_plato);
                                        %que no se sume dos veces el fondo l�quido
                                    else                            
                                        Uh = Uh + self.platos(iterx).salidaL * HLj(num_plato);
                                    end
                                end
                            end
                            self.qb = Fh - Uh  - Wh - Vh1 - LhEND - sum(self.perfil_q) - self.qc;
                            respaldoqc= self.qc;
                            respaldoqb = self.qb;
                            self.respaldoiter = iter;
                        end
                    
                        valI = sum(self.Fk.^2);
    %                     valI = valI + ((self.perfil_t(num_plato) - Tbackup(num_plato))/self.perfil_t(num_plato))^2;
    %                     if self.perfil_v(num_plato) > 0.00001 .* self.tol
    %                         valI = valI + ((self.perfil_v(num_plato) - Vbackup(num_plato))/self.perfil_v(num_plato))^2;
    %                     end
    %                     if self.perfil_l(num_plato) > 0.00001.* self.tol
    %                         valI = valI + ((self.perfil_l(num_plato) - Lbackup(num_plato))/self.perfil_l(num_plato))^2;
    %                     end
                    %fprintf(1, 'Tbackup = %f ', Tbackup(1));
                    %for i = 2:self.etapas - 1
                    %    fprintf(1, ' %f ', Tbackup(i));
                    %end
                    %fprintf(1, ' %f \n', Tbackup(end));
                    if isempty(self.qc)
                        self.qc = 0;
                    end
                    if isempty(self.qb)
                        self.qb = 0;
                    end
                    display(valI);
                    self.actualiter = iter;
                    self.actualvalI = valI;
                    self.respaldoV(self.actualiter + 1, :) = self.perfil_v;
                    self.respaldoL(self.actualiter + 1, :) = self.perfil_l;
                    self.respaldovi(self.actualiter.*self.etapas+1:self.actualiter.*self.etapas + self.etapas, :) = self.perfil_vi;
                    self.respaldoli(self.actualiter.*self.etapas+1:self.actualiter.*self.etapas + self.etapas, :) = self.perfil_li;
                    self.respaldot(self.actualiter+1, :) = self.perfil_t;
                    self.respaldoqc(self.actualiter+1) = self.qc;
                    self.respaldoqb(self.actualiter+1) = self.qb;
                    self.respaldovalI(self.actualiter + 1) = valI;
                

                    if valI < self.tol
                        self.laststep = logical(1);
                    end

                end
                
                
                if isempty(self.salidas) || ~isa(self.salidas, 'cell')
                    self.qc = -self.perfil_v(1)*(HVj(1)) - HLj(1)*self.perfil_l(1) + self.perfil_v(2)*HVj(2);
                else
                    if self.salidas{1} == 1 && self.salidas{2} == 1 
                        self.qc = -self.perfil_v(1)*(HVj(1)) - HLj(1)*self.perfil_l(1) - sj(1)*self.perfil_l(1)*HLj(1) + self.perfil_v(2)*HVj(2);
                    elseif self.salidas{1} == 1 && self.salidas{2} == 0
                        self.qc = -self.perfil_v(1)*(HVj(1)) - HLj(1)*self.perfil_l(1) - sj(1)*self.perfil_l(1)*HLj(1) + self.perfil_v(2)*HVj(2);
                    else
                        self.qc = -self.perfil_v(1)*(HVj(1)) - HLj(1)*self.perfil_l(1) + self.perfil_v(2)*HVj(2);
                    end
                end
                Fh = 0;
                Vh1 = 0;
                Uh = 0; %Salidas Laterales Liquidas 1
                Wh = 0; %Salidas Laterales Vapor 1
                Lh = 0; %Liquido que va al plato 2
                Vh = 0; %Pudiera haber vapor fuga del plato 1
                Vh2 = 0; %Vapor que entra al plato 1 del plato 2
                LhEND = 0;
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
                    if ~isempty(self.platos(iterx).salidaV)
                        if iterx == 1
                            Vh1 = Vh1 +  self.perfil_v(1) * HVj(iterx);  % que no se sume
                                %2 veces el destilado vapor
                        else    
                            Wh = Wh + self.platos(iterx).salidaV * HVj(iterx);
                        end
                    end
                    if ~isempty(self.platos(iterx).salidaL)
                        if iterx == self.etapas
                            LhEND = LhEND  + self.platos(iterx).salidaL * HLj(num_plato);
                            %que no se sume dos veces el fondo l�quido
                        else                            
                            Uh = Uh + self.platos(iterx).salidaL * HLj(num_plato);
                        end
                    end
                end
                self.qb = Fh - Uh  - Wh - Vh1 - LhEND - sum(self.perfil_q) - self.qc;
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
                self.bsubi = respaldobi;
                self.dsubi = respaldodi;
                self.dflujo = sum(respaldodi);
                self.bflujo = sum(respaldobi);
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
%             etapaentrada = zeros(1, floor(length(self.entradas)/2));
%             for i = 1:3:length(self.entradas)
%                 etapaentrada(floor(i/2+1)) = self.entradas{i};
%             end
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
                if ~isempty(etapasalida)
                    plato_salida = find(etapasalida == itei, 1, 'first');
                else 
                    plato_salida = {};
                end
                if itei == 1
                    if isempty(plato_salida)
                        if self.dflujo > 1e-3
                            plato_salida = 1;
                        end
                    end
                end                
                if isempty(plato_aliment) && isempty(plato_salida)
                    self.platos(1, itei) = Plato(itei, self.perfil_t(itei), self.perfil_p(itei), [], [], self.xsubj(itei,:), self.ysubj(itei,:), self.perfil_v(itei), [], self.perfil_k(:, itei)', [], [], [], [], []);
                elseif ~isempty(plato_aliment) && isempty(plato_salida)
                    self.platos(1, itei) = Plato(itei, self.perfil_t(itei), self.perfil_p(itei), [], [], self.xsubj(itei,:), self.ysubj(itei,:), self.perfil_v(itei), [], self.perfil_k(:, itei)',  self.alimentaciones{plato_aliment*2}, [], [], [], []);
                elseif ~isempty(plato_salida) && isempty(plato_aliment)
                    if isempty(self.salidas)
                        salidasy = self.dflujo;
                        self.platos(1, itei) = Plato(itei, self.perfil_t(itei), self.perfil_p(itei), [], [], self.xsubj(itei,:), self.ysubj(itei,:), self.perfil_v(itei), [], self.perfil_k(:, itei)', [], salidasy, [], [], []);
                    else
                        if abs(self.salidas{plato_salida*3-1}-1) <= 1e-4 %Si es liquido
                            salidasx = self.salidas{plato_salida*3};
                            salidasy = self.dflujo;
                            self.platos(1, itei) = Plato(itei, self.perfil_t(itei), self.perfil_p(itei), [], [], self.xsubj(itei,:), self.ysubj(itei,:), self.perfil_v(itei), [], self.perfil_k(:, itei)', [], salidasy, salidasx, [], []);
                        elseif abs(self.salidas(plato_salida*3-1)) <= 1e-4 %Si es vapor
                            salidasy = self.salidas{plato_salida*3};
                            self.platos(1, itei) = Plato(itei, self.perfil_t(itei), self.perfil_p(itei), [], [], self.xsubj(itei,:), self.ysubj(itei,:), self.perfil_v(itei),[],self.perfil_k(:, itei)',  [], salidasy,[], [], []);
                        end
                    end
                end
            end
            if self.dflujo > 1e-3
                self.perfil_v(1) = self.dflujo;
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
                        self.platos(1, itei) = Plato(itei, self.perfil_t(itei), self.perfil_p(itei), [], [], self.xsubj(:, itei), self.ysubj(:,itei), self.perfil_v(itei), self.perfil_l(itei), self.perfil_k(:, itei), [],salidasy , salidasx,[], []);
                    elseif abs(self.salidas(plato_salida*3-1)) <= 1e-4  %Si es vapor
                        salidasy = self.salidas{plato_salida*3};
                        self.platos(1, itei) = Plato(itei, self.perfil_t(itei), self.perfil_p(itei), [], [], self.xsubj(:, itei), self.ysubj(:,itei), self.perfil_v(itei), self.perfil_l(itei),self.perfil_k(:, itei),  [], salidasy,[], [], []);
                    end
                end
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
                valI = 1e308;
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
%                     for l = self.hk + 1:self.num_sust
%                         for m = 3:3:tamano(2)
%                             bsubi(l) = bsubi(l) + self.concfeed(m/3, l).*self.entradas{m};
%                         end
%                         bsubi(l) = bsubi(l) - salidasi(l);
%                     end
%                     if ~isempty(self.salidas) && isa(self.salidas, 'cell')
%                         if self.salidas{1} == 1 && self.salidas{2} == 1
%                             if abs(self.dflujo./self.salidas{3}) < 1e-3
%                                 recover = abs((self.salidas{3} - sum(dsubi))/(Flk));
%                             else
%                                 recover = abs((self.salidas{3}+ self.dflujo - sum(dsubi))/Flk);
%                             end
%                         elseif self.salidas{1} == 1 && self.salidas{2} == 0
%                             recover = abs((self.salidas{3}+ self.dflujo - sum(dsubi))/Flk);
%                         else
%                             recover = abs((self.dflujo - sum(dsubi))/(Flk));
%                         end
%                     else
%                         recover = abs((self.dflujo - sum(dsubi))/(Flk));
%                     end
%                     if ~isempty(self.salidas) && isa(self.salidas, 'cell')
%                         if self.salidas{1} == 1 && self.salidas{2} == 1
%                             if all((self.dflujo + self.salidas{3} - sum(dsubi) - recover.*7.50.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0 )
%                                 cociente = 7.50.*(recover*Flk) / ((recover*Flk) + Fhk);
%                             elseif all((self.dflujo + self.salidas{3} - sum(dsubi) - recover.*5.50.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0 )
%                                 cociente = 5.50.*(recover*Flk) / ((recover*Flk) + Fhk);
%                             elseif all((self.dflujo + self.salidas{3} - sum(dsubi) - recover.*4.00.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0 )
%                                 cociente = 4.00.*(recover*Flk) / ((recover*Flk) + Fhk);
%                             elseif all((self.dflujo + self.salidas{3} - sum(dsubi) - recover.*3.25.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0 )
%                                 cociente = 3.25.*(recover*Flk) / ((recover*Flk) + Fhk);
%                             elseif all((self.dflujo + self.salidas{3} - sum(dsubi) - recover.*2.75.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0 )
%                                 cociente = 2.75.*(recover*Flk) / ((recover*Flk) + Fhk);
%                             elseif all((self.dflujo + self.salidas{3} - sum(dsubi) - recover.*2.00.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0 )
%                                 cociente = 2.00.*(recover*Flk) / ((recover*Flk) + Fhk);
%                             elseif all((self.dflujo + self.salidas{3} - sum(dsubi) - recover.*1.35.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0 )
%                                 cociente = 1.35.*(recover*Flk) / ((recover*Flk) + Fhk);
%                             else
%                                 cociente = (recover*Flk) / ((recover*Flk) + Fhk);
%                             end
%                         elseif self.salidas{1} == 1 && self.salidas{2} == 0
%                             if all((self.dflujo + self.salidas{3} - sum(dsubi) - recover.*7.50.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0 )
%                                 cociente = 7.50.*(recover*Flk) / ((recover*Flk) + Fhk);
%                             elseif all((self.dflujo + self.salidas{3} - sum(dsubi) - recover.*5.50.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0 )
%                                 cociente = 5.50.*(recover*Flk) / ((recover*Flk) + Fhk);
%                             elseif all((self.dflujo + self.salidas{3} - sum(dsubi) - recover.*4.00.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0 )
%                                 cociente = 4.00.*(recover*Flk) / ((recover*Flk) + Fhk);
%                             elseif all((self.dflujo + self.salidas{3} - sum(dsubi) - recover.*3.25.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0 )
%                                 cociente = 3.25.*(recover*Flk) / ((recover*Flk) + Fhk);
%                             elseif all((self.dflujo + self.salidas{3} - sum(dsubi) - recover.*2.75.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0 )
%                                 cociente = 2.75.*(recover*Flk) / ((recover*Flk) + Fhk);
%                             elseif all((self.dflujo + self.salidas{3} - sum(dsubi) - recover.*2.00.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0 )
%                                 cociente = 2.00.*(recover*Flk) / ((recover*Flk) + Fhk);
%                             elseif all((self.dflujo + self.salidas{3} - sum(dsubi) - recover.*1.35.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0 )
%                                 cociente = 1.35.*(recover*Flk) / ((recover*Flk) + Fhk);
%                             else
%                                 cociente = (recover*Flk) / ((recover*Flk) + Fhk);
%                             end
%                         else
%                             if all((self.dflujo - sum(dsubi) - recover.*7.50.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0)
%                                 cociente = 7.50.*(recover*Flk) / ((recover*Flk) + Fhk);
%                             elseif all((self.dflujo - sum(dsubi) - recover.*5.50.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0)
%                                 cociente = 5.50.*(recover*Flk) / ((recover*Flk) + Fhk);
%                             elseif all((self.dflujo - sum(dsubi) - recover.*4.00.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0)
%                                 cociente = 4.00.*(recover*Flk) / ((recover*Flk) + Fhk);
%                             elseif all((self.dflujo - sum(dsubi) - recover.*3.25.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0)
%                                 cociente = 3.25.*(recover*Flk) / ((recover*Flk) + Fhk);
%                             elseif all((self.dflujo - sum(dsubi) - recover.*2.75.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0)
%                                 cociente = 2.75.*(recover*Flk) / ((recover*Flk) + Fhk);
%                             elseif all((self.dflujo - sum(dsubi) - recover.*2.00.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0)
%                                 cociente = 2.00.*(recover*Flk) / ((recover*Flk) + Fhk);
%                             elseif all((self.dflujo - sum(dsubi) - recover.*1.35.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0)
%                                 cociente = 1.35.*(recover*Flk) / ((recover*Flk) + Fhk);
%                             else
%                                 cociente = (recover*Flk) / ((recover*Flk) + Fhk);
%                             end
%                         end
%                     else
%                         if all((self.dflujo - sum(dsubi) - recover.*7.50.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0)
%                             cociente = 7.50.*(recover*Flk) / ((recover*Flk) + Fhk);
%                         elseif all((self.dflujo - sum(dsubi) - recover.*5.50.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0)
%                             cociente = 5.50.*(recover*Flk) / ((recover*Flk) + Fhk);
%                         elseif all((self.dflujo - sum(dsubi) - recover.*4.00.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0)
%                             cociente = 4.00.*(recover*Flk) / ((recover*Flk) + Fhk);
%                         elseif all((self.dflujo - sum(dsubi) - recover.*3.25.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0)
%                             cociente = 3.25.*(recover*Flk) / ((recover*Flk) + Fhk);
%                         elseif all((self.dflujo - sum(dsubi) - recover.*2.75.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0)
%                             cociente = 2.75.*(recover*Flk) / ((recover*Flk) + Fhk);
%                         elseif all((self.dflujo - sum(dsubi) - recover.*2.00.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0)
%                             cociente = 2.00.*(recover*Flk) / ((recover*Flk) + Fhk);
%                         elseif all((self.dflujo - sum(dsubi) - recover.*1.35.*(recover*Flk) / ((recover*Flk) + Fhk).*Flk) >0)
%                             cociente = 1.35.*(recover*Flk) / ((recover*Flk) + Fhk);
%                         else
%                             cociente = (recover*Flk) / ((recover*Flk) + Fhk);
%                         end
%                     end
%                     if ~isempty(self.salidas) && isa(self.salidas, 'cell')
%                         if self.salidas{1} == 1 && self.salidas{2} == 1
%                             if (self.dflujo - self.salidas{3}) > 1e-3
%                                 dsubi(self.hk) = (self.dflujo - self.salidas{3}) - sum(dsubi);
%                             else
%                                 dsubi(self.hk) = self.dflujo - sum(dsubi);
%                             end
%                         elseif self.salidas{1} == 1 && self.salidas{2} == 0
%                             dsubi(self.hk) = self.dflujo - sum(dsubi);
%                         else
%                             dsubi(self.hk) = self.dflujo - sum(dsubi);
%                         end
%                     else
                    if isempty(self.salidas) || ~isa(self.salidas, 'cell')
                        dsubi(self.lk) = self.dflujo - sum(dsubi);
                    elseif isempty(self.dflujo) && self.salidas{1} == 1
                        dsubi(self.lk) = self.salidas{3} - sum(dsubi);
                        self.dflujo = 0;
%                     end 
                    elseif ~isempty(self.dflujo) && self.salidas{1} == 1
                        if self.salidas{3} ~= self.dflujo
                            dsubi(self.hk-1) = self.salidas{3}  + self.dflujo - sum(dsubi);
                        else
                            self.dflujo = 0;
                            dsubi(self.hk-1) = self.salidas{3}  + self.dflujo - sum(dsubi);
                        end
                    end
%                    bsubi(self.hk) = Fhk - dsubi(self.hk);
                    bsubi = self.Fi - dsubi;
                    if any(bsubi < 0)
                        bsubi(bsubi < 0) = 0;
                    end
                    
                    %dsubi = [160, 365.39, 4.61, 1e-10, 1e-14];
                    %bsubi = [1e-11, 4.61, 235.39, 25, 5 ]; 
                    if isempty(self.perfil_t) || all(self.perfil_t == 0)
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
                        if self.dflujo-self.salidas{3} < 1e-3
                            self.perfil_v(1, 2:end) = self.perfil_v(1, 2:end) + ones(1, self.etapas-1).*(self.reflujo + 1).*self.salidas{3};
                        else
                            self.perfil_v(1, 2:end) = self.perfil_v(1, 2:end) + ones(1, self.etapas-1).*(self.reflujo + 1).*(self.salidas{3} + (self.dflujo - self.salidas{3}));
                            self.perfil_v(1,1) = self.dflujo - self.salidas{3};
                        end
                    elseif self.salidas{1} == 1 && self.salidas{2} == 0
                        self.perfil_v(1, 2:end) = self.perfil_v(1, 2:end) +  ones(1, self.etapas-1).*(self.reflujo+1).*((self.dflujo - self.salidas{3})+ self.salidas{3});
                        self.perfil_v(1,1) = (self.dflujo);
                    else
                        self.perfil_v(1, 2:end) = self.perfil_v(1, 2:end) +  ones(1, self.etapas-1).*(self.reflujo+1).*self.dflujo;
                        self.perfil_v(1,1) = self.dflujo;
                    end
                elseif all(self.perfil_v == zeros(1, self.etapas))
                    self.perfil_v(1, 2:end) = self.perfil_v(1, 2:end) +  ones(1, self.etapas-1).*(self.reflujo+1).*self.dflujo;
                    self.perfil_v(1,1) = self.dflujo;
                end
            elseif any(self.dsubi ~= 0)
                self.bsubi = self.Fi - self.dsubi;
                if any(self.bsubi < 0)
                    self.bsubi(self.bsubi < 0) = 0;
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
                self.dsubi = self.Fi - self.bsubi;
                if any(self.dsubi < 0)
                    self.dsubi(self.dsubi < 0) = 0;
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
            if isempty(self.bflujo)
                self.bflujo = self.F - self.dflujo;
            end
            self.perfil_l(end) = self.bflujo;
            
            for op = 2:2:length(self.alimentaciones)
                if self.alimentaciones{op}.q <= 1 && self.alimentaciones{op}.q  >= 0
                    self.perfil_v(self.alimentaciones{op-1}:end) = self.perfil_v(self.alimentaciones{op-1}:end) - self.alimentaciones{op}.molF*(1 - self.alimentaciones{op}.q);
                    if any(self.perfil_v < 0)
                        warning('La condicion de reflujo espeficificada es muy peque�a, pruebe una mayor')
                    end
                elseif self.alimentaciones{op}.q < 0
                    mezclaX = Mezcla(self.sust, self.xsubj(self.alimentaciones{op-1},:), self.alimentaciones{2}.mezcla.kij);
                    mezclaY = Mezcla(self.sust, self.ysubj(self.alimentaciones{op-1},:), self.alimentaciones{2}.mezcla.kij);
                    Hvap = self.MEdE.entalpia(self.perfil_t(self.alimentaciones{op-1}), self.perfil_p(self.alimentaciones{op-1}), mezclaX, 'liq') - self.MEdE.entalpia(self.perfil_t(self.alimentaciones{op-1}), self.perfil_p(self.alimentaciones{op-1}), mezclaY, 'liq');
                    self.perfil_l(self.alimentaciones{op-1}:end-1) =  self.perfil_l(self.alimentaciones{op-1}:end-1) + ( self.alimentaciones{op}.q*self.alimentaciones{op}.molF*Hvap) / Hvap;
                    self.perfil_v(self.alimentaciones{op-1}:end) = self.perfil_v(self.alimentaciones{op-1}:end) - self.alimentaciones{op}.molF + ( self.alimentaciones{op}.q*self.alimentaciones{op}.molF*Hvap) / Hvap;
                    if any(self.perfil_v < 0)
                        error('La condicion de reflujo espeficificada es muy peque�a, pruebe una mayor')
                    end
                else
                    mezclaX = Mezcla(self.sust, self.xsubj(self.alimentaciones{op-1},:), self.alimentaciones{2}.mezcla.kij);
                    mezclaY = Mezcla(self.sust, self.ysubj(self.alimentaciones{op-1},:), self.alimentaciones{2}.mezcla.kij);
                    Hvap = self.MEdE.entalpia(self.perfil_t(self.alimentaciones{op-1}), self.perfil_p(self.alimentaciones{op-1}), mezclaX, 'liq') - self.MEdE.entalpia(self.perfil_t(self.alimentaciones{op-1}), self.perfil_p(self.alimentaciones{op-1}), mezclaY, 'liq');
                    self.perfil_l(self.alimentaciones{op-1}:end-1) =  self.perfil_l(self.alimentaciones{op-1}:end-1) + ( (self.alimentaciones{op}.q - 1)*self.alimentaciones{op}.molF*Hvap) / Hvap;
                    self.perfil_v(self.alimentaciones{op-1}:end) = self.perfil_v(self.alimentaciones{op-1}:end) + ( (self.alimentaciones{op}.q - 1)*self.alimentaciones{op}.molF*Hvap) / Hvap;
                end
            end
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
        function avanza1paso2(self)
            if isempty(self.perfil_k) || isempty(self.perfil_vi) || isempty(self.perfil_t) || isempty(self.perfil_li)
                self = self.balanmasa();
            else
                self.generar_etapas();
            end
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
            self.perfil_v(1) = self.dflujo;
            if ~isempty(self.salidas) && isa(self.salidas, 'cell') 
                if self.salidas{1} == 1 && self.salidas{2} == 1
                    self.dflujo = self.dflujo + self.salidas{3};
                elseif self.salidas{1} == 1 && self.salidas{2} == 0
                    self.dflujo = self.dflujo + self.salidas{3};
                end
            end
            if isempty(self.tol)
                self.tol = self.convergence;
            end
            if isempty(self.damping) || ~isa(self.damping, 'double') 
                self.damping = 1;
            elseif isa(self.damping, 'double') && (self.damping < 0  || self.damping > 2) 
                self.damping = 1;
            end
            nuevas_varX = zeros((2.*self.num_sust + 1)*self.etapas, 1);
            self.varX = nuevas_varX;
            Vbackup = zeros(1, self.etapas);
            Tbackup = zeros(1, self.etapas);
            Lbackup = zeros(1, self.etapas);
            for iitter0 = 1 : self.etapas
                self.varX(1+(iitter0-1)*(2*(self.num_sust)+1):1+(iitter0-1)*(2*(self.num_sust)+1) + self.num_sust - 1) = self.perfil_vi(iitter0,:);
                self.varX(1+self.num_sust+(iitter0-1)*(2*(self.num_sust)+1):1+(iitter0-1)*(2*(self.num_sust)+1) + 2*self.num_sust-1) = self.perfil_li(iitter0,:);
                self.varX(1+(iitter0-1)*(2*(self.num_sust)+1)+2*self.num_sust) = self.perfil_t(iitter0);
            end
            if isempty(self.respaldovalI)
                valI = 1e308;
            else
                valI = self.respaldovalI;
            end
            iter = 0;
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
                if isempty(self.respaldovalI)
                    self.respaldovalI = zeros(1, self.etapas);
                end
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
                if isempty(self.respaldovalI)
                    self.respaldovalI = zeros(1, self.etapas);
                end
                
                    valI = 0;
                    iter = iter + 1;
                    respaldandoiter = iter;
                    if iter == 1
                        for i = 1:self.etapas
                            self.platos(i).L = self.perfil_l(i);
                            self.platos(i).V = self.perfil_v(i);
                            self.platos(i).y_i = self.ysubj(i,:);
                            self.platos(i).x_i = self.xsubj(i,:);
                            self.platos(i).v_i = self.perfil_vi(i,:);
                            self.platos(i).l_i = self.perfil_li(i,:);
                        end
                    end
                    self.Bij = zeros(2*self.num_sust +1, (self.etapas)*(2*self.num_sust + 1));
                    self.Aij = zeros(2*self.num_sust +1, (self.etapas-1)*(2*self.num_sust + 1));
                    self.Cij = zeros(2*self.num_sust +1, (self.etapas-1)*(2*self.num_sust + 1));
                    self.Fk = zeros((2*self.num_sust+1).*self.etapas, 1);
                    tamano = size(self.entradas);
                    HVj = zeros(1, self.etapas);
                    HLj = zeros(1, self.etapas);
                    sj = zeros(1, self.etapas);
                    Sj = zeros(1, self.etapas);
                    for etapa = 1:self.etapas 
                        if etapa < self.etapas
                            for dMidli = 1:self.num_sust
                                self.Cij(dMidli+1, (etapa-1)*(2*self.num_sust + 1)+dMidli) = -1;
                            end
                        end
                        if etapa > 1
                            for dMidli = 1:self.num_sust
                                self.Aij(dMidli+1, self.num_sust + (etapa-2)*(2*self.num_sust + 1)+dMidli) = -1;
                            end
                        end
                    end
                    self.platos_entradas = zeros(1, length(self.entradas)/3);
                    alimen = Corriente.empty(0, length(self.entradas)/3 );
                    for itere=1:3:length(self.entradas)
                        self.platos_entradas((itere-1)/3+1) = self.entradas{itere};
                    end
                    for itere = 1:2:length(self.alimentaciones)
                        alimen((itere-1)/2+1) = self.alimentaciones{itere + 1};
                    end

                    for num_plato = 1:self.etapas
                        if ~isempty(self.platos(num_plato).salidaL);
                                if num_plato ~= self.etapas
                                    if self.platos(num_plato).salidaL > 1e-4
                                        sj(num_plato) = self.platos(num_plato).salidaL/self.platos(num_plato).L;
                                    else
                                        sj(num_plato) = 0;
                                    end
                                else
                                    sj(num_plato) = 0;
                                end
                        else
                            sj(num_plato) = 0;
                        end
                        if num_plato == 1
                            if sj(num_plato) ~= 0 
                                self.Bij(1,  1:self.num_sust) = -self.reflujo-1;
                                self.Bij(1 , 1+self.num_sust + (num_plato-1)*(2*self.num_sust +1):  1+self.num_sust + (num_plato-1)*(2*self.num_sust +1) + self.num_sust - 1) = sj(num_plato)*(-self.reflujo -1);
                                self.Cij(1, 1:self.num_sust) = 1;
                            else
                                self.Bij(1,  1:self.num_sust) = -self.reflujo;
                                self.Bij(1 , 1+self.num_sust + (num_plato-1)*(2*self.num_sust +1):  1+self.num_sust + (num_plato-1)*(2*self.num_sust +1) + self.num_sust - 1) = 1;
                            end
                        end
                        % La derivada de la funci�n de reemplazo OW es siempre +1
                        % para todos los Li y -reflujo (osea -R) para todos los vi
                        if num_plato == self.etapas
                            self.Bij(1, 1+self.num_sust+(num_plato-1)*(2*self.num_sust + 1):1+self.num_sust+(num_plato-1)*(2*self.num_sust + 1) + self.num_sust-1) = 1;
                        end
                        for l = 1:self.num_sust
                            for m = 3:3:tamano(2)
                                if num_plato == self.entradas{m-2}
                                    self.Fk(1+l+(num_plato - 1)*(2.*self.num_sust + 1)) = self.Fk(1+l+(num_plato - 1)*(2.*self.num_sust + 1)) + self.entradas{m}.*self.concfeed(m/3, l);
                                end
                            end
                        end
                        if num_plato == 1
                            T = self.varX(1+(num_plato-1)*(2*(self.num_sust)+1)+2*self.num_sust);   % Temperatura Alternativamente self.platos(num_plato).T
                            P = self.platos(num_plato).P;
                            mezclaY = Mezcla(self.sust, self.platos(num_plato).y_i, self.alimentaciones{2}.mezcla.kij);
                            HgiV = 0;
                            Href = zeros(1, self.num_sust);
                            for i = 1:self.num_sust
                                Href(i) = self.sust(i).href;
                                try 
                                    cp = self.sust(i).cp_gi{1};
                                catch ME
                                    error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                                end
                                deltaHL = integral(@(t) cp(t), 273.15, T);
                                deltaHV = deltaHL;
                                HgiV = HgiV + deltaHV*mezclaY.conc(i) + Href(i)*mezclaY.conc(i);
                            end
                            Hdep_dewY = self.MEdE.entalpia(T, P, mezclaY, 'vap');
                            Hdep_refY = self.MEdE.entalpia(273.15, 101.325, mezclaY, 'liq');
                            HVj(num_plato) = HgiV - Hdep_dewY + Hdep_refY;  % HVj
                        end
                        mezclaX = Mezcla(self.sust, self.platos(num_plato).x_i, self.alimentaciones{2}.mezcla.kij);
                            HgiL = 0;
                            Href = zeros(1, self.num_sust);

                            for i = 1:self.num_sust
                                Href(i) = self.sust(i).href;
                                try 
                                    cp = self.sust(i).cp_gi{1};
                                catch ME
                                    error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                                end
                                deltaHL = integral(@(t) cp(t), 273.15, T);
                                HgiL =  HgiL + deltaHL*mezclaX.conc(i) + Href(i)*mezclaX.conc(i);
                            end
                            Hdep_bubX = self.MEdE.entalpia(T, P, mezclaX, 'liq');
                            Hdep_refX = self.MEdE.entalpia(273.15, 101.325, mezclaX, 'liq');
                            HLj(num_plato) = HgiL - Hdep_bubX + Hdep_refX; % HLj 
                        if num_plato < self.etapas
                            TVp1 = self.varX(1+(num_plato)*(2*(self.num_sust)+1)+2*self.num_sust);   % Temperatura Alternativamente self.platos(num_plato).T
                            PVp1 = self.platos(num_plato+1).P;
                            Ycj = self.ysubj(num_plato + 1,:);
                            mezclaYp1 = Mezcla(self.sust, Ycj, self.alimentaciones{2}.mezcla.kij);
                            if num_plato > 1
                                TLm1 =  self.varX(1+(num_plato-2)*(2*(self.num_sust)+1)+2*self.num_sust);
                                PLm1 = self.platos(num_plato-1).P;
                                Xcjm1 = self.xsubj(num_plato - 1,:);
                                mezclaXm1 = Mezcla(self.sust, Xcjm1, self.alimentaciones{2}.mezcla.kij);
                            end
                            HgiVp1 =0;
                            HgiL = 0;
                            HgiV = 0;
                            Href = zeros(1, self.num_sust);
                            for i = 1:self.num_sust
                                Href(i) = self.sust(i).href;
                                try 
                                    cp = self.sust(i).cp_gi{1};
                                catch ME
                                    error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                                end
                                deltaHV = integral(@(t) cp(t), 273.15, TVp1);
                                HgiVp1 = HgiVp1 + deltaHV*mezclaYp1.conc(i) + Href(i)*mezclaYp1.conc(i);
                            end
                            Hdep_dewY = self.MEdE.entalpia(TVp1, PVp1, mezclaYp1, 'vap');
                            Hdep_refY = self.MEdE.entalpia(273.15, 101.325, mezclaYp1, 'liq');
                            HVj(num_plato+1) = HgiVp1 - Hdep_dewY + Hdep_refY;  % HVj
                        end
                        if ~isempty(self.platos(num_plato).salidaV)
                            if num_plato ~=1
                                if self.platos(num_plato).salidaV > 1e-4
                                    Sj(num_plato) = self.platos(num_plato).salidaV/self.platos(num_plato).V;
                                else
                                    Sj(num_plato) = 0;
                                end
                            else
                                Sj(num_plato) = 0;
                            end
                        else
                            Sj(num_plato) = 0;
                        end
                        HgiLm1T = 0;
                        HgiLp1T = 0;

                        HgiVm1T = 0;
                        HgiVp1T = 0;
                        delta = 1e-8;


                        mezclaY.conc = self.platos(num_plato).y_i;
                        HgiVm1T = 0;
                        HgiVp1T = 0;
                        for i = 1:self.num_sust
                            try 
                                cp = self.sust(i).cp_gi{1};
                            catch ME
                                error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                            end
                            if num_plato > 1                            
                                deltaHLm1T = integral(@(t) cp(t), 273.15, T - delta);
                                deltaHLp1T = integral(@(t) cp(t), 273.15, T + delta);
                                HgiLm1T = HgiLm1T + deltaHLm1T*mezclaX.conc(i) + Href(i) * mezclaX.conc(i);
                                HgiLp1T = HgiLp1T + deltaHLp1T*mezclaX.conc(i) + Href(i) * mezclaX.conc(i);
                            else 
                                HgiLm1T = 0;
                                HgiLp1T = 0;
                            end
                            deltaHVm1T = integral(@(t) cp(t), 273.15, T - delta);
                            deltaHVp1T = integral(@(t) cp(t), 273.15, T + delta);

                            HgiVm1T = HgiVm1T + deltaHVm1T*mezclaY.conc(i) + Href(i) * mezclaY.conc(i);
                            HgiVp1T = HgiVp1T + deltaHVp1T*mezclaY.conc(i) + Href(i) * mezclaY.conc(i);
                        end
                        if num_plato > 1
                            Hdep_bubXm1 = self.MEdE.entalpia(T - delta, P, mezclaX, 'liq');
                            Hdep_bubXp1 = self.MEdE.entalpia(T + delta, P, mezclaX, 'liq');
                        else 
                            Hdep_bubXm1 =  0;
                            Hdep_bubXp1 = 0;
                        end
                            Hdep_dewYm1 = self.MEdE.entalpia(T - delta, P, mezclaY, 'vap');
                            Hdep_dewYp1 = self.MEdE.entalpia(T + delta, P, mezclaY, 'vap');
                            dHVjdT = (HgiVp1T - HgiVm1T)/(2*delta) -(Hdep_dewYp1 - Hdep_dewYm1)/(2*delta);

                        dHLjdT = (HgiLp1T - HgiLm1T)/(2*delta) - (Hdep_bubXp1 - Hdep_bubXm1)/(2*delta);
                        HgiLm1m1T = 0;
                        HgiLm1p1T = 0;       % Hgasideal para calculo L etapa j-1 -delta
                        HgiVp1m1T = 0;
                        HgiVp1p1T = 0;
                        HgiLp1T = 0;
                        HgiLm1T = 0;
                        for i = 1:self.num_sust
                            try 
                                cp = self.sust(i).cp_gi{1};
                            catch ME
                                error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                            end
                            if num_plato > 1
                                deltaHVm1T = integral(@(t) cp(t), 273.15, TLm1 - delta);
                                deltaHVp1T = integral(@(t) cp(t), 273.15, TLm1 + delta);

                                HgiLm1T = HgiLm1T +  deltaHVm1T*mezclaXm1.conc(i)+ Href(i) * mezclaXm1.conc(i);
                                HgiLp1T = HgiLp1T + deltaHVp1T*mezclaXm1.conc(i)+ Href(i) * mezclaXm1.conc(i);
                            else
                                HgiLm1T = 0;
                                HgiLp1T = 0;
                            end
                        end

                        if num_plato > 1
                            Hdep_bubXm1 = self.MEdE.entalpia(TLm1 - delta, PLm1, mezclaXm1, 'liq');
                            Hdep_bubXp1 = self.MEdE.entalpia(TLm1 + delta, PLm1, mezclaXm1, 'liq');
                        else
                            Hdep_bubXp1 = 0;
                            Hdep_bubXm1 = 0;
                        end
                        dHLm1dT = (HgiLp1T - HgiLm1T)/(2*delta) - (Hdep_bubXp1 - Hdep_bubXm1)/(2*delta);   %Diferencial de HL etapa j-1 
                        HgiVp1m1T = 0;
                        HgiVp1p1T = 0;
                        if num_plato < self.etapas
                            for i = 1:self.num_sust
                                try 
                                    cp = self.sust(i).cp_gi{1};
                                catch ME
                                    error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                                end
                                deltaHVp1m1T = integral(@(t) cp(t), 273.15, TVp1 - delta);
                                deltaHVp1p1T = integral(@(t) cp(t), 273.15, TVp1 + delta);
                                HgiVp1m1T = HgiVp1m1T + deltaHVp1m1T*mezclaYp1.conc(i);
                                HgiVp1p1T = HgiVp1p1T + deltaHVp1p1T*mezclaYp1.conc(i);
                            end
                        else
                            HgiVp1m1T = 0;
                            HgiVp1p1T = 0;
                        end
                        if num_plato < self.etapas
                            Hdep_dewYp1 = self.MEdE.entalpia(TVp1 + delta, PVp1, mezclaYp1, 'vap');
                            Hdep_dewYm1 = self.MEdE.entalpia(TVp1 - delta, PVp1, mezclaYp1, 'vap');
                        else
                            Hdep_dewYp1 = 0;
                            Hdep_dewYm1 = 0;
                        end
                        dHVp1dT = (HgiVp1p1T - HgiVp1m1T)/(2*delta) - (Hdep_dewYp1 - Hdep_dewYm1)/(2*delta);   %Diferencial de HV etapa j+1

                        dQjdTjp1 = 0;
                        dQjdTjm1 = 0;
                        if num_plato > 1 && num_plato < self.etapas
                            self.Aij(1, (num_plato - 2)*((2*self.num_sust + 1))+2.*self.num_sust+1) =  -sum(self.platos(num_plato-1).l_i(:))*(dHLm1dT) - dQjdTjm1;
                        end
                        if num_plato < self.etapas && num_plato > 1
                            self.Cij(1, (num_plato - 1)*((2*self.num_sust + 1))+2.*self.num_sust+1) =  -sum(self.platos(num_plato+1).v_i(:))*(dHVp1dT) - dQjdTjp1;
                        end
                        for dFidli = 1:self.num_sust
                            %Funcion Fj siendo Hj
                            if num_plato < self.etapas && num_plato > 1
                                self.Bij(1, self.num_sust + (num_plato - 1)*(2*self.num_sust + 1)+dFidli) = HLj(num_plato)*(1+sj(num_plato));
                            end
                            if num_plato > 1 && num_plato < self.etapas
                                self.Aij(1, self.num_sust + (num_plato-2)*(2*self.num_sust + 1)+dFidli) = -HLj(num_plato - 1);
                            end
                        end

                        if num_plato ~= 1  && num_plato ~= self.etapas
                            self.Bij(1, (num_plato - 1)*((2*self.num_sust + 1))+2*self.num_sust+1) = dHLjdT*(1+sj(num_plato))*sum(self.platos(num_plato).l_i(:)) + dHVjdT.*(1+Sj(num_plato))*sum(self.platos(num_plato).v_i(:));
                        end
                        for dFidvi = 1:self.num_sust
                            %Funcion Fj siendo Hj
                            if num_plato ~= 1 && num_plato ~= self.etapas
                                self.Bij(1, (num_plato-1)*(2*self.num_sust + 1) + dFidvi) = HVj(num_plato)*(1+Sj(num_plato));
                            end
                            if num_plato < self.etapas && num_plato > self.etapas
                                self.Cij(1, (num_plato-1)*(2*self.num_sust + 1)+dFidvi) = -HVj(num_plato + 1);
                            end

                        end
                        mezclaX = mezclaY;
                        mezclaX.conc = self.platos(num_plato).x_i;
                        %por diferenciacion numerica de 3 puntos centrada de Ej = 0 = Kij*lij*(sum(vj))/(sum(lj))
                        % siendo la funci�n de la temperatura solo dependiente de Kij
                        % Kijp1(Tj) = Kij(T0+dif) = Kij(T0 + 1E-7)
                        fGijp1 = self.MEdE.fugF(T+delta, P, mezclaY, 'vap');
                        fLijp1 = self.MEdE.fugF(T+delta, P, mezclaX, 'liq');
                        fGijm1 = self.MEdE.fugF(T-delta, P, mezclaY, 'vap');
                        fLijm1 = self.MEdE.fugF(T-delta, P, mezclaX, 'liq');
                    %                 fGij = MEdE.fugF(Tj, P, mezcla, 'vap');
                    %                 fLij = MEdE.fugF(Tj, P, mezcla, 'liq');
                        Kijp1 = fLijp1 ./ fGijp1;
                        Kijm1 = fLijm1 ./ fGijm1;
                        fGij = self.MEdE.fugF(T, P, mezclaY, 'vap');
                        fLij = self.MEdE.fugF(T, P, mezclaX, 'liq');
                        self.perfil_k(:,num_plato) = fLij ./ fGij;
                        for dFidvi = 1:self.num_sust
                            %Funcion Fj siendo Ej
                            dEjdTj = ((Kijp1(dFidvi) - Kijm1(dFidvi))./(2*delta))*(self.platos(num_plato).l_i(dFidvi)*(sum(self.platos(num_plato).v_i(:)))/(sum(self.platos(num_plato).l_i(:))));
                            self.Bij(self.num_sust + 1 + dFidvi, (num_plato-1)*(2*self.num_sust + 1) + 2*self.num_sust + 1) = dEjdTj;

                            self.Bij(dFidvi + 1, (num_plato-1)*(2*self.num_sust + 1) + dFidvi) = 1+Sj(num_plato);

                            self.Bij(self.num_sust + 1 + dFidvi, (num_plato-1)*(2*self.num_sust + 1) + dFidvi) = self.platos(num_plato).K(dFidvi)*self.platos(num_plato).l_i(dFidvi)*(1)/(sum(self.platos(num_plato).l_i(:))) - 1;
                        end
                        for dFidvi = self.num_sust:2
                            for dFidvj = 1:1:self.num_sust
                                self.Bij(self.num_sust + dFidvi, (num_plato-1)*(2*self.num_sust + 1) + dFidvi) =  self.platos(num_plato).K(dFidvj)*self.platos(num_plato).l_i(dFidvj)*(1)/(sum(self.platos(num_plato).l_i(:)));
                            end
                        end
                        for dFidvi = 1:self.num_sust
                            for dFidvj = 1:1:self.num_sust
                                if dFidvj ~=  dFidvi 
                                    self.Bij(self.num_sust + 1 + dFidvj, (num_plato-1)*(2*self.num_sust + 1) + dFidvi) = self.platos(num_plato).K(dFidvj)*self.platos(num_plato).l_i(dFidvj)*(1)/(sum(self.platos(num_plato).l_i(:)));
                                end
                            end
                        end
                        for dFidli = 1:self.num_sust
                            %Funcion Fj siendo Mij
                            self.Bij(dFidli + 1, (num_plato-1)*(2*self.num_sust + 1) + dFidli + self.num_sust) = 1+sj(num_plato);
                            %fprintf(1, '%f\n', self.platos(num_plato).K(num_plato, dFidli)*((sum(self.platos(num_plato).l_i(num_plato, :))-self.platos(num_plato).l_i(num_plato, dFidli))*(sum(self.platos(num_plato).v_i(num_plato, :))))/(sum(self.platos(num_plato).l_i(num_plato, :))^2))
                            self.Bij(self.num_sust + 1 + dFidli, self.num_sust + (num_plato - 1)*(2*self.num_sust + 1)+dFidli) = self.platos(num_plato).K(dFidli)*((sum(self.platos(num_plato).l_i(:))-self.platos(num_plato).l_i(dFidli))*(sum(self.platos(num_plato).v_i(:))))/(sum(self.platos(num_plato).l_i(:))^2);
                        end
                        for dFidli = 1:self.num_sust
                            for dFidlj = self.num_sust:-1:1
                                %Funcion Fj siendo Ej
                                if dFidlj ~= dFidli
                                    self.Bij(1 + dFidlj + self.num_sust, self.num_sust + (num_plato-1)*(2*self.num_sust + 1) + dFidli) = -self.platos(num_plato).K(dFidlj)*((sum(self.platos(num_plato).v_i(:)))*(self.platos(num_plato).l_i(dFidlj)))/(sum(self.platos(num_plato).l_i(:))^2);
                                end
                            end
                        end
                        for i = 1:self.num_sust
                            self.Fk(1+i+(num_plato - 1)*(2.*self.num_sust + 1)) = self.Fk(1+i+(num_plato - 1)*(2.*self.num_sust + 1)) - self.platos(num_plato).l_i(i)*(1+sj(num_plato)) - self.platos(num_plato).v_i(i)*(1+Sj(num_plato));
                            if num_plato > 1
                                self.Fk(1+i+(num_plato - 1)*(2.*self.num_sust + 1)) =  self.Fk(1+i+(num_plato - 1)*(2.*self.num_sust + 1)) + self.platos(num_plato-1).l_i(i);
                            end
                            if num_plato < self.etapas
                                self.Fk(1+i+(num_plato - 1)*(2.*self.num_sust + 1)) =  self.Fk(1+i+(num_plato - 1)*(2.*self.num_sust + 1)) + self.platos(num_plato+1).v_i(i);
                            end      
                        end

                        for i = 1:self.num_sust
                            self.Fk(1+self.num_sust+i +(num_plato - 1)*(2.*self.num_sust + 1)) = self.Fk(1+self.num_sust+i +(num_plato - 1)*(2.*self.num_sust + 1)) - self.platos(num_plato).K(i)*self.platos(num_plato).l_i(i)*((sum(self.platos(num_plato).v_i(:)))./(sum(self.platos(num_plato).l_i(:))))+self.platos(num_plato).v_i(i);     
                        end
                        indicee = 1;
                        if  any(num_plato == self.platos_entradas(:))    %Entalpia de alimentaci�n HF
                            if length(self.platos_entradas) == 1
                                indicee = 1;
                            else
                                indicee = find((self.platos_entradas(:) == num_plato) ~= 0, 1, 'first');
                            end
                            HF(num_plato) = alimen(indicee).H;
                        else
                            HF(num_plato) = 0;
                        end
                        if num_plato == 1
                            if sj(num_plato) == 0
                                self.Fk(1) = -self.platos(1).L  + self.reflujo * self.platos(1).V;
                            else
                                self.Fk(1) = -sum(self.varX(2*self.num_sust + 1+1:2*self.num_sust + 1 + self.num_sust))  + self.reflujo * (sum(self.varX(1:self.num_sust)) + sj(num_plato)*(sum(self.varX(self.num_sust+1:2*self.num_sust)))) + sj(num_plato)*(sum(self.varX(self.num_sust+1:2*self.num_sust))) + sum(self.varX(1:self.num_sust));
                            end
                        elseif num_plato == self.etapas
                            self.Fk(1+(num_plato - 1)*(2.*self.num_sust + 1)) = -(sum(self.varX(1 + self.num_sust + (num_plato - 1)*(2*self.num_sust + 1):1+self.num_sust +(num_plato - 1)*(2*self.num_sust + 1)+self.num_sust-1)) - self.bflujo);
                        else 
                            self.Fk(1+(num_plato - 1)*(2.*self.num_sust + 1)) =  - HLj(num_plato)*(1+sj(num_plato))*(sum(self.platos(num_plato).l_i(:))) - HVj(num_plato)*(1+Sj(num_plato))*(sum(self.platos(num_plato).v_i(:))) + HF(num_plato)*alimen(indicee).molF + self.perfil_q(num_plato-1) + HLj(num_plato - 1) * sum(self.platos(num_plato-1).l_i(:)) + HVj(num_plato + 1) * sum(self.platos(num_plato + 1).v_i(:));
                        end
                        if num_plato > 1 && num_plato < self.etapas
                            self.Cij(1,1+(num_plato-1)*(2*self.num_sust + 1):1+(num_plato-1)*(2*self.num_sust + 1)+self.num_sust -1) = -HVj(num_plato + 1);
                        end
                        T = TVp1;
                        P = PVp1;
                        mezclaY = mezclaYp1;
                    end

                    [deltaXvar, self.matriz_tridiag ] = tridiagThomas(self.Bij, self.Aij, self.Cij, self.Fk);
                    condition = 0;
                    for num_plato = 1:self.etapas
                        nuevas_varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust) = self.varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust) + self.damping.*deltaXvar(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust);   
                        if any(nuevas_varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust-1) < 0)
                            nuevas_varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust - 1) = self.varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust-1).*exp((self.damping.*deltaXvar(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust-1))./(self.varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust-1)));
                        end
                        if num_plato == self.etapas && condition == 1
                            condition = 0;
                        end
                        Tbackup(num_plato) = self.perfil_t(num_plato);
                        Vbackup(num_plato) = self.perfil_v(num_plato);
                        Lbackup(num_plato) = self.perfil_l(num_plato);
                        self.perfil_t(num_plato) = nuevas_varX(1+(num_plato-1)*(2*(self.num_sust)+1)+2*self.num_sust);
                        if sum(nuevas_varX(1+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + self.num_sust-1)) > self.tol
                            self.perfil_v(num_plato) = sum(nuevas_varX(1+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + self.num_sust-1));
                        end
                        self.perfil_l(num_plato) = sum(nuevas_varX(1+self.num_sust+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + 2*self.num_sust-1));
                        
                        if (sum(self.Fk.^2)) < respaldovalI
                            respaldovalI = (sum(self.Fk.^2));
                            respaldobi = self.bsubi;
                            respaldodi = self.dsubi;
                            respaldoV = self.perfil_v;
                            respaldoL = self.perfil_l;
                            respaldovi = self.perfil_vi;
                            respaldoli = self.perfil_li;
                            respaldot = self.perfil_t;
                            respaldok =  self.perfil_k;
                            respaldoplatos = self.platos;
                            respaldoxi =  self.xsubj;
                            respaldoyi = self.ysubj;
                            if isempty(self.salidas) || ~isa(self.salidas, 'cell')
                                self.qc = -self.perfil_v(1)*(HVj(1)) - HLj(1)*self.perfil_l(1) + self.perfil_v(2)*HVj(2);
                            else
                                if self.salidas{1} == 1 && self.salidas{2} == 1 
                                    self.qc = -self.perfil_v(1)*(HVj(1)) - HLj(1)*self.perfil_l(1) - sj(1)*self.perfil_l(1)*HLj(1) + self.perfil_v(2)*HVj(2);
                                elseif self.salidas{1} == 1 && self.salidas{2} == 0
                                    self.qc = -self.perfil_v(1)*(HVj(1)) - HLj(1)*self.perfil_l(1) - sj(1)*self.perfil_l(1)*HLj(1) + self.perfil_v(2)*HVj(2);
                                else
                                    self.qc = -self.perfil_v(1)*(HVj(1)) - HLj(1)*self.perfil_l(1) + self.perfil_v(2)*HVj(2);
                                end
                            end
                            Fh = 0;
                            Vh1 = 0;
                            Uh = 0; %Salidas Laterales Liquidas 1
                            Wh = 0; %Salidas Laterales Vapor 1
                            Lh = 0; %Liquido que va al plato 2
                            Vh = 0; %Pudiera haber vapor fuga del plato 1
                            Vh2 = 0; %Vapor que entra al plato 1 del plato 2
                            LhEND = 0;
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
                                if ~isempty(self.platos(iterx).salidaV)
                                    if iterx == 1
                                        Vh1 = Vh1 +  self.perfil_v(1) * HVj(iterx);  % que no se sume
                                            %2 veces el destilado vapor
                                    else    
                                        Wh = Wh + self.platos(iterx).salidaV * HVj(iterx);
                                    end
                                end
                                if ~isempty(self.platos(iterx).salidaL)
                                    if iterx == self.etapas
                                        LhEND = LhEND  + self.platos(iterx).salidaL * HLj(num_plato);
                                        %que no se sume dos veces el fondo l�quido
                                    else                            
                                        Uh = Uh + self.platos(iterx).salidaL * HLj(num_plato);
                                    end
                                end
                            end
                            self.qb = Fh - Uh  - Wh - Vh1 - LhEND - sum(self.perfil_q) - self.qc;
                            respaldoqc= self.qc;
                            respaldoqb = self.qb;
                            self.respaldoiter = respaldandoiter;
                        end
                        
                        valI = sum(self.Fk.^2);
    %                     valI = valI + ((self.perfil_t(num_plato) - Tbackup(num_plato))/self.perfil_t(num_plato))^2;
    %                     if self.perfil_v(num_plato) > 0.00001 .* self.tol
    %                         valI = valI + ((self.perfil_v(num_plato) - Vbackup(num_plato))/self.perfil_v(num_plato))^2;
    %                     end
    %                     if self.perfil_l(num_plato) > 0.00001.* self.tol
    %                         valI = valI + ((self.perfil_l(num_plato) - Lbackup(num_plato))/self.perfil_l(num_plato))^2;
    %                     end
                        self.varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust) = nuevas_varX(1 + (num_plato - 1)*(2*self.num_sust + 1):1+(num_plato - 1)*(2*self.num_sust + 1)+2*self.num_sust);

                        self.perfil_vi(num_plato, :) = self.varX(1+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + self.num_sust-1)';
                        self.perfil_li(num_plato, :) = self.varX(1+self.num_sust+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + 2*self.num_sust-1)';
                        if num_plato == 1
                            if ~isempty(self.salidas) && isa(self.salidas, 'cell') 
                                if self.salidas{1} == 1 && self.salidas{2} == 1 && self.perfil_v(num_plato)/(self.perfil_v(num_plato)+self.perfil_l(num_plato)) < 1e-4
                                    self.ysubj(num_plato,:) = self.varX(1+(num_plato)*(2*(self.num_sust)+1):1+(num_plato)*(2*(self.num_sust)+1) + self.num_sust-1)'./sum(self.varX(1+(num_plato)*(2*(self.num_sust)+1):1+(num_plato)*(2*(self.num_sust)+1) + self.num_sust-1));
                                else
                                    self.ysubj(num_plato,:) = self.varX(1:self.num_sust)'./sum(self.varX(1:self.num_sust));
                                end
                            else
                                self.ysubj(num_plato,:) = self.varX(1:self.num_sust)'./sum(self.varX(1:self.num_sust));
                            end
                        else
                            self.ysubj(num_plato,:) = self.varX(1+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + self.num_sust-1)'./sum(self.varX(1+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + self.num_sust-1));
                        end
                        self.xsubj(num_plato,:) = self.varX(1+self.num_sust+(num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + 2*self.num_sust-1)'./sum(self.varX(1+self.num_sust + (num_plato-1)*(2*(self.num_sust)+1):1+(num_plato-1)*(2*(self.num_sust)+1) + 2*self.num_sust-1));                    
                        self.ysubj(num_plato,:) = self.ysubj(num_plato,:);
                        self.xsubj(num_plato,:) = self.xsubj(num_plato,:);
                        self.perfil_t(num_plato) = self.varX(1+(num_plato-1)*(2*(self.num_sust)+1)+2*self.num_sust);

                    end

                    for iiiiter = 1:self.etapas
                        if iiiiter == self.etapas
                            fprintf(1, ' %f \n', self.perfil_t(iiiiter));
                            break
                        end
                        if iiiiter == 1
                            fprintf(1, '\n %f ', self.perfil_t(iiiiter));
                        else
                            fprintf(1, ' %f ', self.perfil_t(iiiiter));
                        end

                    end

                    for iteru = 1:self.etapas

                        self.platos(iteru).setV(self.perfil_v(iteru));
                        self.platos(iteru).setL(self.perfil_l(iteru));
                        self.platos(iteru).setyi(self.ysubj(iteru,:));                    
                        self.platos(iteru).setxi(self.xsubj(iteru,:));
                        self.platos(iteru).K = self.perfil_k(:, iteru);
                    end
                    self.bsubi = self.xsubj(end,:).*(self.perfil_l(end));
    %                 self.perfil_t(1) = Backupt1;

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


                    for i = 1:self.etapas
                        self.platos(i).L = self.perfil_l(i);
                        self.platos(i).V = self.perfil_v(i);
                        self.platos(i).y_i = self.ysubj(i,:);
                        self.platos(i).x_i = self.xsubj(i,:);
                        self.platos(i).v_i = self.perfil_vi(i,:);
                        self.platos(i).l_i = self.perfil_li(i,:);
                    end
                    self.bflujo = sum(self.perfil_li(end,:));
                    if isempty(self.salidas) 
                        self.dflujo = sum(self.perfil_vi(1,:)); 
                    else
                        if isa(self.salidas, 'cell')
                            if self.salidas{1} == 1 && self.salidas{2} == 1
                                self.dflujo =  sum(self.perfil_vi(1,:)) + self.F - self.UmW - sum(self.perfil_li(end,:));
                            end
                        else
                            self.dflujo = sum(self.perfil_vi(1,:)); 
                        end
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
                    self.respaldovalI(self.actualiter) = valI;

                    if valI < self.tol
                        self.laststep = logical(1);
                    end
                
                if isempty(self.salidas) || ~isa(self.salidas, 'cell')
                    mezclaX = Mezcla(self.sust, self.xsubj(1,:), self.alimentaciones{2}.mezcla.kij);
                    [~, yy ] = self.MEdE.BubbleT(mezclaX, self.perfil_p(1), self.alimenciones{2}.mezcla.kij);
                    [Tb, ~ ] = self.MEdE.BubbleT(mezclaX, self.perfil_p(1), self.alimenciones{2}.mezcla.kij, yy);
                    mezclaY = mezclaX;
                    mezclaY.conc = self.ysubj(1,:);
                    for ittwee = 1:self.num_sust
                        try 
                            cp = self.sust(i).cp_gi{1};
                        catch ME
                            error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                        end
                        deltaHV1 = integral(@(t) cp(t), 273.15, Tb);
                        HgiV1 = deltaHV1*mezclaY.conc(i) + Href(ittwee) * mezclaY.conc(i);                
                    end
                    HVj(1) = HgiV1 - self.MEdE.entalpia( Tb, self.perfil_p(1), mezclaY, 'vap') + self.MEdE.entalpia( 273.15, 101.325, mezclaY, 'liq');
                    self.qc = -self.perfil_v(1)*(HVj(1)) - HLj(1)*self.perfil_l(1) + self.perfil_v(2)*HVj(2);
                else
                    if self.salidas{1} == 1 && self.salidas{2} == 1 
                        self.qc = -self.perfil_v(1)*(HVj(1)) - HLj(1)*self.perfil_l(1) - sj(1)*self.perfil_l(1)*HLj(1) + self.perfil_v(2)*HVj(2);
                    elseif self.salidas{1} == 1 && self.salidas{2} == 0
                        mezclaX = Mezcla(self.sust, self.xsubj(1,:), self.alimentaciones{2}.mezcla.kij);
                        [~, yy ] = self.MEdE.BubbleT(mezclaX, self.perfil_p(1), self.alimenciones{2}.mezcla.kij);
                        [Tb, ~ ] = self.MEdE.BubbleT(mezclaX, self.perfil_p(1), self.alimenciones{2}.mezcla.kij, yy);
                        mezclaY = mezclaX;
                        mezclaY.conc = self.ysubj(1,:);
                        for ittwee = 1:self.num_sust
                            try 
                                cp = self.sust(i).cp_gi{1};
                            catch ME
                                error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                            end
                            deltaHV1 = integral(@(t) cp(t), 273.15, Tb);
                            HgiV1 = deltaHV1*mezclaY.conc(i) + Href(ittwee) * mezclaY.conc(i);                
                        end
                        HVj(1) = HgiV1 - self.MEdE.entalpia( Tb, self.perfil_p(1), mezclaY, 'vap') + self.MEdE.entalpia( 273.15, 101.325, mezclaY, 'liq');
                        self.qc = -self.perfil_v(1)*(HVj(1)) - HLj(1)*self.perfil_l(1) - sj(1)*self.perfil_l(1)*HLj(1) + self.perfil_v(2)*HVj(2);
                    else
                         mezclaX = Mezcla(self.sust, self.xsubj(1,:), self.alimentaciones{2}.mezcla.kij);
                        [~, yy ] = self.MEdE.BubbleT(mezclaX, self.perfil_p(1), self.alimenciones{2}.mezcla.kij);
                        [Tb, ~ ] = self.MEdE.BubbleT(mezclaX, self.perfil_p(1), self.alimenciones{2}.mezcla.kij, yy);
                        mezclaY = mezclaX;
                        mezclaY.conc = self.ysubj(1,:);
                        for ittwee = 1:self.num_sust
                            try 
                                cp = self.sust(i).cp_gi{1};
                            catch ME
                                error('Sustancia.cp_gi: Un compuesto no tiene un function_handle de cp_gi. Agregue uno a la clase Sustancia.m correspondiente');
                            end
                            deltaHV1 = integral(@(t) cp(t), 273.15, Tb);
                            HgiV1 = deltaHV1*mezclaY.conc(i) + Href(ittwee) * mezclaY.conc(i);                
                        end
                        HVj(1) = HgiV1 - self.MEdE.entalpia( Tb, self.perfil_p(1), mezclaY, 'vap') + self.MEdE.entalpia( 273.15, 101.325, mezclaY, 'liq');
                        self.qc = -self.perfil_v(1)*(HVj(1)) - HLj(1)*self.perfil_l(1) + self.perfil_v(2)*HVj(2);
                    end

                end
                Fh = 0;
                Vh1 = 0;
                Uh = 0; %Salidas Laterales Liquidas 1
                Wh = 0; %Salidas Laterales Vapor 1
                Lh = 0; %Liquido que va al plato 2
                Vh = 0; %Pudiera haber vapor fuga del plato 1
                Vh2 = 0; %Vapor que entra al plato 1 del plato 2
                LhEND = 0;
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
                    if ~isempty(self.platos(iterx).salidaV)
                        if iterx == 1
                            Vh1 = Vh1 +  self.perfil_v(1) * HVj(iterx);  % que no se sume
                                %2 veces el destilado vapor
                        else    
                            Wh = Wh + self.platos(iterx).salidaV * HVj(iterx);
                        end
                    end
                    if ~isempty(self.platos(iterx).salidaL)
                        if iterx == self.etapas
                            LhEND = LhEND  + self.platos(iterx).salidaL * HLj(num_plato);
                            %que no se sume dos veces el fondo l�quido
                        else                            
                            Uh = Uh + self.platos(iterx).salidaL * HLj(num_plato);
                        end
                    end
                end
                self.qb = Fh - Uh  - Wh - Vh1 - LhEND - sum(self.perfil_q) - self.qc;
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
                self.bsubi = respaldobi;
                self.dsubi = respaldodi;
                self.dflujo = sum(respaldodi);
                self.bflujo = sum(respaldobi);
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
    end
end