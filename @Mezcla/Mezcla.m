classdef Mezcla < handle
    %Mezcla Objeto que maneja la agrupación de clases Sustancia.m y guarda
    %propiedades como el número de sustancias, la composición en fracción
    %molar y másica y los parámetros de interacción binaria entre ellas.
    %
    %Mezcla(-opc- componentes, -opc- concentraciones, -opc- kij)
    %
    %   Luis Jesús Díaz Manzo
    
    properties
        comp = false;
        num_sust = false;
        conc = false;
        kij = false;
        kbinario = false;
        kmatrix = false;
    end
    
    methods
        function self = Mezcla(componentes, concent, interaccion_bin)
            if nargin > 0
                tamano = size(componentes);
                if tamano(1) > tamano(2)
                    componentes = componentes';
                end
                if isa(componentes, 'cell')
                    com = Sustancia.empty(0, length(componentes));
                    for i = 1:length(componentes)
                        com(i) = componentes{i};
                    end
                    componentes = com;
                end
                if isa(componentes, 'Sustancia')
                    self.comp = componentes;
                    tamano = size(self.comp);
                    self.num_sust = tamano(2);
                    self.kij = zeros(1, factorial(self.num_sust)./(2.*factorial(self.num_sust - 2)));
                    tamano = size(self.kij);
                    self.kij = cell(tamano(1), tamano(2));
                    self.kbinario = self.kij;
                    i = 0;
                    while i < tamano(2)
                        for j = i+1: self.num_sust
                            for k = j+1 : self.num_sust
                                i = i+1;
                                n = strcat('k', int2str(j), int2str(k));
                                self.kbinario{i} =  n;
                                self.kij{i} = 0; 
                            end
                        end
                    end
                end
            end
            if nargin > 1
                tamano = size(concent);
                if tamano(1) > tamano(2)
                    concent = abs(concent'); %vector fila de concentraciones positivas
                else 
                    concent = abs(concent);
                end
                tamano = size(concent);
                suma = sum(concent);
                if tamano(2) == self.num_sust - 1 %Si se provee las primeras N-1 concentraciones independientes calcula la dependiente
                    if suma > 100
                        warning('Concentraciones no son fracciones molares en base unitaria. Se supone un(varios) componente(s) no presente(s)');
                        self.conc = concent./suma;
                        self.conc(end + 1) = 0;
                    elseif suma > 1
                        self.conc = concent./suma;
                        self.conc(end + 1) = 0;
                    else
                        self.conc = concent;
                        self.conc(end + 1) = 1 - suma;
                    end
                elseif tamano(2) == self.num_sust %Si se proveen las N concentraciones
                    if (suma - 1) > 1e-6 
                        warning('Las concentraciones deben proveerse en fracciones molares en base unitaria. Se procede a su conversión');
                        self.conc = concent./suma;
                    else
                        self.conc = concent;
                    end
                else 
                    warning('Mala especificación de las concentraciones molares. Sobran o faltan componentes');
                end
            end
            if nargin > 2
                for i = 1:length(interaccion_bin)
                    self.kij{i} = interaccion_bin(i);
                end
            end
            self.kij = cell2mat(self.kij);
            self.kmatrix = self.matrixkij();
        end
        function matrixk = matrixkij(self)
            ki = self.kij;
            kbin = self.kbinario;
            n_sust = self.num_sust;
            matrixk = zeros(n_sust, n_sust);
            for i = 1: n_sust
                for j = 1:n_sust
                    if i > j
                        continue
                    end
                    k = 'k';
                    k = strcat(k, int2str(i));
                    k = strcat(k, int2str(j));
                    for l = 1:length(ki)
                        try
                            if strcmpi(kbin{l},k)
                                matrixk(i, j) = ki(l);
                                matrixk(j, i) = ki(l);
                            end
                        catch
                        end
                    end
                end
            end
        end
        function self = definek(self, interaccion_bin)
            self.kij = zeros(1, factorial(self.num_sust)./(2.*factorial(self.num_sust - 2)));
            tamano = size(self.kij);
            self.kij = cell(tamano(1), tamano(2));
            self.kbinario = self.kij;
            i = 0;
            while i < tamano(2)
                for j = i+1: self.num_sust
                    for k = j+1 : self.num_sust
                        i = i+1;
                        n = strcat('k', int2str(j), int2str(k));
                        self.kbinario{i} =  n;
                        self.kij{i} = 0; 
                    end
                end
            end
            for i = 1:length(interaccion_bin)
                self.kij{i} = interaccion_bin(i);
            end
            kijaux = zeros(1, self.num_sust);
            for i = 1:length(self.kij)
                kijaux(1,i) = self.kij{i};
            end
            self.kij = kijaux;
            self.kmatrix = self.matrixkij();
        end
    end
    
end

