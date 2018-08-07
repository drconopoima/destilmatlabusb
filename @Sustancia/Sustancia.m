classdef Sustancia < handle
    %Sustancia es una clase que contiene propiedades de sustancias puras
    %   Luis Jesús Díaz Manzo
    %   %
    %   Referencias:
    %
    %   % R.L. Rowley, W.V. Wilding, J.L. Oscarson, Y. Yang, N.A. Zundel, T.E. Daubert, R.P. Danner, (2007)
    %   DIPPR(r) Data Compilation of Pure Chemical Properties, Design, Institute for Physical Properties. AIChE., New York. 
    
    properties        
        id = false; 
        molmass = false;
        unit_molmass = 'kg/kg-mol';
        tcri = false;
        unit_tcri = 'K';
        pcri = false;
        unit_pcri = 'kPa';
        w_acent = false;        
        psat = false;
        unit_psat = 'kPa';
        tsat = false; %function handle @tsat(P, T) que con fzero puede hallar la temperatura de saturación
        unit_tsat = 'K';
        href = 0;  %entalpía de referencia a 0.01°C (273.01 K) default es 0
        cp_gi = false;
        Hvap = false;
        antoine = false;
        formul = false;
        indice = false;
    end
    
    methods
        function self = Sustancia(identificador, varargin)
            dbstop if error
            if nargin>0 && ~isempty(identificador)
               self.id = identificador;
            end
            if nargin > 1 && ~isempty(varargin{1})
                pesomolar = varargin{1};
                if isempty(pesomolar) == 0
                    self.molmass = pesomolar;
                end
            end
            if nargin > 2   
                tc = varargin{2};
                if isempty(tc) == 0
                    self.tcri = tc;
                end
            end
            if nargin > 3
                pc = varargin{3};
                if isempty(pc) == 0
                    self.pcri = pc;
                end
            end
            if nargin > 4
            omega = varargin{4};
                if isempty(omega) == 0
                    self.w_acent = omega;
                end
            end
            if nargin > 5         
                pvap = varargin{5};    
                if isempty(pvap) == 0
                    self.psat = pvap;
                end
            end
            if nargin > 6
                tvap = varargin{6};
                if isempty(tvap) == 0
                    self.tsat = tvap;
                end
            end
            if nargin > 7   
                hf = varargin{7};
                if isempty(hf) == 0
                    self.href = hf;
                end
            end
            if nargin > 8
                cp_gasi = varargin{8};
                if isempty(cp_gasi) == 0
                    self.cp_gi = cp_gasi;
                end
            end
            if nargin > 9
                form = varargin{9};
                if isempty(form) == 0
                    self.formul = form;
                end
            end
            if nargin < 12 && nargin > 0
                try
                    if islogical(self.indice)
                        propiedades = get_properties(identificador);
                        if islogical(self.formul)
                            self.formul = propiedades{1};
                        end
                        if islogical(self.molmass)
                            self.molmass = propiedades{2};
                        end
                        if islogical(self.tcri)
                            self.tcri = propiedades{3};
                        end
                        if islogical(self.pcri)
                            self.pcri = propiedades{4};
                        end
                        if islogical(self.w_acent)
                            self.w_acent = propiedades{5};
                        end
                        self.indice = propiedades{6};
                        try 
                            if (islogical(self.psat)||islogical(self.tsat))...
                                    || ((~isa(self.psat, 'function_handle')...
                                    && ~isa(self.psat, 'double')) &&...
                                    (~isa(self.tsat, 'function_handle')...
                                    && ~isa(self.tsat, 'double')))
                                salida = get_vapour_pressure(identificador);
                                self.psat = {salida{1}, salida{3}, salida{4}};
                                self.tsat = {salida{2}, salida{5}, salida{6}};
                            end
                        catch
                            warning('El compuesto no está presente en la base de datos DIPPR. Puede agregarlo con funciones GET_INDICES.m y base de datos propiedades.dat Y al archivo vapourpressures.dat')
                        end
                        try 
                            if islogical(self.cp_gi) || (~isa(self.cp_gi, 'function_handle')...
                                    && ~isa(self.cp_gi, 'double'))
                                cpgasideal = get_ig_cp(identificador);
                                self.cp_gi = cpgasideal;
                            end
                        catch
                            warning('El compuesto no está presente en la base de datos de cp de gas ideal. \nPuede agregarlo manualmente o con funciones GET_INDICES.m y base de datos properties.dat\nY al archivo cpgasideal.dat')
                        end
                        try
                            if islogical(self.Hvap) || (~isa(self.Hvap, 'function_handle')...
                                    && ~isa(self.Hvap, 'double'))
                                latent_heat = get_latent_heat(identificador);
                                self.Hvap = latent_heat;
                            end
                        catch
                            warning('El compuesto no está presente en la base de datos de cp de gas ideal. \nPuede agregarlo manualmente o con funciones GET_INDICES.m y base de datos properties.dat\nY al archivo cpgasideal.dat')
                        end    
                    end
                catch ME
                    exception = ME;
                    msgString = getReport(exception);
                    warning('%sEl compuesto no está presente en la base de datos. \nPuede agregarlo con funciones GET_INDICES.m y base de datos properties.dat \nO puede definir manualmente la sustancia.\n', msgString)
                end  
            end
        end
        function salida = reset(self)
            self.id = false; 
            self.molmass = false;
            self.unit_molmass = 'kg/kg-mol';
            self.tcri = false;
            self.unit_tcri = 'K';
            self.pcri = false;
            self.unit_pcri = 'kPa';
            self.w_acent = false;        
            self.psat = false;
            self.unit_psat = 'kPa';
            self.tsat = false; %function handle @tsat(P, T) que con fzero puede hallar la temperatura de saturación
            self.unit_tsat = 'K';
            self.href = 0;  %entalpía de referencia a 0.01°C (273.01 K) default es 0
            self.cp_gi = false;
            self.Hvap = false;
            self.antoine = false;
            self.formul = false;
            self.indice = false;
            salida = self;
        end
    end
    
    
end

