classdef PRGEdE < EdECubica
    %PRGEdE Constructor de la ecuación de estado de Peng-Robinson-Gasem
	%
	%PRGEdE(sigma, epsilon, OMEGA, PSI) es una clase que calcula la ecuación 
    %de estado de Peng-Robinson-Gasem con las variables sigma, epsilon, OMEGA, PSI
    %de acuerdo a la forma genérica planteada en Smith y Van Ness (referencia)
	%
    %Todos los argumentos son opcionales, se utilizarán los default de la 
    %referencia si no se proveen. Se puede llamar con argumento vacío si se
    %desea dejar el valor por defecto de PR, como argumento "[]"
    %
    %Opcionalmente pueden proveerse function handle para sustituir los alfa/m 
    %default de PR que se muestran a continuación. 
    %
    % A = 2.00
    % B = 0.836
    % m = @(omega) 0.134 + 0.508.*omega - 0.0467.*omega.^2
    % alfa = @(t, tcri, m) exp[(A + B.*(t./tcri)).*(1 - (t./tcri).^m;)]
        %
    %Gasem = PRGEdE([],[],[],[], alfa, m);
    %
    %Sustituirá los function handle @PRGEdE.alfa_fun.m y @PRGEdE.m_fun.m por 
    %los definitidos alfa(t, tcri, m) y m(omega). 
    %
    %Si los valores provistos no son function handle utilizará el default.
    %
    %Referencia:
    %
    %Gasem, Gao, Pan & Robinson, Fluid Phase Equilibria, 181, 113-125 (2001).
    %
    %Figueira, Freddy (2005). "Desarrollo de ecuaciones de estado del tipo Van
    %der Waals para fluidos puros polares y no polares". Universidad Simón
    %Bolívar. Venezuela
    %
    %Smith, Van Ness, Abbott. Introduction to Chemical Engineering Thermodynamics. 
    % 7th edition 
    %
    %   Luis Jesús Díaz Manzo
    
    properties
        id = 'Peng-Robinson-Gasem';
        propied = [1 + sqrt(2), 1 - sqrt(2), 0.07780, 0.45724];
        %parámetros sigma, épsilon, OMEGA y PSI de ecuación cúbica
        %generalizada  (ver referencia)
        alfa_inp = false;
        m_inp = false;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %coeficientes de correlación entre alfar/Tr con Pr/Tr para 
        %toda la región de saturación (propuesta por Carreira
        %y Rodriguez y referenciada por Figueira)
        coefE = false;
        coefH0 = false; %No previstos en la bibliografía
        coefH19 = false;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		derivada_alfa = false;
    end
    
    methods
        function self = PRGEdE(varargin)
            if nargin > 0
                sigma = varargin{1};
                if ~isempty(sigma)
                    self.propied(1) = sigma;
                end
            end
            if nargin > 1                
                epsilon = varargin{2};
                if ~isempty(epsilon)
                    self.propied(2) = epsilon;
                end
            end
            if nargin > 2                
                OMEGA = varargin{3};
                if ~isempty(OMEGA)
                    self.propied(3) = OMEGA;
                end 
            end
            if nargin > 3
                PSI = varargin{4};
                if ~isempty(PSI)
                    self.propied(4) = PSI;
                end
            end
            if nargin > 4
                alfa = varargin{5};
                if  ~isempty(alfa)
                    if isa(alfa, 'function_handle')
                        self.alfa_inp = alfa;
                    else
                        warning(['El input alfa(w) = "' alfa '" no es de clase function handle.']);
                    end
                end
            end
            if nargin > 5
                m = varargin{6};
                if ~isempty(m)
                    if isa(m, 'function_handle')
                        self.m_inp = m;
                    else
                        warning(['El input m(w) = "' m '" no es de clase function handle.']);
                    end
                end
            end
        end
        function [ a_factor ] = a_fun(self, T, Comp)
            if ~islogical(self.alfa_inp)
                if ~islogical(self.m_inp)
                    a_factor = self.a_handle(T, Comp, self.alfa_inp, self.m_inp);
                else
                    a_factor = self.a_handle(T, Comp, self.alfa_inp);
                end
            elseif ~islogical(self.m_inp)
                a_factor = self.a_handle(T, Comp, [], self.m_inp);
            else 
                a_factor = self.a_handle(T, Comp);
            end
        end
    end
    
end

