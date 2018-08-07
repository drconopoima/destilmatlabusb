classdef SRKEdE < EdECubica
    %SRKEdE Constructor de la ecuación de estado de Soave-Redlich-Kwong
	%
	%SRKEdE(sigma, epsilon, OMEGA, PSI) es una clase que calcula la ecuación 
    %de estado de Soave-Redlich-Kwong con las variables sigma, epsilon, OMEGA, PSI
    %de acuerdo a la forma genérica planteada en Smith y Van Ness (referencia)
	%
    %Todos los argumentos son opcionales, se utilizarán los default de la 
    %referencia si no se proveen. Se puede llamar con argumento vacío si se
    %desea dejar el valor por defecto de SRK, como argumento "[]"
    %
    %Opcionalmente pueden proveerse function handle para sustituir los alfa/m 
    %default de SRK que se muestran a continuación. 
    %
    %alfa = @(t, tcri, m) (1 + (m.*(1 - (t./tcri).^0.5))).^2;
    %m = @(omega) 0.480 + 1.574.*omega - 0.176.*omega.^2);
    %
    %SRKEdE([],[],[],[], alfa, m);
    %
    %Sustituirá¡ los function handle @SRKEdE.alfa_fun.m y @SRKEdE.m_fun.m por 
    %los definitidos alfa(t, tcri, m) y m(omega). 
    %
    %Si los valores provistos no son function handle utilizará el default.
    %
    %Referencia:
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
        id = 'Soave-Redlich-Kwong';
        propied = [1, 0, 0.08664, 0.42748];
        %parÃ¡metros sigma, Ã©psilon, OMEGA y PSI de ecuaciÃ³n cÃºbica
        %generalizada  (ver referencia)
        alfa_inp = false;
        m_inp = false;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %coeficientes de correlaciÃ³n entre alfar/Tr con Pr/Tr para 
        %toda la regiÃ³n de saturaciÃ³n (propuesta por Carreira
        %y Rodriguez y referenciada por Figueira)
        coefE = 1.014;
        coefH0 = 3.2748e-1;
        coefH19 = [3.4376954e-1, 1.0596403e-2, -3.7538497e-3, 6.3257197e-4,...
            -6.6481e-5, 4.5517956e-6, -1.9796921e-7, 4.9748592e-9, ...
            -5.5024228e-11];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		derivada_alfa = false;
    end
    
    methods
        function self = SRKEdE(varargin)
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