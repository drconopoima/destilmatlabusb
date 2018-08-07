classdef PREdE < EdECubica
    %PREdE Constructor de la ecuaci�n de estado de Peng-Robinson
	%
	%PREdE(sigma, epsilon, OMEGA, PSI) es una clase que calcula la ecuaci�n 
    %de estado de Peng-Robinson con las variables sigma, epsilon, OMEGA, PSI
    %de acuerdo a la forma gen�rica planteada en Smith y Van Ness (referencia)
	%
    %Todos los argumentos son opcionales, se utilizar�n los default de la 
    %referencia si no se proveen. Se puede llamar con argumento vac�o si se
    %desea dejar el valor por defecto de PR, como argumento "[]"
    %
    %Opcionalmente pueden proveerse function handle para sustituir los alfa/m 
    %default de PR que se muestran a continuaci�n. 
    %
    %alfa = @(t, tcri, m) (1 + (m.*(1 - (t./tcri).^0.5))).^2;
    %m = @(omega) 0.37464 + 1.54226.*omega - 0.26992.*omega.^2);
    %
    %Peng = PREdE([],[],[],[], alfa, m);
    %
    %Sustituir� los function handle @PREdE.alfa_fun.m y @PREdE.m_fun.m por 
    %los definitidos alfa(t, tcri, m) y m(omega). 
    %
    %Si los valores provistos no son function handle utilizar� el default.
    %
    %Referencia:
    %
    %Figueira, Freddy (2005). "Desarrollo de ecuaciones de estado del tipo Van
    %der Waals para fluidos puros polares y no polares". Universidad Sim�n
    %Bol�var. Venezuela
    %
    %Smith, Van Ness, Abbott. Introduction to Chemical Engineering Thermodynamics. 
    % 7th edition 
    %
    %   Luis Jes�s D�az Manzo
    
    properties
        id = 'Peng-Robinson';
        propied = [1 + sqrt(2), 1 - sqrt(2), 0.07780, 0.45724];
        %par�metros sigma, �psilon, OMEGA y PSI de ecuaci�n c�bica
        %generalizada  (ver referencia)
        alfa_inp = false;
        m_inp = false;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %coeficientes de correlaci�n entre alfar/Tr con Pr/Tr para 
        %toda la regi�n de saturaci�n (propuesta por Carreira
        %y Rodriguez y referenciada por Figueira)
        coefE = 0.00041;
        coefH0 = 2.9803582e-1;
        coefH19 = [1.5003698e-2, -4.7527103e-3, 8.036716e-4,...
            -8.9548695e-5, 6.8691611e-6, -3.6067317e-7, 1.2409205e-8, ...
            -2.5222671e-10, 2.2955503e-12];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		derivada_alfa = false;
    end
    
    methods
        function self = PREdE(varargin)
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

