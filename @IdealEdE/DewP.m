function [P, x, K, flag, val, iter] = DewP(self, T, mezcla, sup_x, varargin)
%[P, x, K, flag, val] = DewP(T, mezcla, sup_x, tolerancia, max_iter, P_semilla)
%
%Calcula la presion de rocío (P) de una mezcla en fase vapor así como la
%composición de la primera gota de líquido (x)
%
%Parámetros
%
%Input:
%T: Temperatura de la mezcla en K
%mezcla: Clase Mezcla.m que contiene las propiedades de la mezcla como
%composición, y clases Sustancia.m que definen las propiedades de sustancia
%pura
%sup_x: Vector fila que contiene la suposición inicial de la primera gota
%
%Output:
%
%P: Presión de rocío
%x: composición de la primera gota de rocío (vector fila)
%K: Volatilidades de los componentes (vector fila)
%
%flag: flag del método de resolución numérico:
%     1  Encontró un cero X.
%     2  Cálculos convergieron a la región de una fase
%     3  Se alcanzó el máximo número de iteraciones
%val: Salida de la función handle de fzero la cual hay que minimizar.

    Gas = mezcla;
    num_sust = Gas.num_sust;
    Comp = mezcla.comp;
    kij = mezcla.kij;
    y = mezcla.conc;
    if (nargin < 4) || (isempty(sup_x))
        sup_x = zeros(1, num_sust);
        sup_x(end) = 1;
    end
    Liquido = Mezcla(Comp, sup_x, kij);
    if nargin < 7
        P_semilla = zeros(1,num_sust);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Estimado inicial de T (usa antoine si puede o base de datos DIPPR
        for i=1:num_sust
            try
                psat = Comp(i).psat{1};
                P_semilla(i) = psat(T);
            catch
                P_semilla(i) = psat(T);
            end
        end    

        y = y(P_semilla ~= 0);
        y = y./sum(y);
        P_semilla = sum(P_semilla(P_semilla ~= 0).*y(P_semilla~=0)); %Si algún compuesto no se calculó el Tsat lo ignora
        P = P_semilla;
        P(P == 0) = 1000;
    else
        P = varargin{3};
    end
    dP = 0.04 .* P;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    y = mezcla.conc;
    if nargin < 5 
        tol = 5e-6;
    else
        tol = varargin{1};
        if isempty(tol)
            tol = 5e-6;
        end
    end
    if nargin < 6
        max_iter = 1000;
    else
        max_iter = varargin{2};
        if isempty(max_iter)
            max_iter = 1000;
        end
    end
    valI = 1;
    error_ = 1;
    iter = 0;
    x_calc_ant = sup_x;    
    K = zeros(1, num_sust);
    x_calc = zeros(1, num_sust); 
    while ((valI > tol) && (iter < max_iter) && (abs(error_) > tol)) %Itera mientras no cumpla la tolerancia o
                                            % no sea superado el número de iteraciones
        iter = iter + 1;
        fiG = self.fug(T,P,Gas,'vap'); %fugacidades de vapor
        fiL = self.fug(T,P,Liquido,'liq');  %Obtiene la fugacidad de líquido
        for i = 1:num_sust %Coeficientes de distribución de líquido y gas (volatilidades) y composiciones calculadas de líquido
            K(i) = fiL(i)/fiG(i);
            x_calc(i) = y(i)/K(i);
        end
        valI = 0;
        for i = 1:num_sust
            valI = valI + abs(x_calc_ant(i) - x_calc(i)); %Compara la composición de líquido calculada de 2 iteraciones sucesivas
        end
        if iter > 1
            error_ant = error_;
        end
        error_ = sum(x_calc) - 1;
        if abs(error_) > tol
            if iter > 1   
                if sign(error_) ~= sign(error_ant)
                    Psup = (P - P_ant)./(error_-error_ant).*(-error_ant)+ P_ant;
                    dP = dP ./ 4;
                else
                    Psup = P - sign(error_)*dP;
                end
            else
                Psup = P - sign(error_)*dP;
            end
            error_ant = error_;
            P_ant = P;
            P = Psup;    
        end
        x_calc_ant = x_calc;
        Liquido.conc = x_calc./sum(x_calc);
    end
%Si se alcanzó la región de una fase, donde no puede haber equilibrio, muestra un warning
    single_phase = 0;
    if (all(abs(K - 1)) < 1e-2)
        single_phase = 1;
    end
    if single_phase == 1
        warning('MATLAB:IdealEdE. Se alcanzó la región de solo una fase para la mezcla ');
        flag = 2;
    else        
        flag = 1;
    end
    if (iter == max_iter) %Si el máximo número de iteraciones se alcanzó, muestra un warning, 
        s = sprintf('%f',valI);
        warning('MATLAB:IdealEdE', ['Error de convergencia en DewPoint. El valor final de la función handle fue: ' s '.']);
        flag = 3;
    end
    if (P < 0)
        warning('Matlab:IdealEdE. Error de cálculo de presión de rocío. El valor es negativo. Región de 1 fase');
    end
    
    val = sum(x_calc) - 1; %Calcula la suma Xi - 1 (debe hacerse cero) y regresa el resultado
    x = Liquido.conc; %
end