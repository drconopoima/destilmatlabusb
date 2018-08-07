function [P, y, K, flag, val, iter] = BubbleP(self, T, mezcla, sup_y, varargin)
%[P, y, K, flag, val] = BubbleP(T, mezcla, sup_y, tolerancia, max_iter, P_semilla)
%
%Calcula la presión de burbuja (P) de una mezcla en fase líquida así como la
%composición de la primera burbuja de vapor (y)
%
%Parámetros
%
%Input:
%T: Temperatura de la mezcla en K
%mezcla: Clase Mezcla.m que contiene las propiedades de la mezcla como
%composición, x clases Sustancia.m que definen las propiedades de sustancia
%pura
%sup_y: Vector fila que contiene la suposición inicial de la primera gota
%
%Output:
%
%P: Presión de burbuja
%y: composición de la primera gota de rocío (vector fila)
%K: Volatilidades de los componentes (vector fila)
%
%flag: flag del método de resolución numérico:
%     1  Encontró un cero X.
%     2  Cálculos convergieron a la región de una fase
%     3  Se alcanzó el máximo número de iteraciones
%val: Salida de la función handle de fzero la cual hax que minimizar.
%iter: Devuelve el número de iteraciones que fueron requeridas
    
    Liquido = mezcla;
    num_sust = Liquido.num_sust;
    Comp = mezcla.comp;
    kij = mezcla.kij;
    x = mezcla.conc;
    if (nargin < 4) || (isempty(sup_y))
        sup_y = zeros(1, num_sust);
        sup_y(1) = 1;
    end
    Gas = Mezcla(Comp, sup_y, kij);
    if nargin < 7
        P_semilla = zeros(1,num_sust);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Estimado inicial de T (usa antoine si puede o base de datos DIPPR
        for i=1:num_sust
            x11 = x;
            if T > Comp(i).tcri
                if x == x11
                    x(i) = 0.06*x(i);
                end
            end
            x = x./sum(x);
            try
                psat = Comp(i).psat{1};
                P_semilla(i) = psat(T);
            catch
                P_semilla(i) = 0;
            end
        end    
        x = x11 / sum(x);
        x = x(P_semilla ~= 0);
        x = x./sum(x);
        P_semilla = sum(P_semilla(P_semilla ~= 0).*x(P_semilla~=0)); %Si algún compuesto no se calculó el Psat lo ignora
        P = P_semilla;
        P(P==0) = 1000;
    else
        P = varargin{3};
    end
    dP = 0.04 .* P;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    x = mezcla.conc;
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
    y_calc_ant = sup_y;    
    K = zeros(1, num_sust);
    y_calc = zeros(1, num_sust); 
    while ((valI > tol) && (iter < max_iter) && (abs(error_) > tol))  %Itera mientras no cumpla la tolerancia o
                                            % no sea superado el número de iteraciones
        iter = iter + 1;
        fiG = self.fug(T,P,Gas,'vap'); %fugacidades de vapor
        fiL = self.fug(T,P,Liquido,'liq');  %Obtiene la fugacidad de líquido
        for i = 1:num_sust %Coeficientes de distribución de líquido y gas (volatilidades) y composiciones calculadas de líquido
            K(i) = fiL(i)/fiG(i);
            y_calc(i) = K(i)*x(i);
        end
        valI = 0;
        for i = 1:num_sust
            valI = valI + abs(y_calc_ant(i) - y_calc(i)); %Compara la composición de líquido calculada de 2 iteraciones sucesivas
        end
        if iter > 1
            error_ant = error_;
        end
        error_ = sum(y_calc) - 1;
        if abs(error_) > tol
            if iter > 1   
                if sign(error_) ~= sign(error_ant)
                    Psup = (P - P_ant)./(error_-error_ant).*(-error_ant)+ P_ant;
                    dP = dP ./ 4;
                else
                    Psup = P + sign(error_)*dP;
                end
            else
                Psup = P + sign(error_)*dP;
            end
            error_ant = error_;
            P_ant = P;
            P = Psup;    
        end
        y_calc_ant = y_calc;
        Gas.conc = y_calc./sum(y_calc);
    end
    
%Si se alcanzó la región de una fase, donde no puede haber equilibrio, muestra un warning
    single_phase = 0;
    if (all(abs(K - 1)) < 1e-2)
        single_phase = 1;
    end
    if single_phase == 1
        warning('MATLAB:IdealEdE. Se alcanzó la región de solo una fase para la mezcla');
        flag = 2;
    else        
        flag = 1;
    end
    if (iter == max_iter) %Si el máximo número de iteraciones se alcanzó, muestra un warning, 
        s = sprintf('%f',valI);
        warning('MATLAB:IdealEdE', ['Error de convergencia en DewPoint. El valor final de la función handle fue: ' s '.']);
        flag = 3;
    end
    val = sum(y_calc) - 1; %Calcula la suma Xi - 1 (debe hacerse cero) x regresa el resultado
    y = Gas.conc;
end