function [T, x, K, flag, val, out_iter] = DewT(self, P, mezcla, sup_x, varargin)
%[T, x, K, flag, val] = DewT(P, mezcla, sup_x, tolerancia, max_iter, T_semilla)
%
%Calcula la temperatura de rocío (T) de una mezcla en fase vapor así como la
%composición de la primera gota de líquido (x)
%
%Parámetros
%
%Input:
%P: Presión de la mezcla en kPa
%mezcla: Clase Mezcla.m que contiene las propiedades de la mezcla como
%composición, y clases Sustancia.m que definen las propiedades de sustancia
%pura
%sup_x: Vector fila que contiene la suposición inicial de la primera gota
%
%Output:
%
%T: Temperatura de rocío
%x: composición de la primera gota de rocío (vector fila)
%K: Volatilidades de los componentes (vector fila)
%
%flag: flag del método de resolución numérico:
%     1  Encontró un cero X.
%     2  Cálculos convergieron a la región de una fase
%     3  Se alcanzó el máximo número de iteraciones
%val: Salida de la función handle de fzero la cual hay que minimizar.
%iter: Devuelve el número de iteraciones que fueron requeridas

    Gas = mezcla;
    num_sust = Gas.num_sust;
    Comp = mezcla.comp;
    kij = mezcla.kij;
    y = mezcla.conc;
    Ttotal = 10;
    single_phase = 0;
    single_phase_ant = 0;
    t_semilla = 0;
    out_iter = 0;
    if nargin > 3 && ~isempty(sup_x)
        Liquido = Mezcla(Comp, sup_x, kij);
    else
        Liquido = Mezcla(Comp, y, kij);
    end
    if nargin > 7
        loop = varargin{4};
    else
        loop = 0;
    end
    flag = 1;    
    variable = 3;
    while abs(Ttotal) > 1e-2 && single_phase_ant < 7 && out_iter < 50  
            out_iter = out_iter + 1;
        if (single_phase_ant < 4) || ~isempty(loop)  
            soluc_trivial = 0;
            for i = 1:num_sust
                soluc_trivial = soluc_trivial + (abs(Liquido.conc(i) - y(i)));
            end
            if ((nargin < 7) || (isempty(varargin{3}))) || ((soluc_trivial < 5e-5 && out_iter > 1 ) && t_semilla == 0)
                t_semilla = zeros(1,num_sust);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Estimado inicial de T (usa antoine si puede o base de datos DIPPR
                for i=1:num_sust
                    try
                        tsat = Comp(i).tsat{1};
                        try 
                            antoineT = Comp(i).tant{1};
                            pmin = Comp(i).tant{2};
                            pmax = Comp(i).tant{3};
                            tmin = Comp(i).pant{2};
                            tmax = Comp(i).pant{3};
                            t_ant = antoineT(P);
                            if t_ant > 0
                                if t_ant < tmax && t_ant > tmin
                                    antT = t_ant;
                                elseif t_ant < tmin
                                    if P < pmin
                                        antT = t_ant;
                                    else 
                                        antT = tmin;
                                    end
                                else 
                                    if P > pmax && P < comp(i).pcri
                                        antT = t_ant;
                                    else
                                        antT = tmax;
                                    end
                                end                        
                            end
                            options = optimset('Displax', 'none');
                            t_semilla(i) = fzero(@(T) tsat(T, P), antT, options);
                        catch me
                            antT = Comp(i).tcri - 0.5;
                            options = optimset('Display', 'none');
                            t_semilla(i) = fzero(@(T) tsat(T, P), antT, options);
                        end
                    catch
                        t_semilla(i) = 0;
                    end
                end
                
                t_semilla = sum(t_semilla(t_semilla~=0).*y(t_semilla~=0)); %Si algún compuesto no se calculó el Tsat lo ignora
                t_semilla(t_semilla == 0) = 380;
                soluc_trivial = 0;
                for i = 1:num_sust
                    soluc_trivial = soluc_trivial + (abs(Liquido.conc(i) - y(i)));
                end
                if Ttotal == 10 
                    T = t_semilla;
                else
                    t_semilla = T;
                end
            else
                t_semilla = varargin{3};
                if Ttotal == 10
                    T = t_semilla;
                else
                    t_semilla = T;
                end
            end
            dT = 0.015 .* T;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            x = mezcla.conc;
            if nargin < 5 && out_iter < 2 
                tol = 5e-8;
            elseif out_iter < 2
                tol = varargin{1};
                if isempty(tol)
                    tol = 5e-8;
                end
            end
            if nargin < 6 && out_iter < 2
                max_iter = 100;
            elseif out_iter < 2
                max_iter = varargin{2};
                if isempty(max_iter)
                    max_iter = 100;
                end
            end
            soluc_trivial = 0;
            for i = 1:num_sust
                soluc_trivial = soluc_trivial + (abs(Liquido.conc(i) - x(i)));
            end

            if ((nargin < 4) || (isempty(sup_x)) ) || (soluc_trivial < 5e-5)
                antP = zeros(1, num_sust);
                K = zeros(1, num_sust);
                sup_x = zeros(1, num_sust);
                fiG = P; %fugacidades de vapor
                for i = 1 : num_sust
                    try
                        antP1 = Comp(i).psat{1};
                        antP(i) = antP1(T);
                    catch
                        antP(i) = 1000;
                    end
                end
                fiL = antP;  %Obtiene la fugacidad de líquido
                for i = 1:num_sust %Coeficientes de distribución de líquido y gas (volatilidades) y composiciones calculadas de líquido
                    K(i) = fiL(i)/fiG;
                    sup_x(i) = Gas.conc(i)./K(i);
                end
                sup_x = sup_x ./ (sum(sup_x));
                sup_x2 = sup_x;
                Liquido = Mezcla(Comp, sup_x, kij);
            else
                soluc_trivial = 0;
                for i = 1:num_sust
                    soluc_trivial = soluc_trivial + (abs(sup_x(i) - x(i)));
                end
                if soluc_trivial > 5e-5
                    sup_x2 = sup_x;
                else
                    antP = zeros(1, num_sust);
                    K = zeros(1, num_sust);
                    sup_x2 = zeros(1, num_sust);
                    fiG = P; %fugacidades de vapor
                    for i = 1 : num_sust
                        antP1 = Comp(i).psat{1};
                        antP(i) = antP1(T);
                    end
                    fiL = antP;  %Obtiene la fugacidad de líquido
                    for i = 1:num_sust %Coeficientes de distribución de líquido y gas (volatilidades) y composiciones calculadas de líquido
                        K(i) = fiL(i)/fiG;
                        sup_x2(i) = Gas.conc(i)./K(i);
                    end
                    sup_x2 = sup_x2 ./ (sum(sup_x2));
                    Liquido = Mezcla(Comp, sup_x, kij);
                end
            end

            valI = 1;
            error_ = 1;
            iter = 0;
            x_calc_ant = sup_x;    
            K = zeros(1, num_sust);
            x_calc = zeros(1, num_sust); 
            while (((valI > tol) || (abs(error_)>tol))) && (iter < max_iter)  %Itera mientras no cumpla la tolerancia o
                                                    % no sea superado el número de iteraciones
                iter = iter + 1;
                fiG = self.fug(T,P,Gas,'vap'); %fugacidades de vapor
                fiL = self.fug(T,P,Liquido,'liq');  %Obtiene la fugacidad de líquido
                for i = 1:num_sust %Coeficientes de distribución de líquido y gas (volatilidades) y composiciones calculadas de líquido
                    K(i) = fiL(i)/fiG(i);
                    x_calc(i) = Gas.conc(i)/K(i);
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
                            Tsup = (T - T_ant)./(error_-error_ant).*(-error_ant)+ T_ant;
                            dT = dT ./ 4;
                        else
                            Tsup = T + sign(error_)*dT;
                        end
                    else
                        Tsup = T + sign(error_)*dT;
                    end
                    error_ant = error_;
                    T_ant = T;
                    T = Tsup;    
                end
                x_calc_ant = x_calc;
                Liquido.conc = x_calc./sum(x_calc);

            end
            if (all((abs(K - 1)) < 5e-6)) && single_phase == 0
                single_phase = 1;
            elseif single_phase == 1
                if (all((abs(K - 1)) > 5e-6))
                    single_phase = 0;
                    single_phase_ant = 0;
                else
                    single_phase_ant = single_phase_ant + 1;
                end
            end
            Ttotal = T - t_semilla;
            if single_phase == 1
            t_semilla2 = zeros(1, num_sust);
            for i=1:num_sust
                    tsat = Comp(i).tsat{1};
                    try 
                        antoineT = Comp(i).tant{1};
                        pmin = Comp(i).tant{2};
                        pmax = Comp(i).tant{3};
                        tmin = Comp(i).pant{2};
                        tmax = Comp(i).pant{3};
                        t_ant = antoineT(P);
                        if t_ant > 0
                            if t_ant < tmax && t_ant > tmin
                                antT = t_ant;
                            elseif t_ant < tmin
                                if P < pmin
                                    antT = t_ant;
                                else 
                                    antT = tmin;
                                end
                            else 
                                if P > pmax && P < comp(i).pcri
                                    antT = t_ant;
                                else
                                    antT = tmax;
                                end
                            end                        
                        end
                        options = optimset('Displax', 'none');
                        t_semilla2(i) = fzero(@(T) tsat(T, P), antT, options);
                    catch me
                        antT = Comp(i).tcri - 0.5;
                        options = optimset('Display', 'none');
                        t_semilla2(i) = fzero(@(T) tsat(T, P), antT, options);
                    end
            end

                sup_x2 = sup_x2(t_semilla2 ~= 0);
                sup_x2 = sup_x2./sum(sup_x2);
                t_semilla2 = sum(t_semilla2(t_semilla2 ~= 0).*sup_x2); %Si algún compuesto no se calculó el Tsat lo ignora
            end
            if single_phase == 1;
                if out_iter > 1 && sign(T - t_semilla2) == signo
                    variable = variable ./ 3;
                    T = T + out_iter.*variable;
                    Ttotal = T - t_semilla;
                    signo = sign(T - t_semilla2);
                elseif out_iter == 1
                    T = T + out_iter.*variable;
                    Ttotal = T - t_semilla;
                    signo = sign(T - t_semilla2);
                end
            else
                signo = sign(T - t_semilla);
            end
        else
            P1 = P - 0.1.*P; 
            if loop <=5
                [t_semilla, sup_x, K, flag] = self.DewT(P1, mezcla, [], tol, max_iter, [], loop+1);   
            else
                flag = 10;
                break
            end

            if flag == 1
                emergency_t = t_semilla;
            end
                Ttotal = t_semilla - T;
                if flag ~= 1
                    P1 = P1 - 0.1.*P1;
                    if loop <=5
                        [t_semilla, sup_x, K, flag] = self.DewT(P1, mezcla, [], tol, max_iter, [], loop+1);
                    else
                        flag = 10;
                        break
                    end
                            
                    if flag == 1
                        P1 = P1/0.9;
                        if loop <=5
                            [t_semilla, sup_x, K, flag] = self.DewT(P1, mezcla, sup_x, tol, max_iter, t_semilla, loop+1);
                        else
                            flag = 10;
                            break
                        end
                            
                        if flag == 1
                            if loop <=5
                                [t_semilla, sup_x, K, flag] = self.DewT(P, mezcla, sup_x, tol, max_iter, t_semilla, loop+1);
                            else
                                flag = 10;
                                break
                            end
                            
                            if flag == 1
                                single_phase = 0;
                            end
                            Ttotal = t_semilla - T;
                            T = t_semilla;
                        end
                    else
                        single_phase_ant = single_phase_ant + 1;
                        Ttotal = T - t_semilla;
                        T = t_semilla;
                    end
                else
                    if loop <=5
                        [t_semilla, sup_x, K, flag] = self.DewT(P, mezcla, sup_x, tol, max_iter, t_semilla, loop+1);
                    else
                        flag = 10;
                        break
                    end
                    if flag ~= 1
                        single_phase_ant = single_phase_ant + 1;
                    else
                        single_phase = 0;
                        T = t_semilla;
                    end
                    Ttotal = t_semilla - T;
                     
                end
                Liquido.conc = sup_x./sum(sup_x);
                soluc_trivial = 0;
                for i = 1:num_sust
                    soluc_trivial = soluc_trivial + (abs(Liquido.conc(i) - x(i)));
                end
                if (soluc_trivial) < 1e-5
                    single_phase = 1;
                end

        end
    end
    if flag ~= 10
        if single_phase == 1 && isempty(loop)
            try
                T = emergency_t;
            catch
            end
            warning('MATLAB:IdealEdE. Se alcanzó la región de solo una fase para el componente.');
            flag = 2;
        else        
            flag = 1;
        end
        if (out_iter == max_iter) && isempty(loop) %Si el máximo número de iteraciones se alcanzó, muestra un warning, 
            s = sprintf('%f',valI);
            warning('MATLAB:IdealEdE', ['Error de convergencia en DewPoint. El valor final de la función handle fue: ' s '.']);
            flag = 3;
        end
    else
        s = sprintf('%f',valI);
        warning('MATLAB:IdealEdE', ['Error de convergencia en DewPoint. El valor final de la función handle fue: ' s '.']);
    end
    val = sum(x_calc) - 1; %Calcula la suma Xi - 1 (debe hacerse cero) y regresa el resultado
    x = Liquido.conc;
end