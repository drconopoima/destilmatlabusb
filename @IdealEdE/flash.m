function [x, y, beta, K, flag, valI , iter] = flash(self, T, P, mezcla, sup_x, sup_y, sup_beta, varargin)
%[x, y, beta, K, flag, val] = flash(T, P, mezcla, sup_x, sup_y, sup_beta, tol, max_iter)
% Encuentra BETA para condiciones T, P y mezcla composicion Z
%   BETA es la fracci�n vaporizada de la mezcla de l�quido vapor.
%   Si se err� en la temperatura, a la cual no pudiese haber mezcla, los
%   valores de beta del resultado no estar�n comprendidos entre 0 y 1.
% Luis Jes�s D�az Manzo
%
%Par�metros
%
%Entrada:
% P: Presi�n real de la mezcla
% T: Temperatura real de la mezcla
%mezcla: Clase Mezcla.m que contiene las propiedades de la mezcla como
%composici�n, x clases Sustancia.m que definen las propiedades de sustancia
%pura
%sup_x, sup_y: suposiciones iniciales de las composiciones l�quida y vapor
%sup_beta: suposici�n inicial de la fracci�n vaporizada. 
%tol: Tolerancia del c�lculo sucesivo de 2 composiciones de l�quido y de vapor
%(default 1e-6)
%max_iter: n�mero m�ximo de iteraciones del algoritmo (default = 100)
%
%Salida:
%x: Composici�n molar de l�quido de la mezcla
%y: Composici�n molar de vapor de la mezcla
%beta: Fracci�n vaporizada (calidad) de la mezcla l�quida y vapor
%K: volatilidades de los componentes de la mezcla, como fugL/fugV o y/x
%(coeficientes de distribuci�n de fases l�quido y gas)
%flag: flag del m�todo de resoluci�n num�rico:
%     1  Encontr� un cero X.
%     2  C�lculos convergieron a la regi�n de una fase
%    -1  Algorithm terminated by output function.
%    -3  NaN or Inf function value encountered during search for an interval
%         containing a sign change.
%    -4  Complex function value encountered during search for an interval 
%         containing a sign change.
%    -5  FZERO may have converged to a singular point.
%    -6  FZERO can not detect a change in sign of the function.

num_sust = mezcla.num_sust;
if (nargin < 5) || (isempty(sup_x))
    sup_x = zeros(1, num_sust);
    sup_x(end) = 1;
end
if (nargin < 6) || (isempty(sup_y))
    sup_y = zeros(1, num_sust);
    sup_y(1) = 1;
end
if nargin > 7
	if ~isempty(varargin{1})
		tol = varargin{1};
	end
else
    tol = 1e-6;
end
if nargin > 8
	if ~isempty(varargin{2})
		max_iter = varargin{2};
    end
else
    max_iter = 100;
end

options = optimset('Display', 'none', 'TolX', tol, 'TolFun', tol, 'MaxIter', max_iter, 'MaxFunEvals', max_iter);

y = sup_y;
x = sup_x;
z = mezcla.conc;
beta = sup_beta;
comp = mezcla.comp;
kij = mezcla.kij;
Gas = Mezcla(comp, y, kij);
Liquido = Mezcla(comp, x, kij);
valI = 1; %
iter = 0;

while ((valI > tol) && (iter < max_iter))
	iter = iter + 1;
	x_ant = x;
	y_ant = y;
	fiG = self.fug(T, P, Gas, 'vap'); %fugacidades de vapor
	fiL = self.fug(T, P, Liquido, 'liq'); %Obtiene la fugacidad de l�quido
	K = zeros(1, num_sust);
	for i = 1:num_sust %Coeficientes de distribuci�n de l�quido y gas (volatilidades)
		K(i) = fiL(i)/fiG(i);
	end
	[beta, ~, flag] = fsolve(@(b) RachfordRice(b, z, K), beta, options);
	if beta < 0
		beta = 1e-8;
	elseif beta > 1
		beta = 1 - 1e-8;
	end
	
	%Recalcula las composiciones de l�quido y de vapor
	x = z./(1+beta.*(K-1));
	y = K.*Liquido.conc;
	
	Gas.conc = y/sum(y);
	Liquido.conc = x/sum(x);
	
	%Calcula la diferencia entre iteraciones sucesivas
	valI = 0;
	for i =1: num_sust
		valI = valI + abs(x(i) - x_ant(i)) + abs(y(i) - y_ant(i));
	end
end

single_phase = 0;
if ((1 - 1.1e-8) < beta )|| (1.1e-8 > beta)
    single_phase = 1;
end
if single_phase == 1
    flag = 2;
end

if (iter == max_iter) %Si se alcanz� el m�ximo n�mero de iteraciones muestra warning
    s = sprintf('%f',valI);
    warning('MATLAB:IdealEdEflash',['Error de convergencia en flash. El valor final de las funciones de convergencia: ' s '']);
end

x=x/sum(x);
y=y/sum(y);
end