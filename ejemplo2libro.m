%Defino un conjunto de sustancias C2H6, C3H8, nC4H10, nC5H10, nC6H12
sust = [Sustancia('ethane'), Sustancia('Propane'), Sustancia('Butane'), Sustancia('Pentane'), Sustancia('Hexane')];

% la concentracion de la alimentación 1 es:
conc1 = [0.5342, 0.3163, 0.0612, 0.0701, 0.0183];

%la concentracion de la alimentacion 2 es:
conc2 = [0.2236, 0.5194, 0.1119, 0.1230, 0.0221];

%kij ideales, también se hubiera podido usar los de Peng y Robinson arrojados por PRO/II
kij = 0;

mezcla1 = Mezcla(sust, conc1, kij);
mezcla2 = Mezcla(sust, conc2, kij);

%Se utiliza la ecuacion de estado de Peng y Robinson
 
EdE = PREdE();
%Se usan las reglas de mezclado de Van der Waals
 
MEdE = RMVdW(EdE);

% Teniendo todo se puede calcular la primera corriente, a 1320 kPa y en su
% punto de burbuja y un flujo molar de 154.3 kg-mol/h, de la forma:
 
P = 1320;
vargin1 = 'P';  %P significa el argumento introducido es una presión en kPa
x = 0;
vargin2 = 'x'; %x significa utilizar fracción vaporizada de argumento
flujo1 = 154.3;
varginflujo = 'm'; %m significa utilizar un flujo en base molar

corriente1 = Corriente(mezcla1, P, vargin1, x, vargin2, flujo1, varginflujo, MEdE, 'S-01');

%El flujo de la corriente 2 es de 64.02 kg-mol

flujo2 = 64.02;

corriente2 = Corriente(mezcla2, P, vargin1, x, vargin2, flujo2, varginflujo, MEdE, 'S-02');

alim = {4, corriente1, 6, corriente2};

% El número de etapas teóricas es 9, incluyendo rehervidor y el condensador
etapas = 9;

%La columna es prácticamente isobárica a 1320 kPa salvo por el condensador
%a menor presión (1260 kPa) y el rehervidor a la mayor presión (1380 kPa)

presion_por_plato = ones(1,9).*P;
presion_por_plato(1) = 1260;
presion_por_plato(end) = 1380;

%Se desea condensador total, es decir, una salida de destilado líquido en
% el tope, siendo esta de 79.8 kg-mol/hr. De ser el destilado vapor, no es
% necesario definir 
 
salidatope = {1,1,79.8};

%Un par de corrientes placeholders en las cuales se definirían los
%productos (aún no implementado)
 
dest = Corriente.empty(0,1);
 
fondo = Corriente.empty(0,1);

% El calor intercambiado por plato es despreciable, la columna es adiabática
 
calor = 0;
%El reflujo (L/D) es de 2
 
reflujo = 2.5;

% La otra condición es dada por el balance de masa global, con la cual al
% entrar 218.32 kg-mol alimentados y salir en el tope 79.8 kg-mol, en el
% fondo van a salir 138.52 kg-mol y se le suministra 'B' para indicar el
% total de fondo
 
flujofondo = 138.52;
condflujo = 'B';

%Se le restringirá que el máximo a realizar son 200 iteraciones
 
iteraciones = 200;

%Ahora sí teniendo todos los datos requeridos se procede a definir la torre

tiempo = tic;
torreBP = BubblePoint(etapas, alim, dest, fondo, salidatope, presion_por_plato, calor, reflujo, flujofondo, condflujo, iteraciones);

%Se realiza el balance de materia

torreBP.balanmasa();

%Inicio del método Bubble Point de Wang y Henke
 
torreBP.wanghenke();
toc(tiempo)
%Para la torre de Newton Raphson se define un amortiguamiento bajo de 0.45
 
damping = 0.40;

tiempo = tic;
torreNS = NaphtaliSandholm(etapas, alim, dest, fondo, salidatope, presion_por_plato, calor, iteraciones, [], [], [], [], damping,  reflujo, 'R', flujofondo, condflujo,  torreBP.Fi - torreBP.dsubi, 'bi', torreBP.perfil_t, 'Tj'); 

% Se realiza una distribución ideal de los compuestos mediante un balance de masa

torreNS.balanmasa();

%Puesto que no se posee generador de estimados iniciales, se procede a utilizar como %estimados de Newton Raphson los valores finales de la convergencia de BubblePoint

torreNS.perfil_t = torreBP.perfil_t; %Estimados de temperatura
torreNS.perfil_v = torreBP.perfil_v;  %Estimados de flujo molar de vapor por plato
torreNS.perfil_l = torreBP.perfil_l; %Estimados de flujo molar de líquido por plato
torreNS.perfil_li = torreBP.perfil_li; %Estimados de flujo molar líquido por componente por plato
torreNS.perfil_vi = torreBP.perfil_vi; %Estimados de flujo molar vapor por componente por plato
torreNS.perfil_k = torreBP.perfil_k; %Estimados de coeficiente de reparto por componente por plato
torreNS.xsubj = torreBP.xsubj; %Estimados de fracción molar de fase líquida por componente por plato
torreNS.ysubj = torreBP.ysubj;  %Estimados de fracción molar de fase vapor por componente por plato

%Debido a que el flujo de vapor destilado debe restringirse a cero, para tener un estimado %consistente se debe restringir a 0 el destilado vapor
torreNS.dflujo = 0;
torreNS.perfil_v(1) = 0;

%El estimado de flujo de residuo líquido por la última etapa es el mismo que con el método %de BubblePoint
torreNS.bflujo = torreBP.bflujo;

%Finalmente teniendo un set de estimados iniciales que es consistente y que debería %generarse acoplando un algoritmo de Youseff para la estimación de perfiles por plato y por %componente mediante un cálculo corto

torreNS.newtonraphson2(); %El 2 hace referencia a la especificación de reflujo y flujo desilado
toc(tiempo)