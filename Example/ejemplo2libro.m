%Defino un conjunto de sustancias C2H6, C3H8, iC4H10, nC4H10, iC5H10,
%nC5H10
sust = [Sustancia('ethane'), Sustancia('Propane'), Sustancia('2-methylpropane'), Sustancia('Butane'), Sustancia('2-methylbutane') , Sustancia('Pentane')];

%la concentracion de la alimentación 1 es:
conc1 = [0.5246, 0.3106, 0.0601, 0.0688, 0.018, 0.0101];

%la concentracion de la alimentacion 2 es:

conc2 = [1.2e-4, 0.5646 0.14501, 0.15361, 0.076867, 0.05976];

%Según PRO/II los siguientes parámetros de interacción binaria son para
%Soave (1972)
kij = [-2.2e-3, -1e-2, 6.7e-3, 0, 5.6e-3,-1e-2, 0, 7.8e-3, 0.0233, 1.1e-3, 0, 0, 0, 0.0204, 0];


mezcla1 = Mezcla(sust, conc1, kij);

mezcla2 = Mezcla(sust, conc2, kij);

%Se utilizará la ecuacion de estado de Soave-Redlich-Kwong

EdE = SRKEdE();

%Se usarán las reglas de mezclado de Van der Waals

MEdE = RMVdW(EdE);

% Teniendo todo se puede calcular la primera corriente, a 1400 kPa y en su
% punto de burbuja y un flujo molar de 154.3 kg-mol/h, de la forma

P = 1400;
vargin1 = 'P'; %P significa el argumento introducido es una presión en kPa
x = 0;
vargin2 = 'x'; %x significa utilizar fracción vaporizada de argumento
flujo1 = 154.3;
varginflujo = 'm'; %m significa utilizar un flujo en base molar


corriente1 = Corriente(mezcla1, P, vargin1, x, vargin2, MEdE, flujo, varginflujo, MEdE, 'S-01');


%De la misma manera, una segunda corriente a la misma presión de 1400 kPa y
%su punto de burbuja de 64.02 kg-mol/hr se define de la forma:

flujo2 = 64.02;
corriente2 = Corriente(mezcla2, P, vargin1, x, vargin2, MEdE, flujo2, varginflujo, MEdE, 'S-01');


% La primera alimentación entra en la etapa 6 y la segunda en la 11, así
% que el "cell" de alimentación, alim, se proveerá de la forma:

alim = {6, corriente1, 11, corriente2};

% El número de etapas teóricas de la columna es de 21, incluyendo rehervidor 
% y condensador

etapas = 21;

%La columna es prácticamente isobárica a 1380 kPa salvo por el condensador
%a menor presión (1320 kPa) y el rehervidor a la mayor presión (1440 kPa)

presion_por_plato = ones(1,21).*1380;
presion_por_plato(1) = 1320;
presion_por_plato(end) = 1440;


%Se desea condensador total, es decir, una salida de destilado líquido en
% el tope, siendo esta de 79.8 kg-mol/hr. De ser el destilado vapor, no es
% necesario definir este cell

salidatope = {1,1,79.8};


%Un par de corrientes placeholders en las cuales se definirían los
%productos (aún no implementado)

dest = Corriente.empty(0,1);

fondo = Corriente.empty(0,1);


% El calor intercambiado por plato es despreciable, la columna es adiabática

calor = 0;

%El reflujo (L/D) es de 2

reflujo = 2;

% La otra condición es dada por el balance de masa global, con la cual al
% entrar 218.32 kg-mol alimentados y salir en el tope 79.8 kg-mol, en el
% fondo van a salir 138.52 kg-mol y se le suministra 'B' para indicar el
% total de fondo

flujofondo = 138.52 kg-mol
condflujo = 'B';

%Se le permitirán realizar 200 iteraciones

iteraciones = 200;

%Se tomará el tiempo a la corrida

tiempo = tic;

torreBP = BubblePoint(etapas, alim, dest, fondo, salidatope, presion_por_plato, calor, reflujo, flujofondo, condflujo, iteraciones);

%se muestra el tiempo total de corrida

toc(tiempo)

%Para la torre de Newton Raphson se usará un amortiguamiento bajo de 0.4

damping = 0.4;


NaphtaliSandholm(etapas, alim, dest, fondo, salidatope, presion_por_plato, calor, iteraciones, [], [], damping, [], , 2, 'R', torrey.bsubi, 'bi', torrey.perfil_t, 'tj')
