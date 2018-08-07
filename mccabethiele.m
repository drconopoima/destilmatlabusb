function varargout = mccabethiele(varargin)
% MCCABETHIELE Código MATLAB para la ejecución de mccabethiele.fig
%
%      Función de control de la interfaz gráfica para el método gráfico binario
%      de McCabe y Thiele.      
%
%
%      MCCABETHIELE, by itself, creates a new MCCABETHIELE or raises the existing
%      singleton*.
%
%      H = MCCABETHIELE returns the handle to a new MCCABETHIELE or the handle to
%      the existing singleton*.
%
%      MCCABETHIELE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MCCABETHIELE.M with the given input arguments.
%
%      MCCABETHIELE('Property','Value',...) creates a new MCCABETHIELE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mccabethiele_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mccabethiele_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mccabethiele

% Last Modified by GUIDE v2.5 01-May-2015 21:38:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mccabethiele_OpeningFcn, ...
                   'gui_OutputFcn',  @mccabethiele_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before mccabethiele is made visible.
function mccabethiele_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mccabethiele (see VARARGIN)

% Choose default command line output for mccabethiele
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mccabethiele wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mccabethiele_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Ejecuta cambios en la selección del compuestoLK.
function compuestoLK_Callback(hObject, eventdata, handles)
% Crea una clase sustancia que contiene las propiedades de la sust. deseada 
%almacena la clase sustancia en la variable sustLK111


% hObject    handle to compuestoLK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns compuestoLK contents as cell array
%        contents{get(hObject,'Value')} returns selected item from compuestoLK
global sustLK111;  global ident_sustLK;
global Tc111; global Pc111; global W111;
    comp1=get(handles.compuestoLK,'Value');
    if comp1 ~= 2 && comp1 ~= 1
        load('indices1.mat'); %Carga la base de datos de sustancias
        ident_sustLK = database1(comp1-2);
        sustLK111 = Sustancia(ident_sustLK); 
        %%% Pasan a mostrarse los datos críticos de la sustancia
        %%% seleccionada
        set(handles.text22, 'ForegroundColor', [0,0,0]);
        set(handles.text23, 'ForegroundColor', [0,0,0]);
        set(handles.text24, 'ForegroundColor', [0,0,0]);
        set(handles.Pc1,'String',sustLK111.pcri)
        set(handles.Tc1,'String',sustLK111.tcri)
        set(handles.W1,'String',sustLK111.w_acent)
        %%% Se inicializan en una variable global
        Pc111  = sustLK111.pcri;
        Tc111=sustLK111.tcri;
        W111=sustLK111.w_acent;
        set(handles.livianotag,'ForegroundColor',[0 0 0])
    elseif comp1 == 1
        set(handles.livianotag,'ForegroundColor',[1 0 0])
        set(handles.text22, 'ForegroundColor', [0,0,0]);
        set(handles.text23, 'ForegroundColor', [0,0,0]);
        set(handles.text24, 'ForegroundColor', [0,0,0]);
        clearvars sustLK111; clearvars Pc111; clearvars Tc111; clearvars W111;
        set(handles.Pc1, 'String', ''); set(handles.Pc1, 'String', ''); set(handles.W1, 'String', '');
        
    else
        clearvars ident_sustLK; clearvars Tc111; clearvars Pc111; clearvars W111;
        set(handles.modelotermo, 'ForegroundColor', [1,0,0]);
        set(handles.livianotag,'ForegroundColor',[1 0 0])
        clearvars sustLK111; clearvars Pc111; clearvars Tc111; clearvars W111;
        set(handles.Pc1, 'String', ''); set(handles.Tc1, 'String', ''); set(handles.W1, 'String', '');
        set(handles.text22, 'ForegroundColor', [1, 0, 0]); set(handles.text23, 'ForegroundColor', [1, 0, 0]); set(handles.text24,'ForegroundColor', [1,0,0]);
    end
if all(get(handles.livianotag, 'ForegroundColor') == [0,0,0]) && all(get(handles.pesadotag, 'ForegroundColor') == [0,0,0])
    set(handles.modelotermo, 'ForegroundColor', [0,0,0]);
end

% --- Executes during object creation, after setting all properties.
function compuestoLK_CreateFcn(hObject, eventdata, handles)
% hObject    handle to compuestoLK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Ejecuta cambios en la selección del compuestoHK.
function compuestoHK_Callback(hObject, eventdata, handles)
% hObject    handle to compuestoHK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns compuestoHK contents as cell array
%        contents{get(hObject,'Value')} returns selected item from compuestoHK
global sustHK111; global ident_sustHK;
global Tc112; global Pc112; global W112; 
    comp2=get(handles.compuestoHK,'Value');
    if comp2 ~= 2 && comp2 ~= 1
        load('indices1.mat'); %Carga la base de datos de sustancias
        ident_sustHK = database1(comp2-2);
        sustHK111 = Sustancia(ident_sustHK); 
        %%% Pasan a mostrarse los datos críticos de la sustancia
        %%% seleccionada
        set(handles.text25, 'ForegroundColor', [0,0,0]);
        set(handles.text26, 'ForegroundColor', [0,0,0]);
        set(handles.text27, 'ForegroundColor', [0,0,0]);
        set(handles.Pc2,'String',sustHK111.pcri)
        set(handles.Tc2,'String',sustHK111.tcri)
        set(handles.W2,'String',sustHK111.w_acent)
        %%% Se inicializan en una variable global
        Pc112  = sustHK111.pcri;
        Tc112=sustHK111.tcri;
        W112=sustHK111.w_acent;
        set(handles.pesadotag,'ForegroundColor',[0 0 0])
    elseif comp2 == 1
        set(handles.pesadotag,'ForegroundColor',[1 0 0])
        set(handles.text25, 'ForegroundColor', [0,0,0]);
        set(handles.text26, 'ForegroundColor', [0,0,0]);
        set(handles.text27, 'ForegroundColor', [0,0,0]);
        clearvars sustHK111; clearvars Pc112; clearvars Tc112; clearvars W112;
    else
        clearvars ident_sustHK; clearvars Tc112; clearvars Pc112; clearvars W112;
        set(handles.modelotermo, 'ForegroundColor', [1,0,0]);
        set(handles.pesadotag,'ForegroundColor',[1 0 0])
        clearvars sustHK111; clearvars Pc112; clearvars Tc112; clearvars W112;
        set(handles.Pc2, 'String', ''); set(handles.Tc2, 'String', ''); set(handles.W2, 'String', '');
        set(handles.text25, 'ForegroundColor', [1, 0, 0]); set(handles.text26, 'ForegroundColor', [1, 0, 0]); set(handles.text27,'ForegroundColor', [1,0,0]);

    end
if all(get(handles.livianotag, 'ForegroundColor') == [0,0,0]) && all(get(handles.pesadotag, 'ForegroundColor') == [0,0,0])
    set(handles.modelotermo, 'ForegroundColor', [0,0,0]);
end

% --- Executes during object creation, after setting all properties.
function compuestoHK_CreateFcn(hObject, eventdata, handles)
% hObject    handle to compuestoHK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function P_Callback(hObject, eventdata, handles)
% Define la presión del sistema a la que se calcula el equilibrio Liq-Vap

% hObject    handle to P (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global vP111;
vP111=get(handles.P,'string');
vP111=str2double(vP111);
cnt=vP111;
if isnan(cnt)==1
    errordlg('Debe introducir un valor para la Presión del Sistema',' Error en inserción de datos ');
    set(handles.presiontag,'ForegroundColor',[1 0 0])
else
    set(handles.presiontag,'ForegroundColor',[0 0 0])
end

% --- Executes during object creation, after setting all properties.
function P_CreateFcn(hObject, eventdata, handles)
% hObject    handle to P (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kij_Callback(hObject, eventdata, handles)
%Permite introducir un kij de interacción binaria en la mezcla

% hObject    handle to kij (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global kij111
    kij111 = get(handles.kij, 'String');
    kij111 = str2double(kij111);
    cnt = kij111;
    if isnan(cnt)==1
        errordlg('Debe introducir un valor para la Presión del Sistema',' Error en inserción de datos ');
        set(handles.kijtag,'ForegroundColor',[0 0 1])
        set(handles.kij, 'String', eventdata.Source.Value);
    else
        set(handles.kijtag, 'ForegroundColor', [0, 0, 0]);
    end


% --- Executes during object creation, after setting all properties.
function kij_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kij (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function N_Callback(hObject, eventdata, handles)
%N_Callback determina el número de elementos en la curva de equilibrio


% hObject    handle to N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global N111;

    N111 = get(handles.N, 'String');
    N111 = str2double(N111);
% Hints: get(hObject,'String') returns contents of N as text
%        str2double(get(hObject,'String')) returns contents of N as a double



% --- Executes during object creation, after setting all properties.
function N_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in caleq.
function caleq_Callback(hObject, eventdata, handles)
%Ejecución del cálculo de la curva de equilibrio y llenado de la tabla
%uitable1

% hObject    handle to caleq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sustHK111; global sustLK111; global N111; global vP111;
global MEdE111; global kij111; global cmpbinario111
    correr=get(handles.modelotermo,'Value');
    if correr == 3 || correr == 2 || correr == 4 || correr == 5
        vP111=get(handles.P,'string');
        vP111=str2double(vP111);
        try
            if  ~isempty(N111)
                x = linspace(0, 1, 2*N111-1);
            else
                N111 = 41;
                x = linspace(0, 1, 2*N111-1);
            end
        catch
            N111 = 41;
            x = linspace(0, 1, 2*N111-1);
        end
        T = x;
        y = zeros(2, length(x));
        y(:,end) = [1;1];
            try
                T(1) = MEdE111.EdE.isofugT(vP111, sustHK111, 1e-7);
                T(end) = MEdE111.EdE.isofugT(vP111, sustLK111, 1e-7);
            catch ME
                tsat = sustHK111.tsat;
                tsat = tsat{1};
                options = optimset('Display', 'none');
                T(1) = fzero(@(T) tsat(T, vP111), (sustHK111.tcri - 1), options);
                tsat = sustLK111.tsat;
                tsat = tsat{1};
                T(end) = fzero(@(T) tsat(T, vP111), (sustLK111.tcri - 1), options);
            end
            %consigo la curva de equilibrio para 39 puntos
            cmpbinario111= [sustLK111, sustHK111];
            for i = 2:2:length(x)-1
                try
                    if ~isempty(kij111)
                        mezcla = Mezcla(cmpbinario111, x(i), kij111);
                    else
                        kij111 = 0;
                        mezcla = Mezcla(cmpbinario111, x(i), kij111);
                    end
                catch
                    kij111 = 0;
                    mezcla = Mezcla(cmpbinario111, x(i), kij111);
                end
                T(i) = (T(end) - T(1))*x(i) + T(1);
                try
                    [T(i), y(:,i)] = MEdE111.BubbleT(vP111, mezcla, [], [], [], T(i));
                catch
                    errordlg('Probablemente se encuentran LK y HK invertidos en volatilidad')
                end
            end
            y = y(1,:);
            % interpolo el resto de los puntos mediante Lagrange cuadrático
            for i = 3:2:length(x)-1
                y(i) = y(i-2)*(((x(i) - x(i-1))*(x(i)-x(i+1)))/((x(i-2) - x(i-1))*(x(i-2) - x(i+1)))) + y(i-1)*(((x(i) - x(i-2))*(x(i)-x(i+1)))/((x(i-1) - x(i-2))*(x(i-1) - x(i+1)))) + y(i+1)*(((x(i) - x(i-2))*(x(i)-x(i-1)))/((x(i+1) - x(i-2))*(x(i+1) - x(i-1))));
                T(i) = T(i-2)*(((x(i) - x(i-1))*(x(i)-x(i+1)))/((x(i-2) - x(i-1))*(x(i-2) - x(i+1)))) + T(i-1)*(((x(i) - x(i-2))*(x(i)-x(i+1)))/((x(i-1) - x(i-2))*(x(i-1) - x(i+1)))) + T(i+1)*(((x(i) - x(i-2))*(x(i)-x(i-1)))/((x(i+1) - x(i-2))*(x(i+1) - x(i-1))));
            end
        set(handles.uitable1,'Data',[x',y', T'])
        global X111; global T111; global Y111; global XE111; global YE111;
        X111 = x; Y111 = y; T111 = T; 
        XE111 = X111; YE111 = Y111;
    end
set(handles.modelotermo, 'ForegroundColor', [0,0,0]);


% --- Ejecuta cambios en la selección de modelotermo.
function modelotermo_Callback(hObject, eventdata, handles)
% Decide el modelo termodinamico a utilizar o define uno el usuario
%
%

% hObject    handle to modelotermo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns modelotermo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from modelotermo
global MEdE111;
correr=get(handles.modelotermo,'Value');
if correr==2 || correr == 3 || correr == 4 || correr == 5
    set(handles.modelotermo, 'Foregroundcolor', [1,0,0]);
    set(handles.volatilpanel,'Visible','off')
    set(handles.EdEpanel,'Visible','on')
    
    %Se indica al usuario introducir datos (color rojo)
    set(handles.livianotag, 'ForegroundColor', [1,0,0]);
    set(handles.compuestoLK, 'Value', 1);
    set(handles.compuestoHK, 'Value', 1);
    set(handles.presiontag, 'ForegroundColor', [1,0,0]);
    set(handles.kijtag, 'ForegroundColor', [0, 0, 1]);
    set(handles.uitable1,'Position',[2.0 0.6923076923076923 52.4 7.53846153846154])   
    set(handles.uitable1,'Visible','on');
    set(handles.crit_data, 'Visible', 'on');
    %%% Se inicializa el modelo termodinámico seleccionado
    if correr == 2
        MEdE111 = IdealEdE();
    elseif correr == 3
        MEdE111 = SRKEdE();
        MEdE111 = RMVdW(MEdE111);
    elseif correr == 4
        MEdE111 = PREdE();
        MEdE111 = RMVdW(MEdE111);
    elseif correr == 5
        MEdE111 = PRGEdE();
        MEdE111 = RMVdW(MEdE111);
    end
global N111;
N111 = 41; 
elseif correr == 6
    set(handles.modelotermo, 'Foregroundcolor', [1,0,0]);
    set(handles.EdEpanel,'Visible','off')
    set(handles.crit_data, 'Visible', 'off')
    set(handles.uitable1,'Position',[2.0 0.6923076923076923 52.4 24.53846153846154])   
    set(handles.uitable1,'Visible','on');
    set(handles.volatilpanel, 'Visible', 'on');
    set(handles.crit_data, 'Visible', 'off');
    %En rojo la data que necesita esta termodinamica
    set(handles.volatilpanel, 'ForegroundColor', [1,0,0]);
    set(handles.modelotermo, 'ForegroundColor', [1,0,0]);
    %En negro todos los titulos de data que no necesita esta termodinamica
    set(handles.livianotag, 'ForegroundColor', [0,0,0]);
    set(handles.pesadotag, 'ForegroundColor', [0,0,0]);
    set(handles.presiontag, 'ForegroundColor', [0,0,0]);
    clearvars MEdE111
end
% --- Executes during object creation, after setting all properties.
function modelotermo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to modelotermo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Pc1_Callback(hObject, eventdata, handles)
global Pc111; global sustLK111; global Tc111; global W111; global ident_sustLK;
Pc111 = str2double(get(handles.Pc1, 'String'));
set(handles.text22, 'ForegroundColor', [0,0,0]);
    try         
        if get(handles.compuestoLK, 'Value') ~= 2
            sustLK111 = Sustancia(ident_sustLK, [], Tc111);
        elseif ~isempty(Pc111) && ~isempty(Tc111) && ~isempty(W111)
            sustLK111 = Sustancia([], [],  Tc111, Pc111, W111);
        end    
        if ~isempty(get(handles.Tc1, 'String')) && ~isempty(get(handles.Pc1, 'String')) && ~isempty(get(handles.W1, 'String'))
            set(handles.livianotag, 'ForegroundColor', [0,0,0]);
        end
    catch
    end
if all(get(handles.text22, 'ForegroundColor') == [0,0,0]) && all(get(handles.text23, 'ForegroundColor') == [0,0,0]) && all(get(handles.text24, 'ForegroundColor') == [0,0,0]) && all(get(handles.text25, 'ForegroundColor') == [0,0,0]) && all(get(handles.text26, 'ForegroundColor') == [0,0,0]) && all(get(handles.text27, 'ForegroundColor') == [0,0,0]) && all(get(handles.livianotag, 'ForegroundColor') == [0,0,0]) && all(get(handles.pesadotag, 'ForegroundColor') == [0,0,0])
    set(handles.modelotermo, 'ForegroundColor', [0,0,0]);
end
function Pc1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Tc1_Callback(hObject, eventdata, handles)
global Tc111; global sustLK111; global Pc111; global W111; global ident_sustLK;
Tc111 = str2double(get(handles.Tc1, 'String'));
set(handles.text23, 'ForegroundColor', [0,0,0]);
    try 
        if get(handles.compuestoLK, 'Value') ~= 2
            sustLK111 = Sustancia(ident_sustLK, [], [], Pc111);
        elseif ~isempty(Pc111) && ~isempty(Tc111) && ~isempty(W111)
            sustLK111 = Sustancia([], [],  Tc111, Pc111, W111);
        end
        if ~isempty(get(handles.Tc1, 'String')) && ~isempty(get(handles.Pc1, 'String')) && ~isempty(get(handles.W1, 'String'))
            set(handles.livianotag, 'ForegroundColor', [0,0,0]);
        end
    catch
    end
if all(get(handles.text22, 'ForegroundColor') == [0,0,0]) && all(get(handles.text23, 'ForegroundColor') == [0,0,0]) && all(get(handles.text24, 'ForegroundColor') == [0,0,0]) && all(get(handles.text25, 'ForegroundColor') == [0,0,0]) && all(get(handles.text26, 'ForegroundColor') == [0,0,0]) && all(get(handles.text27, 'ForegroundColor') == [0,0,0]) && all(get(handles.livianotag, 'ForegroundColor') == [0,0,0]) && all(get(handles.pesadotag, 'ForegroundColor') == [0,0,0])
    set(handles.modelotermo, 'ForegroundColor', [0,0,0]);
end

function Tc1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function W1_Callback(hObject, eventdata, handles)
global W111; global sustLK111; global Tc111; global Pc111; global ident_sustLK;
W111 = str2double(get(handles.W1, 'String'));
set(handles.text24, 'ForegroundColor', [0,0,0]);
    try 
        if get(handles.compuestoLK, 'Value') ~= 2
            sustLK111 = Sustancia(ident_sustLK, [], [], [], W111);
        elseif ~isempty(Pc111) && ~isempty(Tc111) && ~isempty(W111)
            sustLK111 = Sustancia([], [],  Tc111, Pc111, W111);
        end
        if ~isempty(get(handles.Tc1, 'String')) && ~isempty(get(handles.Pc1, 'String')) && ~isempty(get(handles.W1, 'String'))
            set(handles.livianotag, 'ForegroundColor', [0,0,0]);
        end
    catch
    end
if all(get(handles.text22, 'ForegroundColor') == [0,0,0]) && all(get(handles.text23, 'ForegroundColor') == [0,0,0]) && all(get(handles.text24, 'ForegroundColor') == [0,0,0]) && all(get(handles.text25, 'ForegroundColor') == [0,0,0]) && all(get(handles.text26, 'ForegroundColor') == [0,0,0]) && all(get(handles.text27, 'ForegroundColor') == [0,0,0]) && all(get(handles.livianotag, 'ForegroundColor') == [0,0,0]) && all(get(handles.pesadotag, 'ForegroundColor') == [0,0,0])
    set(handles.modelotermo, 'ForegroundColor', [0,0,0]);
end

function W1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Pc2_Callback(hObject, eventdata, handles)
global Pc112; global sustHK112; global Tc112; global W112; global ident_sustHK;
Pc112 = str2double(get(handles.Pc2, 'String'));
set(handles.text25, 'ForegroundColor', [0,0,0]);
    try 
        if get(handles.compuestoHK, 'Value') ~= 2
            sustHK112 = Sustancia(ident_sustHK, [], Tc112);
        elseif ~isempty(Pc112) && ~isempty(Tc112) && ~isempty(W112)
            sustHK112 = Sustancia([], [],  Tc112, Pc112, W112);
        end
        if ~isempty(get(handles.Tc2, 'String')) && ~isempty(get(handles.Pc2, 'String')) && ~isempty(get(handles.W2, 'String'))
            set(handles.pesadotag, 'ForegroundColor', [0,0,0]);
        end
    catch
    end
if all(get(handles.text22, 'ForegroundColor') == [0,0,0]) && all(get(handles.text23, 'ForegroundColor') == [0,0,0]) && all(get(handles.text24, 'ForegroundColor') == [0,0,0]) && all(get(handles.text25, 'ForegroundColor') == [0,0,0]) && all(get(handles.text26, 'ForegroundColor') == [0,0,0]) && all(get(handles.text27, 'ForegroundColor') == [0,0,0]) && all(get(handles.livianotag, 'ForegroundColor') == [0,0,0]) && all(get(handles.pesadotag, 'ForegroundColor') == [0,0,0])
    set(handles.modelotermo, 'ForegroundColor', [0,0,0]);
end
    
function Pc2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Tc2_Callback(hObject, eventdata, handles)
global Tc112; global sustHK112; global Pc112; global W112; global ident_sustHK;
Tc112 = str2double(get(handles.Tc2, 'String'));
set(handles.text26, 'ForegroundColor', [0,0,0]);
    try     
        if get(handles.compuestoHK, 'Value') ~= 2
            sustHK112 = Sustancia(ident_sustHK, [], [], Pc112);
        elseif ~isempty(Pc112) && ~isempty(Tc112) && ~isempty(W112)
            sustHK112 = Sustancia([], [],  Tc112, Pc112, W112);
        end
        if ~isempty(get(handles.Tc2, 'String')) && ~isempty(get(handles.Pc2, 'String')) && ~isempty(get(handles.W2, 'String'))
            set(handles.pesadotag, 'ForegroundColor', [0,0,0]);
        end
    catch
    end
if all(get(handles.text22, 'ForegroundColor') == [0,0,0]) && all(get(handles.text23, 'ForegroundColor') == [0,0,0]) && all(get(handles.text24, 'ForegroundColor') == [0,0,0]) && all(get(handles.text25, 'ForegroundColor') == [0,0,0]) && all(get(handles.text26, 'ForegroundColor') == [0,0,0]) && all(get(handles.text27, 'ForegroundColor') == [0,0,0]) && all(get(handles.livianotag, 'ForegroundColor') == [0,0,0]) && all(get(handles.pesadotag, 'ForegroundColor') == [0,0,0])
    set(handles.modelotermo, 'ForegroundColor', [0,0,0]);
end

function Tc2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function W2_Callback(hObject, eventdata, handles)
global Pc112; global sustHK112; global Tc112; global W112; global ident_sustHK;
W112 = str2double(get(handles.W2, 'String'));
set(handles.text27, 'ForegroundColor', [0,0,0]);
    try 
        if get(handles.compuestoHK, 'Value') ~= 2
            sustHK112 = Sustancia(ident_sustHK, [], [], [], W112);
        elseif ~isempty(Pc112) && ~isempty(Tc112) && ~isempty(W112)
            sustHK112 = Sustancia([], [],  Tc112, Pc112, W112);
        end
        if ~isempty(get(handles.Tc2, 'String')) && ~isempty(get(handles.Pc2, 'String')) && ~isempty(get(handles.W2, 'String'))
            set(handles.pesadotag, 'ForegroundColor', [0,0,0]);
        end
    catch
    end
if all(get(handles.text22, 'ForegroundColor') == [0,0,0]) && all(get(handles.text23, 'ForegroundColor') == [0,0,0]) && all(get(handles.text24, 'ForegroundColor') == [0,0,0]) && all(get(handles.text25, 'ForegroundColor') == [0,0,0]) && all(get(handles.text26, 'ForegroundColor') == [0,0,0]) && all(get(handles.text27, 'ForegroundColor') == [0,0,0]) && all(get(handles.livianotag, 'ForegroundColor') == [0,0,0]) && all(get(handles.pesadotag, 'ForegroundColor') == [0,0,0])
    set(handles.modelotermo, 'ForegroundColor', [0,0,0]);
end
    
function W2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function volatil_Callback(hObject, eventdata, handles)
% hObject    handle to volatil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of volatil as text
%        str2double(get(hObject,'String')) returns contents of volatil as a double

global alfa111
alfa111 = str2double(get(handles.volatil, 'String'));
set(handles.volatilpanel,'Foregroundcolor', [0,0,0]);
set(handles.modelotermo, 'Foregroundcolor', [0,0,0]);

% --- Executes during object creation, after setting all properties.
function volatil_CreateFcn(hObject, eventdata, handles)
% hObject    handle to volatil (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function volatilN_Callback(hObject, eventdata, handles)
% hObject    handle to volatilN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of volatilN as text
%        str2double(get(hObject,'String')) returns contents of volatilN as a double
global N111
N111 = str2double(get(handles.volatilN, 'String'));

% --- Executes during object creation, after setting all properties.
function volatilN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to volatilN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global N111; global alfa111; global X111; global Y111; 
try
    if  ~isempty(N111)
        x = linspace(0, 1, 2*N111-1);
    else
        N111 = 41;
        x = linspace(0, 1, 2*N111-1);
    end
catch
    N111 = 41;
    x = linspace(0, 1, 2*N111-1);
end
try
    y = alfa111.*x./(1 + (alfa111 - 1).*x);
catch
    errordlg('Debe especificar una volatilidad relativa')
end
X111 = x; Y111 = y;
set(handles.uitable1,'Data',[x',y', zeros(length(x),1)]);


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5



function E_Callback(hObject, eventdata, handles)
% hObject    handle to E (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Effic; global X111; global Y111; global XE111; global YE111;
Effic = str2double(get(hObject, 'String'));
if Effic > 100 || Effic < 0
    Effic = [];
    warning('La eficiencia debe estar comprendida entre 0 y 100');
    set(hObject, 'String', '')
else
    XE111 = X111; 
    for iter = 1:length(X111)
        YE111(iter) = (Y111(iter) - X111(iter))*Effic/100  + X111(iter);
    end
end


% Hints: get(hObject,'String') returns contents of E as text
%        str2double(get(hObject,'String')) returns contents of E as a double


% --- Executes during object creation, after setting all properties.
function E_CreateFcn(hObject, eventdata, handles)
% hObject    handle to E (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Re_Callback(hObject, eventdata, handles)
% hObject    handle to Re (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Re111;

    Re111=str2double(get(handles.Re,'String'));
    
if isnan(Re111)==1
    errordlg('Debe introducir un valor para el reflujo (L/D)',' Error en inserción de datos ');
    set(handles.uipanel7,'ForegroundColor',[1 0 0])
else
    set(handles.uipanel7,'ForegroundColor',[0 0 0])
end
% Hints: get(hObject,'String') returns contents of Re as text
%        str2double(get(hObject,'String')) returns contents of Re as a double


% --- Executes during object creation, after setting all properties.
function Re_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Re (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in xRem.
function xRem_Callback(hObject, eventdata, handles)
% hObject    handle to xRem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of xRem



function Xd_Callback(hObject, eventdata, handles)
% hObject    handle to Xd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Xd as text
%        str2double(get(hObject,'String')) returns contents of Xd as a double
global xDrD;
xDrD=get(handles.Xd,'String');
xDrD=str2double(xDrD);
if isnan(xDrD)==1
    errordlg('Debe introducir un valor entre 0 y 1 para la composición',' Error en inserción de datos ');
    set(handles.composic_dest,'ForegroundColor',[1 0 0])
elseif xDrD>1
    errordlg('Debe introducir un valor entre 0 y 1 para la composición',' Error en inserción de datos ');
    set(handles.composic_dest,'ForegroundColor',[1 0 0])
elseif xDrD<0
    errordlg('Debe introducir un valor entre 0 y 1 para la composición',' Error en inserción de datos ');
    set(handles.composic_dest,'ForegroundColor',[1 0 0])
else
    set(handles.composic_dest,'ForegroundColor',[0 0 0])
end

% --- Executes during object creation, after setting all properties.
function Xd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Xw_Callback(hObject, eventdata, handles)
% hObject    handle to Xw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Xw as text
%        str2double(get(hObject,'String')) returns contents of Xw as a double
global xWrW
xWrW=get(handles.Xw,'String');
xWrW=str2double(xWrW);
if isnan(xWrW)==1
    errordlg('Debe introducir un valor entre 0 y 1 para la composición',' Error en inserción de datos ');
    set(handles.composic_res,'ForegroundColor',[1 0 0])
elseif xWrW>1
    errordlg('Debe introducir un valor entre 0 y 1 para la composición',' Error en inserción de datos ');
    set(handles.composic_res,'ForegroundColor',[1 0 0])
elseif xWrW<0
    errordlg('Debe introducir un valor entre 0 y 1 para la composición',' Error en inserción de datos ');
    set(handles.composic_res,'ForegroundColor',[1 0 0])
else
    set(handles.composic_res,'ForegroundColor',[0 0 0])
end


% --- Executes during object creation, after setting all properties.
function Xw_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in composic_dest.
function composic_dest_Callback(hObject, eventdata, handles)
% hObject    handle to composic_dest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns composic_dest contents as cell array
%        contents{get(hObject,'Value')} returns selected item from composic_dest
global condicion_mccabe1

condicion_mccabe1 = get(hObject, 'Value');


% --- Executes during object creation, after setting all properties.
function composic_dest_CreateFcn(hObject, eventdata, handles)
% hObject    handle to composic_dest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in composic_res.
function composic_res_Callback(hObject, eventdata, handles)
% hObject    handle to composic_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns composic_res contents as cell array
%        contents{get(hObject,'Value')} returns selected item from composic_res
global condicion_mccabe2

condicion_mccabe2 = get(hObject, 'Value');

% --- Executes during object creation, after setting all properties.
function composic_res_CreateFcn(hObject, eventdata, handles)
% hObject    handle to composic_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function F_Callback(hObject, eventdata, handles)
% hObject    handle to F (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of F as text
%        str2double(get(hObject,'String')) returns contents of F as a double
global F111
F111=get(handles.F,'String');
F111=str2double(F111);
if isnan(F111)==1
    errordlg('Debe introducir un valor para el flujo de alimentación',' Error en inserción de datos ');
    set(handles.text2,'ForegroundColor',[1 0 0])
else
    set(handles.text2,'ForegroundColor',[0 0 0])
end

% --- Executes during object creation, after setting all properties.
function F_CreateFcn(hObject, eventdata, handles)
% hObject    handle to F (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Zf_Callback(hObject, eventdata, handles)
% hObject    handle to Zf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Zf as text
%        str2double(get(hObject,'String')) returns contents of Zf as a double
global Zf111
Zf111=get(handles.Zf,'String');
Zf111=str2double(Zf111);
if isnan(Zf111)==1
    errordlg('Debe introducir un valor entre 0 y 1 para la composición',' Error en inserción de datos ');
    set(handles.text1,'ForegroundColor',[1 0 0])
elseif Zf111>1
    errordlg('Debe introducir un valor entre 0 y 1 para la composición',' Error en inserción de datos ');
    set(handles.text1,'ForegroundColor',[1 0 0])
elseif Zf111<0
    errordlg('Debe introducir un valor entre 0 y 1 para la composición',' Error en inserción de datos ');
    set(handles.text1,'ForegroundColor',[1 0 0])
else
    set(handles.text1,'ForegroundColor',[0 0 0])
end

% --- Executes during object creation, after setting all properties.
function Zf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Zf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function qf_Callback(hObject, eventdata, handles)
% hObject    handle to qf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of qf as text
%        str2double(get(hObject,'String')) returns contents of qf as a double
global qf111; global Tf111; global vapf111; global condicion2
Tf111 = []; vapf111 = []; qf111 = [];

if condicion2 == 2
    Tf111 = get(handles.qf,'String');
    Tf111 = str2double(Tf111);
    set(handles.segunda_condicion, 'ForegroundColor', [0, 0, 0])
elseif condicion2 == 3
    vapf111 = get(handles.qf,'String');
    str2double(vapf111);
    set(handles.segunda_condicion, 'ForegroundColor', [0, 0, 0])
elseif condicion2 == 4
    qf111=get(handles.qf,'String');
    qf111=str2double(qf111);
    set(handles.segunda_condicion, 'ForegroundColor', [0, 0, 0])
end

% --- Executes during object creation, after setting all properties.
function qf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to qf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in segunda_condicion.
function segunda_condicion_Callback(hObject, eventdata, handles)
% hObject    handle to segunda_condicion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns segunda_condicion contents as cell array
%        contents{get(hObject,'Value')} returns selected item from segunda_condicion
global condicion2
condicion2 = get(hObject, 'Value');
if condicion2 == 1
    set(handles.segunda_condicion, 'ForegroundColor', [1, 0, 0])
end

% --- Executes during object creation, after setting all properties.
function segunda_condicion_CreateFcn(hObject, eventdata, handles)
% hObject    handle to segunda_condicion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)

global MccabeObject; global X111; global Y111; global T111; global XE111; global YE111;
global Corrienteppal; global cmpbinario111; global Zf111; global kij111;
global MEdE111; global vP111; global Effic; global vapvivo; global Re111; global salidas
salidas = {};
if all(get(handles.livianotag,'ForegroundColor')==[1 0 0])
    errordlg('Los datos requeridos no han sido introducidos completamente, verifique items resaltados en color rojo',' Error en inserción de datos ')
elseif all(get(handles.pesadotag,'ForegroundColor')==[1 0 0])
    errordlg('Los datos requeridos no han sido introducidos completamente, verifique items resaltados en color rojo',' Error en inserción de datos ')
elseif all(get(handles.presiontag,'ForegroundColor')==[1 0 0])
    errordlg('Los datos requeridos no han sido introducidos completamente, verifique items resaltados en color rojo',' Error en inserción de datos ')
elseif all(get(handles.text2,'ForegroundColor')==[1 0 0])
    errordlg('Los datos requeridos no han sido introducidos completamente, verifique items resaltados en color rojo',' Error en inserción de datos ')
elseif all(get(handles.text1,'ForegroundColor')==[1 0 0])
    errordlg('Los datos requeridos no han sido introducidos completamente, verifique items resaltados en color rojo',' Error en inserción de datos ')
elseif all(get(handles.segunda_condicion,'ForegroundColor')==[1 0 0])
    errordlg('Los datos requeridos no han sido introducidos completamente, verifique items resaltados en color rojo',' Error en inserción de datos ')
elseif all(get(handles.composic_dest,'ForegroundColor')==[1 0 0]) || all(get(handles.composic_res,'ForegroundColor')==[1 0 0]) || all(get(handles.uipanel7,'ForegroundColor')==[1 0 0])
    errordlg('Los datos requeridos no han sido introducidos completamente, verifique items resaltados en color rojo',' Error en inserción de datos ')

elseif get(handles.modelotermo, 'Value') ~=6 && get(handles.modelotermo, 'Value') ~=1 
    
if get(handles.radiobutton6,'Value') == 1
    set(handles.radiobutton5, 'Value', 0);
    Effic = 100;
end
sustancia = cmpbinario111;
conc = [Zf111, 1 - Zf111];
if ~isempty(kij111)
    mezcla = Mezcla(sustancia, conc, kij111);
else
    kij111 = 0;
    mezcla = Mezcla(cmpbinario111, conc, kij111);
end
if get(handles.segunda_condicion, 'Value') == 2
    T_o_P_o_beta = str2double(get(handles.qf, 'String'));
    topob1 = 'T';
elseif get(handles.segunda_condicion, 'Value') == 3
    T_o_P_o_beta = str2double(get(handles.qf, 'String'));
    topob1 = 'x';
elseif get(handles.segunda_condicion, 'Value') == 4
    T_o_P_o_beta = str2double(get(handles.qf, 'String'));
    T_o_P_o_beta = abs(T_o_P_o_beta - 1);
    topob1 = 'x';
end
flujo = str2double(get(handles.F, 'String'));
m_o_w_o_v = 'm';
Corrienteppal = Corriente(mezcla, T_o_P_o_beta, topob1, vP111, 'P', flujo, m_o_w_o_v , MEdE111, 'C-000');
entradas = [];
try 
    entradas(end + 1) = mezcla.conc(1);
    entradas(end + 1) = Corrienteppal(end).q;
    entradas(end + 1) = flujo;
catch
end

if strcmpi('Alimentacion',handles.uitable3.Data(1,5))
    try
        flujo = str2double(handles.uitable3.Data(1,1));
        conc = [str2double( handles.uitable3.Data(1,2) ), 1 - str2double( handles.uitable3.Data(1,2) )];
        if ~isempty(kij111)
            mezcla = Mezcla(sustancia, conc, kij111);
        else
            kij111 = 0;
            mezcla = Mezcla(cmpbinario111, conc, kij111);
        end
        try 
            T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(1,3) );
            topob1 = 'x';
        Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, vP111, 'P', flujo, m_o_w_o_v , MEdE111, 'C-001')];
        catch ME
            T_o_P_o_beta = str2double( handles.uitable3.Data(1,4) );
            topob1 = 'T';
        Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, vP111, 'P', flujo, m_o_w_o_v , MEdE111, 'C-001')];
        end

        try 
            entradas(end + 1) = conc(1);
            entradas(end + 1) = Corrienteppal(end).q;
            entradas(end + 1) = flujo;
        catch
        end
    catch
    end
elseif strcmpi('Salida',handles.uitable3.Data(1,5))    
    try
        flujo = str2double(handles.uitable3.Data(1,1));
        conc = [str2double( handles.uitable3.Data(1,2) ), 1 - str2double( handles.uitable3.Data(1,2) )];
        try 
            if abs((str2double( handles.uitable3.Data(1,3) )) - 1)<1e-5 || abs((str2double( handles.uitable3.Data(1,3) )))<1e-5
                T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(1,3) );
            else
                error('No es adecuado el valor')
            end
        catch
            errordlg('En el caso de salidas solo puede utilizarse condicion térmica y debe ser 0 o 1');
        end
        salidas = {salidas{:}, conc, flujo, T_o_P_o_beta};
    catch
    end
end

if strcmpi('Alimentacion', handles.uitable3.Data(2,5) )
    try
        flujo = str2double( handles.uitable3.Data(2,1) );
        conc = [str2double( handles.uitable3.Data(2,2) ), 1 - str2double( handles.uitable3.Data(2,2) )];
        if ~isempty(kij111)
            mezcla = Mezcla(sustancia, conc, kij111);
        else
            kij111 = 0;
            mezcla = Mezcla(cmpbinario111, conc, kij111);
        end
        try 
            T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(2,3) );
            topob1 = 'x';
        Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, vP111, 'P', flujo, m_o_w_o_v , MEdE111, 'C-002')];
        catch ME
            T_o_P_o_beta = str2double( handles.uitable3.Data(2,4) );
            topob1 = 'T';
        Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, vP111, 'P', flujo, m_o_w_o_v , MEdE111, 'C-002')];
        end
        
        try 
            entradas(end + 1) = conc(1);
            entradas(end + 1) = Corrienteppal(end).q;
            entradas(end + 1) = flujo;
        catch
        end
    catch
    end
elseif strcmpi('Salida',handles.uitable3.Data(2,5))    
    try
        flujo = str2double(handles.uitable3.Data(2,1));
        conc = [str2double( handles.uitable3.Data(2,2) ), 1 - str2double( handles.uitable3.Data(2,2) )];
        try 
            if abs((str2double( handles.uitable3.Data(2,3) )) - 1)<1e-5 || abs((str2double( handles.uitable3.Data(2,3) )))<1e-5
                T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(2,3) );
            else
                error('No es adecuado el valor')
            end
        catch
            errordlg('En el caso de salidas solo puede utilizarse condicion térmica y debe ser 0 o 1');
        end
        salidas = {salidas{:}, conc, flujo, T_o_P_o_beta};
    catch
    end
end    
if strcmpi('Alimentacion', handles.uitable3.Data(3,5) )
    try
        flujo = str2double( handles.uitable3.Data(3,1) );
        conc = [str2double( handles.uitable3.Data(3,2) ), 1 - str2double( handles.uitable3.Data(3,2) )];
        if ~isempty(kij111)
            mezcla = Mezcla(sustancia, conc, kij111);
        else
            kij111 = 0;
            mezcla = Mezcla(cmpbinario111, conc, kij111);
        end
        try 
            T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(3,3) );
            topob1 = 'x';
        Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, vP111, 'P', flujo, m_o_w_o_v , MEdE111, 'C-002')];
        catch ME
            T_o_P_o_beta = str2double( handles.uitable3.Data(3,4) );
            topob1 = 'T';
        Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, vP111, 'P', flujo, m_o_w_o_v , MEdE111, 'C-002')];
        end
        try 
            entradas(end + 1) = conc(1);
            entradas(end + 1) = Corrienteppal(end).q;
            entradas(end + 1) = flujo;
        catch
        end
    catch
    end
elseif strcmpi('Salida',handles.uitable3.Data(3,5))    
    try
        flujo = str2double(handles.uitable3.Data(3,1));
        conc = [str2double( handles.uitable3.Data(3,2) ), 1 - str2double( handles.uitable3.Data(3,2) )];
        try 
            if abs((str2double( handles.uitable3.Data(3,3) )) - 1)<1e-5 || abs((str2double( handles.uitable3.Data(3,3) )))<1e-5
                T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(3,3) );
            else
                error('No es adecuado el valor')
            end
        catch
            errordlg('En el caso de salidas solo puede utilizarse condicion térmica y debe ser 0 o 1');
        end
        salidas = {salidas{:}, conc, flujo, T_o_P_o_beta};
    catch
    end
end

if strcmpi('Alimentacion', handles.uitable3.Data(4,5) )
    try
        flujo = str2double( handles.uitable3.Data(4,1) );
        conc = [str2double( handles.uitable3.Data(4,2) ), 1 - str2double( handles.uitable3.Data(4,2) )];
        if ~isempty(kij111)
            mezcla = Mezcla(sustancia, conc, kij111);
        else
            kij111 = 0;
            mezcla = Mezcla(cmpbinario111, conc, kij111);
        end
        try 
            T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(4,3) );
            topob1 = 'x';
        Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, vP111, 'P', flujo, m_o_w_o_v , MEdE111, 'C-003')];
        catch ME
            T_o_P_o_beta = str2double( handles.uitable3.Data(4,4) );
            topob1 = 'T';
        Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, vP111, 'P', flujo, m_o_w_o_v , MEdE111, 'C-003')];
        end
        
        try 
            entradas(end + 1) = conc(1);
            entradas(end + 1) = Corrienteppal(end).q;
            entradas(end + 1) = flujo;
        catch
        end
    catch
    end
elseif strcmpi('Salida',handles.uitable3.Data(4,5))    
    try
        flujo = str2double(handles.uitable3.Data(4,1));
        conc = [str2double( handles.uitable3.Data(4,2) ), 1 - str2double( handles.uitable3.Data(4,2) )];
        try 
            if abs((str2double( handles.uitable3.Data(4,3) )) - 1)<1e-5 || abs((str2double( handles.uitable3.Data(4,3) )))<1e-5
                T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(4,3) );
            else
                error('No es adecuado el valor')
            end
        catch
            errordlg('En el caso de salidas solo puede utilizarse condicion térmica y debe ser 0 o 1');
        end
        salidas = {salidas{:}, conc, flujo, T_o_P_o_beta};
    catch
    end
end


B = 0;
D = 0;
xD = str2double(get(handles.Xd, 'String'));
xB = str2double(get(handles.Xw, 'String'));
vapvivo = get(handles.radiobutton10, 'Value');

if get(handles.composic_dest, 'Value') == 1 

    if get(handles.composic_res, 'Value') == 1
        MccabeObject = Hengstebeck(Corrienteppal, 1, [], 2, [], salidas, [], [], [], [], [], vapvivo, get(handles.radiobutton4, 'Value'));
        MccabeObject.xDi =  str2double(get(handles.Xd, 'String'));
        MccabeObject.xBi =  str2double(get(handles.Xw, 'String'));
        for i = 1:length(Corrienteppal)
            B = B + (Corrienteppal(i).molF*(xD - Corrienteppal(i).conc(1)));
            D = D + (Corrienteppal(i).molF*(Corrienteppal(i).conc(1) - xB));
        end
        D = D/(xD - xB);
        B = B/(xD - xB);
        MccabeObject.Di = [D*xD, D*(1 - xD)];
        MccabeObject.Bi = [B*xB, B*(1 - xB)];
    else
        MccabeObject = Hengstebeck(Corrienteppal, 1, [], 2, str2double(get(handles.Xw, 'String')),  salidas, [], [], [], [], [], vapvivo, get(handles.radiobutton4, 'Value'));
        for i = 1:length(Corrienteppal)
            B = B + (Corrienteppal(i).molF);
            D = D + (Corrienteppal(i).molF*(Corrienteppal(i).conc(1) - xB));
        end
        F = B;
        B = xB*D;
        D = D/(xD - xB);
        B = B/(xD - xB);
        MccabeObject.Di = [D*xD, D*(1 - xD)];
        MccabeObject.Bi = [B*xB, B*(1 - xB)];
        MccabeObject.xDi =  str2double(get(handles.Xd, 'String'));
    end
else
    if get(handles.composic_res, 'Value') == 1
        MccabeObject = Hengstebeck(Corrienteppal, 1, str2double(get(handles.Xd, 'String')), 2, [],  salidas, [], [], [], [], [], vapvivo, get(handles.radiobutton4, 'Value'));
        MccabeObject.xBi =  str2double(get(handles.Xw, 'String'));
    else
        MccabeObject = Hengstebeck(Corrienteppal, 1, str2double(get(handles.Xd, 'String')), 2, str2double(get(handles.Xw, 'String')), salidas, Re111, [], [], [], Effic, vapvivo, get(handles.radiobutton4, 'Value'));
    end

end
MccabeObject.entradas = entradas;
MccabeObject.E_x = X111;
MccabeObject.E_y = Y111;
MccabeObject.Tbin = T111;
if ~isempty(Effic)
    MccabeObject.Murphree = Effic/100;
else
    MccabeObject.Murphree = 1;
end
MccabeObject.Mur_x = XE111;
MccabeObject.Mur_y = YE111;
MccabeObject.alter_define(MccabeObject.alimentaciones, 1, MccabeObject.reco_lk, 2, MccabeObject.reco_hk, MccabeObject.E_x);
if get(handles.xRem, 'Value')
    MccabeObject.distributemin();
    MccabeObject.reflux_remin(Re111);
    Re111 = MccabeObject.reflujo;
end
MccabeObject.distribute(Re111);

set(handles.nmin, 'String', MccabeObject.platosaefficmin)
set(handles.eteo, 'String', MccabeObject.platosmin)
if get(handles.radiobutton6, 'Value') == 1
    set(handles.ereal, 'String', ceil(MccabeObject.platosmin./((str2double(get(handles.E, 'String'))/100))))
else
    set(handles.ereal, 'String', ceil(MccabeObject.etapas))
end
set(handles.eopt, 'String', MccabeObject.etapasemin(1)+1);
set(handles.dest, 'String', sum(MccabeObject.Di));
set(handles.resi, 'String', sum(MccabeObject.Bi));
set(handles.re_tope, 'String', MccabeObject.Rmin);
set(handles.retope, 'String', MccabeObject.Rmin.*sum(MccabeObject.Di)./(sum(MccabeObject.Di).*(1 + MccabeObject.Rmin)));
set(handles.re_fondo, 'String', MccabeObject.vap_min);
for i = 1:4
    if strcmpi('Alimentacion',handles.uitable3.Data(i,5))
        handles.uitable2.Data(i,1)= {MccabeObject.etapasemin(i+1) + sum(MccabeObject.etapasemin(1:i))+1};
    elseif strcmpi('Salida', handles.uitable3.Data(i,5))
        handles.uitable2.Data(i,1)= {MccabeObject.etapasemin(i+1) + sum(MccabeObject.etapasemin(1:i))+1};
    end
end
secciones = {'I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX'};
handles.uitable4.Data(1,1) = {secciones{1}};
handles.uitable4.Data(1,2) = {MccabeObject.rectas(1,1)};
handles.uitable4.Data(1,3) = {MccabeObject.rectas(2,1)};

handles.uitable4.Data(2,1) = {'Última'};
handles.uitable4.Data(2,2) = {MccabeObject.rectas(1,2)};
handles.uitable4.Data(2,3) = {MccabeObject.rectas(2,2)};

for i = 2:5
    try
        if length(Mccabethiele.rectas) > 2
            handles.uitable4.Data(i+1,2) = {MccabeObject.rectas(i,1)};
            handles.uitable4.Data(i+1,3) = {MccabeObject.rectas(i,1)};
            handles.uitable4.Data(i+1,1) = {secciones{i}};
        end
    catch
        break
    end
end

else
    if get(handles.radiobutton6,'Value') == 1
        set(handles.radiobutton5, 'Value', 0);
        Effic = 100;
    end
    cmpbinario111 = [Sustancia('Ethane'), Sustancia('Propane')];
    sustancia = cmpbinario111;
    conc = [Zf111, 1 - Zf111];
    if ~isempty(kij111)
        mezcla = Mezcla(sustancia, conc, kij111);
    else
        kij111 = 0;
        mezcla = Mezcla(cmpbinario111, conc, kij111);
    end
    if get(handles.segunda_condicion, 'Value') == 2
        errordlg('No se puede definir con temperatura')
    elseif get(handles.segunda_condicion, 'Value') == 3
        T_o_P_o_beta = str2double(get(handles.qf, 'String'));
        topob1 = 'x';
    elseif get(handles.segunda_condicion, 'Value') == 4
        T_o_P_o_beta = str2double(get(handles.qf, 'String'));
        T_o_P_o_beta = abs(T_o_P_o_beta - 1);
        topob1 = 'x';
    end
    flujo = str2double(get(handles.F, 'String'));
    m_o_w_o_v = 'm';
    Corrienteppal = Corriente(mezcla, T_o_P_o_beta, topob1, 100, 'P', flujo, m_o_w_o_v , MEdE111, 'C-000');
    entradas = [];
    try 
        entradas(end + 1) = mezcla.conc(1);
        entradas(end + 1) = Corrienteppal(end).q;
        entradas(end + 1) = flujo;
    catch
    end
    if strcmpi('Alimentacion',handles.uitable3.Data(1,5))
        try
            flujo = str2double(handles.uitable3.Data(1,1));
            conc = [str2double( handles.uitable3.Data(1,2) ), 1 - str2double( handles.uitable3.Data(1,2) )];
            if ~isempty(kij111)
                mezcla = Mezcla(sustancia, conc, kij111);
            else
                kij111 = 0;
                mezcla = Mezcla(cmpbinario111, conc, kij111);
            end
            try 
                T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(1,3) );
                topob1 = 'x';
            Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1,100, 'P', flujo, m_o_w_o_v , MEdE111, 'C-001')];
            catch ME
                T_o_P_o_beta = str2double( handles.uitable3.Data(1,4) );
                topob1 = 'T';
            Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1,100, 'P', flujo, m_o_w_o_v , MEdE111, 'C-001')];
            end

            try 
                entradas(end + 1) = conc(1);
                entradas(end + 1) = Corrienteppal(end).q;
                entradas(end + 1) = flujo;
            catch
            end
        catch
        end
    end

    if strcmpi('Alimentacion', handles.uitable3.Data(2,5) )
        try
            flujo = str2double( handles.uitable3.Data(2,1) );
            conc = [str2double( handles.uitable3.Data(2,2) ), 1 - str2double( handles.uitable3.Data(2,2) )];
            if ~isempty(kij111)
                mezcla = Mezcla(sustancia, conc, kij111);
            else
                kij111 = 0;
                mezcla = Mezcla(cmpbinario111, conc, kij111);
            end
            try 
                T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(2,3) );
                topob1 = 'x';
            Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1,100, 'P', flujo, m_o_w_o_v , MEdE111, 'C-002')];
            catch ME
                T_o_P_o_beta = str2double( handles.uitable3.Data(2,4) );
                topob1 = 'T';
            Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, 100, 'P', flujo, m_o_w_o_v , MEdE111, 'C-002')];
            end

            try 
                entradas(end + 1) = conc(1);
                entradas(end + 1) = Corrienteppal(end).q;
                entradas(end + 1) = flujo;
            catch
            end
        catch
        end
    end    
    if strcmpi('Alimentacion', handles.uitable3.Data(3,5) )
        try
            flujo = str2double( handles.uitable3.Data(3,1) );
            conc = [str2double( handles.uitable3.Data(3,2) ), 1 - str2double( handles.uitable3.Data(3,2) )];
            if ~isempty(kij111)
                mezcla = Mezcla(sustancia, conc, kij111);
            else
                kij111 = 0;
                mezcla = Mezcla(cmpbinario111, conc, kij111);
            end
            try 
                T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(3,3) );
                topob1 = 'x';
            Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, 100, 'P', flujo, m_o_w_o_v , MEdE111, 'C-002')];
            catch ME
                T_o_P_o_beta = str2double( handles.uitable3.Data(3,4) );
                topob1 = 'T';
            Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, 100, 'P', flujo, m_o_w_o_v , MEdE111, 'C-002')];
            end
            try 
                entradas(end + 1) = conc(1);
                entradas(end + 1) = Corrienteppal(end).q;
                entradas(end + 1) = flujo;
            catch
            end
        catch
        end
    end

    if strcmpi('Alimentacion', handles.uitable3.Data(4,5) )
        try
            flujo = str2double( handles.uitable3.Data(4,1) );
            conc = [str2double( handles.uitable3.Data(4,2) ), 1 - str2double( handles.uitable3.Data(4,2) )];
            if ~isempty(kij111)
                mezcla = Mezcla(sustancia, conc, kij111);
            else
                kij111 = 0;
                mezcla = Mezcla(cmpbinario111, conc, kij111);
            end
            try 
                T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(4,3) );
                topob1 = 'x';
            Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, 100, 'P', flujo, m_o_w_o_v , MEdE111, 'C-003')];
            catch ME
                T_o_P_o_beta = str2double( handles.uitable3.Data(4,4) );
                topob1 = 'T';
            Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, 100, 'P', flujo, m_o_w_o_v , MEdE111, 'C-003')];
            end

            try 
                entradas(end + 1) = conc(1);
                entradas(end + 1) = Corrienteppal(end).q;
                entradas(end + 1) = flujo;
            catch
            end
        catch
        end
    end

    B = 0;
    D = 0;
    xD = str2double(get(handles.Xd, 'String'));
    xB = str2double(get(handles.Xw, 'String'));
    vapvivo = get(handles.radiobutton10, 'Value');

    if get(handles.composic_dest, 'Value') == 1 

        if get(handles.composic_res, 'Value') == 1
            MccabeObject = Hengstebeck(Corrienteppal, 1, [], 2, [], [], [], [], [], [], [], vapvivo, get(handles.radiobutton4, 'Value'));
            MccabeObject.xDi =  str2double(get(handles.Xd, 'String'));
            MccabeObject.xBi =  str2double(get(handles.Xw, 'String'));
            for i = 1:length(Corrienteppal)
                B = B + (Corrienteppal(i).molF*(xD - Corrienteppal(i).conc(1)));
                D = D + (Corrienteppal(i).molF*(Corrienteppal(i).conc(1) - xB));
            end
            D = D/(xD - xB);
            B = B/(xD - xB);
            MccabeObject.Di = [D*xD, D*(1 - xD)];
            MccabeObject.Bi = [B*xB, B*(1 - xB)];
        else
            MccabeObject = Hengstebeck(Corrienteppal, 1, [], 2, str2double(get(handles.Xw, 'String')),  [], [], [], [], [], [], vapvivo, get(handles.radiobutton4, 'Value'));
            for i = 1:length(Corrienteppal)
                B = B + (Corrienteppal(i).molF);
                D = D + (Corrienteppal(i).molF*(Corrienteppal(i).conc(1) - xB));
            end
            F = B;
            B = xB*D;
            D = D/(xD - xB);
            B = B/(xD - xB);
            MccabeObject.Di = [D*xD, D*(1 - xD)];
            MccabeObject.Bi = [B*xB, B*(1 - xB)];
            MccabeObject.xDi =  str2double(get(handles.Xd, 'String'));
        end
    else
        if get(handles.composic_res, 'Value') == 1
            MccabeObject = Hengstebeck(Corrienteppal, 1, str2double(get(handles.Xd, 'String')), 2, [],  [], [], [], [], [], [], vapvivo, get(handles.radiobutton4, 'Value'));
            MccabeObject.xBi =  str2double(get(handles.Xw, 'String'));
        else
            MccabeObject = Hengstebeck(Corrienteppal, 1, str2double(get(handles.Xd, 'String')), 2, str2double(get(handles.Xw, 'String')), [], Re111, [], [], [], Effic, vapvivo, get(handles.radiobutton4, 'Value'));
        end

    end
    MccabeObject.entradas = entradas;
    MccabeObject.E_x = X111;
    MccabeObject.E_y = Y111;
    MccabeObject.Tbin = T111;
    if ~isempty(Effic)
        MccabeObject.Murphree = Effic/100;
    else
        MccabeObject.Murphree = 1;
    end
    MccabeObject.Mur_x = XE111;
    MccabeObject.Mur_y = YE111;
    MccabeObject.E_x = X111;
        MccabeObject.E_y = Y111;
        MccabeObject.Mur_x = X111;
        MccabeObject.Mur_y = (Y111 - X111)* Effic/100+X111;
        MccabeObject.distributemin;
    if get(handles.xRem, 'Value')
        MccabeObject.distributemin();
        MccabeObject.reflux_remin(Re111);
        Re111 = MccabeObject.reflujo;
    end
    MccabeObject.distribute(Re111);

    set(handles.nmin, 'String', MccabeObject.platosaefficmin)
    set(handles.eteo, 'String', MccabeObject.platosmin)
    if get(handles.radiobutton6, 'Value') == 1
        set(handles.ereal, 'String', ceil(MccabeObject.platosmin./((str2double(get(handles.E, 'String'))/100))))
    else
        set(handles.ereal, 'String', ceil(MccabeObject.etapas))
    end
    set(handles.eopt, 'String', MccabeObject.etapasemin(1)+1);
    set(handles.dest, 'String', sum(MccabeObject.Di));
    set(handles.resi, 'String', sum(MccabeObject.Bi));
    set(handles.re_tope, 'String', MccabeObject.Rmin);
    set(handles.retope, 'String', MccabeObject.Rmin.*sum(MccabeObject.Di)./(sum(MccabeObject.Di).*(1 + MccabeObject.Rmin)));
    set(handles.re_fondo, 'String', MccabeObject.vap_min);
    for i = 1:4
        if strcmpi('Alimentacion',handles.uitable3.Data(i,5))
            handles.uitable2.Data(i,1)= {MccabeObject.etapasemin(i+1) + sum(MccabeObject.etapasemin(1:i))+1};
        end
    end
    secciones = {'I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX'};
    handles.uitable4.Data(1,1) = {secciones{1}};
    handles.uitable4.Data(1,2) = {MccabeObject.rectas(1,1)};
    handles.uitable4.Data(1,3) = {MccabeObject.rectas(2,1)};

    handles.uitable4.Data(2,1) = {'Última'};
    handles.uitable4.Data(2,2) = {MccabeObject.rectas(1,2)};
    handles.uitable4.Data(2,3) = {MccabeObject.rectas(2,2)};

    for i = 2:5
        try
            if length(Mccabethiele.rectas) > 2
                handles.uitable4.Data(i+1,2) = {MccabeObject.rectas(i,1)};
                handles.uitable4.Data(i+1,3) = {MccabeObject.rectas(i,1)};
                handles.uitable4.Data(i+1,1) = {secciones{i}};
            end
        catch
            break
        end
    end


end



function nmin_Callback(hObject, eventdata, handles)
% hObject    handle to nmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nmin as text
%        str2double(get(hObject,'String')) returns contents of nmin as a double


% --- Executes during object creation, after setting all properties.
function nmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function re_tope_Callback(hObject, eventdata, handles)
% hObject    handle to re_tope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of re_tope as text
%        str2double(get(hObject,'String')) returns contents of re_tope as a double


% --- Executes during object creation, after setting all properties.
function re_tope_CreateFcn(hObject, eventdata, handles)
% hObject    handle to re_tope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function re_fondo_Callback(hObject, eventdata, handles)
% hObject    handle to re_fondo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of re_fondo as text
%        str2double(get(hObject,'String')) returns contents of re_fondo as a double


% --- Executes during object creation, after setting all properties.
function re_fondo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to re_fondo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function retope_Callback(hObject, eventdata, handles)
% hObject    handle to retope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of retope as text
%        str2double(get(hObject,'String')) returns contents of retope as a double


% --- Executes during object creation, after setting all properties.
function retope_CreateFcn(hObject, eventdata, handles)
% hObject    handle to retope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton6


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
        global MccabeObject; global X111; global Y111; global T111; global XE111; global YE111;
        global Corrienteppal; global cmpbinario111; global Zf111; global kij111;
        global MEdE111; global vP111; global Effic; global vapvivo; global salidas


    if get(handles.modelotermo, 'Value') ~=6 && get(handles.modelotermo, 'Value') ~= 1
        if all(get(handles.livianotag,'ForegroundColor')==[1 0 0])
            errordlg('Los datos requeridos no han sido introducidos completamente, verifique items resaltados en color rojo',' Error en inserción de datos ')
        elseif all(get(handles.pesadotag,'ForegroundColor')==[1 0 0])
            errordlg('Los datos requeridos no han sido introducidos completamente, verifique items resaltados en color rojo',' Error en inserción de datos ')
        elseif all(get(handles.presiontag,'ForegroundColor')==[1 0 0])
            errordlg('Los datos requeridos no han sido introducidos completamente, verifique items resaltados en color rojo',' Error en inserción de datos ')
        elseif all(get(handles.text2,'ForegroundColor')==[1 0 0])
            errordlg('Los datos requeridos no han sido introducidos completamente, verifique items resaltados en color rojo',' Error en inserción de datos ')
        elseif all(get(handles.text1,'ForegroundColor')==[1 0 0])
            errordlg('Los datos requeridos no han sido introducidos completamente, verifique items resaltados en color rojo',' Error en inserción de datos ')
        elseif all(get(handles.segunda_condicion,'ForegroundColor')==[1 0 0])
            errordlg('Los datos requeridos no han sido introducidos completamente, verifique items resaltados en color rojo',' Error en inserción de datos ')
        elseif all(get(handles.composic_dest,'ForegroundColor')==[1 0 0]) || all(get(handles.composic_res,'ForegroundColor')==[1 0 0])
            errordlg('Los datos requeridos no han sido introducidos completamente, verifique items resaltados en color rojo',' Error en inserción de datos ')
        else
        if get(handles.radiobutton6,'Value') == 1
            set(handles.radiobutton5, 'Value', 0);
            Effic = 100;
        end
        sustancia = cmpbinario111;
        conc = [Zf111, 1 - Zf111];
        if ~isempty(kij111)
            mezcla = Mezcla(sustancia, conc, kij111);
        else
            kij111 = 0;
            mezcla = Mezcla(cmpbinario111, conc, kij111);
        end
        if get(handles.segunda_condicion, 'Value') == 2
            T_o_P_o_beta = str2double(get(handles.qf, 'String'));
            topob1 = 'T';
        elseif get(handles.segunda_condicion, 'Value') == 3
            T_o_P_o_beta = str2double(get(handles.qf, 'String'));
            topob1 = 'x';
        elseif get(handles.segunda_condicion, 'Value') == 4
            T_o_P_o_beta = str2double(get(handles.qf, 'String'));
            T_o_P_o_beta = abs(T_o_P_o_beta - 1);
            topob1 = 'x';
        end
        flujo = str2double(get(handles.F, 'String'));
        m_o_w_o_v = 'm';
        Corrienteppal = Corriente(mezcla, T_o_P_o_beta, topob1, vP111, 'P', flujo, m_o_w_o_v , MEdE111, 'C-000');
        entradas = [];
        try 
            entradas(end + 1) = mezcla.conc(1);
            entradas(end + 1) = Corrienteppal(end).q;
            entradas(end + 1) = flujo;
        catch
        end
        if strcmpi('Alimentacion',handles.uitable3.Data(1,5))
            try
                flujo = str2double(handles.uitable3.Data(1,1));
                conc = [str2double( handles.uitable3.Data(1,2) ), 1 - str2double( handles.uitable3.Data(1,2) )];
                if ~isempty(kij111)
                    mezcla = Mezcla(sustancia, conc, kij111);
                else
                    kij111 = 0;
                    mezcla = Mezcla(cmpbinario111, conc, kij111);
                end
                try 
                    T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(1,3) );
                    topob1 = 'x';
                Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, vP111, 'P', flujo, m_o_w_o_v , MEdE111, 'C-001')];
                catch ME
                    T_o_P_o_beta = str2double( handles.uitable3.Data(1,4) );
                    topob1 = 'T';
                Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, vP111, 'P', flujo, m_o_w_o_v , MEdE111, 'C-001')];
                end

                try 
                    entradas(end + 1) = conc(1);
                    entradas(end + 1) = Corrienteppal(end).q;
                    entradas(end + 1) = flujo;
                catch
                end
            catch
            end
        elseif strcmpi('Salida',handles.uitable3.Data(1,5))    
            try
                flujo = str2double(handles.uitable3.Data(1,1));
                conc = [str2double( handles.uitable3.Data(1,2) ), 1 - str2double( handles.uitable3.Data(1,2) )];
                try 
                    if abs((str2double( handles.uitable3.Data(1,3) )) - 1)<1e-5 || abs((str2double( handles.uitable3.Data(1,3) )))<1e-5
                        T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(1,3) );
                    else
                        error('No es adecuado el valor')
                    end
                catch
                    errordlg('En el caso de salidas solo puede utilizarse condicion térmica y debe ser 0 o 1');
                end
                salidas = {salidas{:}, conc, flujo, T_o_P_o_beta};
            catch
            end
        end

        if strcmpi('Alimentacion', handles.uitable3.Data(2,5) )
            try
                flujo = str2double( handles.uitable3.Data(2,1) );
                conc = [str2double( handles.uitable3.Data(2,2) ), 1 - str2double( handles.uitable3.Data(2,2) )];
                if ~isempty(kij111)
                    mezcla = Mezcla(sustancia, conc, kij111);
                else
                    kij111 = 0;
                    mezcla = Mezcla(cmpbinario111, conc, kij111);
                end
                try 
                    T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(2,3) );
                    topob1 = 'x';
                Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, vP111, 'P', flujo, m_o_w_o_v , MEdE111, 'C-002')];
                catch ME
                    T_o_P_o_beta = str2double( handles.uitable3.Data(2,4) );
                    topob1 = 'T';
                Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, vP111, 'P', flujo, m_o_w_o_v , MEdE111, 'C-002')];
                end

                try 
                    entradas(end + 1) = conc(1);
                    entradas(end + 1) = Corrienteppal(end).q;
                    entradas(end + 1) = flujo;
                catch
                end
            catch
            end
        elseif strcmpi('Salida',handles.uitable3.Data(2,5))    
            try
                flujo = str2double(handles.uitable3.Data(2,1));
                conc = [str2double( handles.uitable3.Data(2,2) ), 1 - str2double( handles.uitable3.Data(2,2) )];
                try 
                    if abs((str2double( handles.uitable3.Data(2,3) )) - 1)<1e-5 || abs((str2double( handles.uitable3.Data(2,3) )))<1e-5
                        T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(2,3) );
                    else
                        error('No es adecuado el valor')
                    end
                catch
                    errordlg('En el caso de salidas solo puede utilizarse condicion térmica y debe ser 0 o 1');
                end
                salidas = {salidas{:}, conc, flujo, T_o_P_o_beta};
            catch
            end
        end    
        if strcmpi('Alimentacion', handles.uitable3.Data(3,5) )
            try
                flujo = str2double( handles.uitable3.Data(3,1) );
                conc = [str2double( handles.uitable3.Data(3,2) ), 1 - str2double( handles.uitable3.Data(3,2) )];
                if ~isempty(kij111)
                    mezcla = Mezcla(sustancia, conc, kij111);
                else
                    kij111 = 0;
                    mezcla = Mezcla(cmpbinario111, conc, kij111);
                end
                try 
                    T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(3,3) );
                    topob1 = 'x';
                Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, vP111, 'P', flujo, m_o_w_o_v , MEdE111, 'C-002')];
                catch ME
                    T_o_P_o_beta = str2double( handles.uitable3.Data(3,4) );
                    topob1 = 'T';
                Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, vP111, 'P', flujo, m_o_w_o_v , MEdE111, 'C-002')];
                end
                try 
                    entradas(end + 1) = conc(1);
                    entradas(end + 1) = Corrienteppal(end).q;
                    entradas(end + 1) = flujo;
                catch
                end
            catch
            end
        elseif strcmpi('Salida',handles.uitable3.Data(3,5))    
            try
                flujo = str2double(handles.uitable3.Data(3,1));
                conc = [str2double( handles.uitable3.Data(3,2) ), 1 - str2double( handles.uitable3.Data(3,2) )];
                try 
                    if abs((str2double( handles.uitable3.Data(3,3) )) - 1)<1e-5 || abs((str2double( handles.uitable3.Data(3,3) )))<1e-5
                        T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(3,3) );
                    else
                        error('No es adecuado el valor')
                    end
                catch
                    errordlg('En el caso de salidas solo puede utilizarse condicion térmica y debe ser 0 o 1');
                end
                salidas = {salidas{:}, conc, flujo, T_o_P_o_beta};
            catch
            end
        end

        if strcmpi('Alimentacion', handles.uitable3.Data(4,5) )
            try
                flujo = str2double( handles.uitable3.Data(4,1) );
                conc = [str2double( handles.uitable3.Data(4,2) ), 1 - str2double( handles.uitable3.Data(4,2) )];
                if ~isempty(kij111)
                    mezcla = Mezcla(sustancia, conc, kij111);
                else
                    kij111 = 0;
                    mezcla = Mezcla(cmpbinario111, conc, kij111);
                end
                try 
                    T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(4,3) );
                    topob1 = 'x';
                Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, vP111, 'P', flujo, m_o_w_o_v , MEdE111, 'C-003')];
                catch ME
                    T_o_P_o_beta = str2double( handles.uitable3.Data(4,4) );
                    topob1 = 'T';
                Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, vP111, 'P', flujo, m_o_w_o_v , MEdE111, 'C-003')];
                end

                try 
                    entradas(end + 1) = conc(1);
                    entradas(end + 1) = Corrienteppal(end).q;
                    entradas(end + 1) = flujo;
                catch
                end
            catch
            end
        elseif strcmpi('Salida',handles.uitable3.Data(4,5))    
            try
                flujo = str2double(handles.uitable3.Data(4,1));
                conc = [str2double( handles.uitable3.Data(4,2) ), 1 - str2double( handles.uitable3.Data(4,2) )];
                try 
                    if abs((str2double( handles.uitable3.Data(4,3) )) - 1)<1e-5 || abs((str2double( handles.uitable3.Data(4,3) )))<1e-5
                        T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(4,3) );
                    else
                        error('No es adecuado el valor')
                    end
                catch
                    errordlg('En el caso de salidas solo puede utilizarse condicion térmica y debe ser 0 o 1');
                end
                salidas = {salidas{:}, conc, flujo, T_o_P_o_beta};
            catch
            end
        end
        B = 0;
        D = 0;
        xD = str2double(get(handles.Xd, 'String'));
        xB = str2double(get(handles.Xw, 'String'));
            vapvivo = get(handles.radiobutton10, 'Value');
        if get(handles.composic_dest, 'Value') == 1 

            if get(handles.composic_res, 'Value') == 1
                MccabeObject = Hengstebeck(Corrienteppal, 1, [], 2, [], salidas, [], [], [], [], [], vapvivo, get(handles.radiobutton4, 'Value'));
                MccabeObject.xDi =  str2double(get(handles.Xd, 'String'));
                MccabeObject.xBi =  str2double(get(handles.Xw, 'String'));
                for i = 1:length(Corrienteppal)
                    B = B + (Corrienteppal(i).molF*(xD - Corrienteppal(i).conc(1)));
                    D = D + (Corrienteppal(i).molF*(Corrienteppal(i).conc(1) - xB));
                end
                D = D/(xD - xB);
                B = B/(xD - xB);
                MccabeObject.Di = [D*xD, D*(1 - xD)];
                MccabeObject.Bi = [B*xB, B*(1 - xB)];
            else
                MccabeObject = Hengstebeck(Corrienteppal, 1, [], 2, str2double(get(handles.Xw, 'String')),  salidas, [], [], [], [], [], vapvivo, get(handles.radiobutton4, 'Value'));
                for i = 1:length(Corrienteppal)
                    B = B + (Corrienteppal(i).molF);
                    D = D + (Corrienteppal(i).molF*(Corrienteppal(i).conc(1) - xB));
                end
                F = B;
                B = xB*D;
                D = D/(xD - xB);
                B = B/(xD - xB);
                MccabeObject.Di = [D*xD, D*(1 - xD)];
                MccabeObject.Bi = [B*xB, B*(1 - xB)];
                MccabeObject.xDi =  str2double(get(handles.Xd, 'String'));
            end
        else
            if get(handles.composic_res, 'Value') == 1
                MccabeObject = Hengstebeck(Corrienteppal, 1, str2double(get(handles.Xd, 'String')), 2, [],  salidas, [], [], [], [], [], vapvivo, get(handles.radiobutton4, 'Value'));
                MccabeObject.xBi =  str2double(get(handles.Xw, 'String'));
            else
                MccabeObject = Hengstebeck(Corrienteppal, 1, str2double(get(handles.Xd, 'String')), 2, str2double(get(handles.Xw, 'String')), salidas, Re111, [], [], [], Effic, vapvivo, get(handles.radiobutton4, 'Value'));
            end

        end
        MccabeObject.entradas = entradas;
        MccabeObject.E_x = X111;
        MccabeObject.E_y = Y111;
        MccabeObject.Tbin = T111;
        if ~isempty(Effic)
            MccabeObject.Murphree = Effic/100;
        else
            MccabeObject.Murphree = 1;
        end
        MccabeObject.Mur_x = XE111;
        MccabeObject.Mur_y = YE111;
        MccabeObject.alter_define(MccabeObject.alimentaciones, 1, MccabeObject.reco_lk, 2, MccabeObject.reco_hk, MccabeObject.E_x);
        MccabeObject.distributemin;
        MccabeObject.graficamin();
        set(handles.nmin, 'String', MccabeObject.platosaefficmin)
        set(handles.eteo, 'String', MccabeObject.platosmin)
        if get(handles.radiobutton6, 'Value') == 1
            set(handles.ereal, 'String', ceil(MccabeObject.platosmin./((str2double(get(handles.E, 'String'))/100))))
        else
            set(handles.ereal, 'String', ceil(MccabeObject.etapas))
        end
        set(handles.dest, 'String', sum(MccabeObject.Di));
        set(handles.resi, 'String', sum(MccabeObject.Bi));
        set(handles.re_tope, 'String', MccabeObject.Rmin);
        set(handles.retope, 'String', MccabeObject.Rmin.*sum(MccabeObject.Di)./(sum(MccabeObject.Di).*(1 + MccabeObject.Rmin)));
        set(handles.re_fondo, 'String', MccabeObject.vap_min);
        handles.uitable4.Data(1, 1) = {'Inicial'};
        handles.uitable4.Data(1, 2) = {MccabeObject.rectasminima(1,1)};
        handles.uitable4.Data(1,3) = {MccabeObject.rectasminima(2,1)};
        handles.uitable4.Data(2,1) = {'Final'};
        handles.uitable4.Data(2,2) = {MccabeObject.rectasminima(1,2)};
        handles.uitable4.Data(2,3) = {MccabeObject.rectasminima(2,2)};
        end
    elseif get(handles.modelotermo, 'Value') == 6 
        if get(handles.radiobutton6,'Value') == 1
            set(handles.radiobutton5, 'Value', 0);
            Effic = 100;
        end
        cmpbinario111 = [Sustancia('Ethane'), Sustancia('Propane')];
        sustancia = cmpbinario111;
        conc = [Zf111, 1 - Zf111];
        if ~isempty(kij111)
            mezcla = Mezcla(sustancia, conc, kij111);
        else
            kij111 = 0;
            mezcla = Mezcla(cmpbinario111, conc, kij111);
        end
        if get(handles.segunda_condicion, 'Value') == 2
            errordlg('No se puede definir con temperatura')
        elseif get(handles.segunda_condicion, 'Value') == 3
            T_o_P_o_beta = str2double(get(handles.qf, 'String'));
            topob1 = 'x';
        elseif get(handles.segunda_condicion, 'Value') == 4
            T_o_P_o_beta = str2double(get(handles.qf, 'String'));
            T_o_P_o_beta = abs(T_o_P_o_beta - 1);
            topob1 = 'x';
        end
        flujo = str2double(get(handles.F, 'String'));
        m_o_w_o_v = 'm';
        Corrienteppal = Corriente(mezcla, T_o_P_o_beta, topob1, 100, 'P', flujo, m_o_w_o_v , MEdE111, 'C-000');
        entradas = [];
        try 
            entradas(end + 1) = mezcla.conc(1);
            entradas(end + 1) = Corrienteppal(end).q;
            entradas(end + 1) = flujo;
        catch
        end
        if strcmpi('Alimentacion',handles.uitable3.Data(1,5))
            try
                flujo = str2double(handles.uitable3.Data(1,1));
                conc = [str2double( handles.uitable3.Data(1,2) ), 1 - str2double( handles.uitable3.Data(1,2) )];
                if ~isempty(kij111)
                    mezcla = Mezcla(sustancia, conc, kij111);
                else
                    kij111 = 0;
                    mezcla = Mezcla(cmpbinario111, conc, kij111);
                end
                try 
                    T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(1,3) );
                    topob1 = 'x';
                Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, 100, 'P', flujo, m_o_w_o_v , MEdE111, 'C-001')];
                catch ME
                    T_o_P_o_beta = str2double( handles.uitable3.Data(1,4) );
                    topob1 = 'T';
                Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, 100, 'P', flujo, m_o_w_o_v , MEdE111, 'C-001')];
                end

                try 
                    entradas(end + 1) = conc(1);
                    entradas(end + 1) = Corrienteppal(end).q;
                    entradas(end + 1) = flujo;
                catch
                end
            catch
            end
        elseif strcmpi('Salida',handles.uitable3.Data(1,5))    
            try
                flujo = str2double(handles.uitable3.Data(1,1));
                conc = [str2double( handles.uitable3.Data(1,2) ), 1 - str2double( handles.uitable3.Data(1,2) )];
                try 
                    if abs((str2double( handles.uitable3.Data(1,3) )) - 1)<1e-5 || abs((str2double( handles.uitable3.Data(1,3) )))<1e-5
                        T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(1,3) );
                    else
                        error('No es adecuado el valor')
                    end
                catch
                    errordlg('En el caso de salidas solo puede utilizarse condicion térmica y debe ser 0 o 1');
                end
                salidas = {salidas{:}, conc, flujo, T_o_P_o_beta};
            catch
            end
        end

        if strcmpi('Alimentacion', handles.uitable3.Data(2,5) )
            try
                flujo = str2double( handles.uitable3.Data(2,1) );
                conc = [str2double( handles.uitable3.Data(2,2) ), 1 - str2double( handles.uitable3.Data(2,2) )];
                if ~isempty(kij111)
                    mezcla = Mezcla(sustancia, conc, kij111);
                else
                    kij111 = 0;
                    mezcla = Mezcla(cmpbinario111, conc, kij111);
                end
                try 
                    T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(2,3) );
                    topob1 = 'x';
                Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, 100, 'P', flujo, m_o_w_o_v , MEdE111, 'C-002')];
                catch ME
                    T_o_P_o_beta = str2double( handles.uitable3.Data(2,4) );
                    topob1 = 'T';
                Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, 100, 'P', flujo, m_o_w_o_v , MEdE111, 'C-002')];
                end

                try 
                    entradas(end + 1) = conc(1);
                    entradas(end + 1) = Corrienteppal(end).q;
                    entradas(end + 1) = flujo;
                catch
                end
            catch
            end
        elseif strcmpi('Salida',handles.uitable3.Data(2,5))    
            try
                flujo = str2double(handles.uitable3.Data(2,1));
                conc = [str2double( handles.uitable3.Data(2,2) ), 1 - str2double( handles.uitable3.Data(2,2) )];
                try 
                    if abs((str2double( handles.uitable3.Data(2,3) )) - 1)<1e-5 || abs((str2double( handles.uitable3.Data(2,3) )))<1e-5
                        T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(2,3) );
                    else
                        error('No es adecuado el valor')
                    end
                catch
                    errordlg('En el caso de salidas solo puede utilizarse condicion térmica y debe ser 0 o 1');
                end
                salidas = {salidas{:}, conc, flujo, T_o_P_o_beta};
            catch
            end
        end    
        if strcmpi('Alimentacion', handles.uitable3.Data(3,5) )
            try
                flujo = str2double( handles.uitable3.Data(3,1) );
                conc = [str2double( handles.uitable3.Data(3,2) ), 1 - str2double( handles.uitable3.Data(3,2) )];
                if ~isempty(kij111)
                    mezcla = Mezcla(sustancia, conc, kij111);
                else
                    kij111 = 0;
                    mezcla = Mezcla(cmpbinario111, conc, kij111);
                end
                try 
                    T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(3,3) );
                    topob1 = 'x';
                Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, 100, 'P', flujo, m_o_w_o_v , MEdE111, 'C-002')];
                catch ME
                    T_o_P_o_beta = str2double( handles.uitable3.Data(3,4) );
                    topob1 = 'T';
                Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, 100, 'P', flujo, m_o_w_o_v , MEdE111, 'C-002')];
                end
                try 
                    entradas(end + 1) = conc(1);
                    entradas(end + 1) = Corrienteppal(end).q;
                    entradas(end + 1) = flujo;
                catch
                end
            catch
            end
        elseif strcmpi('Salida',handles.uitable3.Data(3,5))    
            try
                flujo = str2double(handles.uitable3.Data(3,1));
                conc = [str2double( handles.uitable3.Data(3,2) ), 1 - str2double( handles.uitable3.Data(3,2) )];
                try 
                    if abs((str2double( handles.uitable3.Data(3,3) )) - 1)<1e-5 || abs((str2double( handles.uitable3.Data(3,3) )))<1e-5
                        T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(3,3) );
                    else
                        error('No es adecuado el valor')
                    end
                catch
                    errordlg('En el caso de salidas solo puede utilizarse condicion térmica y debe ser 0 o 1');
                end
                salidas = {salidas{:}, conc, flujo, T_o_P_o_beta};
            catch
            end    
        end

        if strcmpi('Alimentacion', handles.uitable3.Data(4,5) )
            try
                flujo = str2double( handles.uitable3.Data(4,1) );
                conc = [str2double( handles.uitable3.Data(4,2) ), 1 - str2double( handles.uitable3.Data(4,2) )];
                if ~isempty(kij111)
                    mezcla = Mezcla(sustancia, conc, kij111);
                else
                    kij111 = 0;
                    mezcla = Mezcla(cmpbinario111, conc, kij111);
                end
                try 
                    T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(4,3) );
                    topob1 = 'x';
                Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, 100, 'P', flujo, m_o_w_o_v , MEdE111, 'C-003')];
                catch ME
                    T_o_P_o_beta = str2double( handles.uitable3.Data(4,4) );
                    topob1 = 'T';
                Corrienteppal = [Corrienteppal, Corriente(mezcla, T_o_P_o_beta, topob1, 100, 'P', flujo, m_o_w_o_v , MEdE111, 'C-003')];
                end

                try 
                    entradas(end + 1) = conc(1);
                    entradas(end + 1) = Corrienteppal(end).q;
                    entradas(end + 1) = flujo;
                catch
                end
            catch
            end
        elseif strcmpi('Salida',handles.uitable3.Data(4,5))    
            try
                flujo = str2double(handles.uitable3.Data(4,1));
                conc = [str2double( handles.uitable3.Data(4,2) ), 1 - str2double( handles.uitable3.Data(4,2) )];
                try 
                    if abs((str2double( handles.uitable3.Data(4,3) )) - 1)<1e-5 || abs((str2double( handles.uitable3.Data(4,3) )))<1e-5
                        T_o_P_o_beta = 1 - str2double( handles.uitable3.Data(4,3) );
                    else
                        error('No es adecuado el valor')
                    end
                catch
                    errordlg('En el caso de salidas solo puede utilizarse condicion térmica y debe ser 0 o 1');
                end
                salidas = {salidas{:}, conc, flujo, T_o_P_o_beta};
            catch
            end
        end
        B = 0;
        D = 0;
        xD = str2double(get(handles.Xd, 'String'));
        xB = str2double(get(handles.Xw, 'String'));
            vapvivo = get(handles.radiobutton10, 'Value');
        if get(handles.composic_dest, 'Value') == 1 

            if get(handles.composic_res, 'Value') == 1
                MccabeObject = Hengstebeck(Corrienteppal, 1, [], 2, [], salidas, [], [], [], [], [], vapvivo, get(handles.radiobutton4, 'Value'));
                MccabeObject.xDi =  str2double(get(handles.Xd, 'String'));
                MccabeObject.xBi =  str2double(get(handles.Xw, 'String'));
                for i = 1:length(Corrienteppal)
                    B = B + (Corrienteppal(i).molF*(xD - Corrienteppal(i).conc(1)));
                    D = D + (Corrienteppal(i).molF*(Corrienteppal(i).conc(1) - xB));
                end
                D = D/(xD - xB);
                B = B/(xD - xB);
                MccabeObject.Di = [D*xD, D*(1 - xD)];
                MccabeObject.Bi = [B*xB, B*(1 - xB)];
            else
                MccabeObject = Hengstebeck(Corrienteppal, 1, [], 2, str2double(get(handles.Xw, 'String')),  salidas, [], [], [], [], [], vapvivo, get(handles.radiobutton4, 'Value'));
                for i = 1:length(Corrienteppal)
                    B = B + (Corrienteppal(i).molF);
                    D = D + (Corrienteppal(i).molF*(Corrienteppal(i).conc(1) - xB));
                end
                F = B;
                B = xB*D;
                D = D/(xD - xB);
                B = B/(xD - xB);
                MccabeObject.Di = [D*xD, D*(1 - xD)];
                MccabeObject.Bi = [B*xB, B*(1 - xB)];
                MccabeObject.xDi =  str2double(get(handles.Xd, 'String'));
            end
        else
            if get(handles.composic_res, 'Value') == 1
                MccabeObject = Hengstebeck(Corrienteppal, 1, str2double(get(handles.Xd, 'String')), 2, [],  salidas, [], [], [], [], [], vapvivo, get(handles.radiobutton4, 'Value'));
                MccabeObject.xBi =  str2double(get(handles.Xw, 'String'));
            else
                MccabeObject = Hengstebeck(Corrienteppal, 1, str2double(get(handles.Xd, 'String')), 2, str2double(get(handles.Xw, 'String')), salidas, Re111, [], [], [], Effic, vapvivo, get(handles.radiobutton4, 'Value'));
            end

        end
        MccabeObject.entradas = entradas;
        MccabeObject.E_x = X111;
        MccabeObject.E_y = Y111;
        MccabeObject.Tbin = T111;
        if ~isempty(Effic)
            MccabeObject.Murphree = Effic/100;
        else
            MccabeObject.Murphree = 1;
        end
        MccabeObject.E_x = X111;
        MccabeObject.E_y = Y111;
        MccabeObject.Mur_x = X111;
        MccabeObject.Mur_y = (Y111 - X111)* Effic/100+X111;
        MccabeObject.distributemin;
        MccabeObject.graficamin();
        set(handles.nmin, 'String', MccabeObject.platosaefficmin)
        set(handles.eteo, 'String', MccabeObject.platosmin)
        if get(handles.radiobutton6, 'Value') == 1
            set(handles.ereal, 'String', ceil(MccabeObject.platosmin./((str2double(get(handles.E, 'String'))/100))))
        else
            set(handles.ereal, 'String', ceil(MccabeObject.etapas))
        end
        set(handles.dest, 'String', sum(MccabeObject.Di));
        set(handles.resi, 'String', sum(MccabeObject.Bi));
        set(handles.re_tope, 'String', MccabeObject.Rmin);
        set(handles.retope, 'String', MccabeObject.Rmin.*sum(MccabeObject.Di)./(sum(MccabeObject.Di).*(1 + MccabeObject.Rmin)));
        set(handles.re_fondo, 'String', MccabeObject.vap_min);
        handles.uitable4.Data(1, 1) = {'Inicial'};
        handles.uitable4.Data(1, 2) = {MccabeObject.rectasminima(1,1)};
        handles.uitable4.Data(1,3) = {MccabeObject.rectasminima(2,1)};
        handles.uitable4.Data(2,1) = {'Final'};
        handles.uitable4.Data(2,2) = {MccabeObject.rectasminima(1,2)};
        handles.uitable4.Data(2,3) = {MccabeObject.rectasminima(2,2)};
    end



function eteo_Callback(hObject, eventdata, handles)
% hObject    handle to eteo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eteo as text
%        str2double(get(hObject,'String')) returns contents of eteo as a double


% --- Executes during object creation, after setting all properties.
function eteo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eteo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ereal_Callback(hObject, eventdata, handles)
% hObject    handle to ereal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ereal as text
%        str2double(get(hObject,'String')) returns contents of ereal as a double


% --- Executes during object creation, after setting all properties.
function ereal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ereal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eopt_Callback(hObject, eventdata, handles)
% hObject    handle to eopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eopt as text
%        str2double(get(hObject,'String')) returns contents of eopt as a double


% --- Executes during object creation, after setting all properties.
function eopt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eopt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dest_Callback(hObject, eventdata, handles)
% hObject    handle to dest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dest as text
%        str2double(get(hObject,'String')) returns contents of dest as a double


% --- Executes during object creation, after setting all properties.
function dest_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function resi_Callback(hObject, eventdata, handles)
% hObject    handle to resi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of resi as text
%        str2double(get(hObject,'String')) returns contents of resi as a double


% --- Executes during object creation, after setting all properties.
function resi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit34_Callback(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit34 as text
%        str2double(get(hObject,'String')) returns contents of edit34 as a double


% --- Executes during object creation, after setting all properties.
function edit34_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
global X111; global T111; global Y111; global XE111; global YE111;
figure(2);
hold on
plot([X111(1), X111(end)], [X111(1), X111(end)], 'Color', 'black');
plot(X111, Y111, '--', 'Color', 'red', 'LineWidth', 1);
legend('', 'Curva de equilibrio', 'Location', 'southeast');
axis([-0.0002 1.0002 -0.0002 1.0002])
