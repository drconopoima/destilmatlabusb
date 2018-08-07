function varargout = Riguroso(varargin)
% RIGUROSO Código MATLAB para la ejecución de Riguroso.fig
%
%      Función de control de la interfaz gráfica para los métodos rigurosos 
%      y multicomponentes.      
%
%      RIGUROSO, by itself, creates a new RIGUROSO or raises the existing
%      singleton*.
%
%      H = RIGUROSO returns the handle to a new RIGUROSO or the handle to
%      the existing singleton*.
%
%      RIGUROSO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RIGUROSO.M with the given input arguments.
%
%      RIGUROSO('Property','Value',...) creates a new RIGUROSO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Riguroso_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Riguroso_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Riguroso

% Last Modified by GUIDE v2.5 04-Sep-2015 11:25:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Riguroso_OpeningFcn, ...
                   'gui_OutputFcn',  @Riguroso_OutputFcn, ...
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


% --- Executes just before Riguroso is made visible.
function Riguroso_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Riguroso (see VARARGIN)

% Choose default command line output for Riguroso
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Riguroso wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Riguroso_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in cual_metodo.
function cual_metodo_Callback(hObject, eventdata, handles)
% hObject    handle to cual_metodo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cual_metodo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cual_metodo
global Torre1; global salidas; global damping

if isempty(salidas)
    salidas = cell.empty(0,1);
end
metodonumero = get(handles.cual_metodo, 'Value');
set(handles.axes1, 'Visible', 'on')
I = imread('torre.jpg');
imshow(I);
if  metodonumero ~= 1
    set(handles.cual_metodo, 'ForegroundColor', [0,0,0]);
    if metodonumero == 2
        Torre1 = BubblePoint();        
        set(handles.newtonraphsondata, 'Visible', 'off');
        set(handles.bubblepointdata, 'Visible', 'on');
        set(handles.panelespecif, 'Visible', 'on');
        set(handles.uipanel3, 'Visible', 'on');
        set(handles.perfil_v_l, 'Visible', 'on');
        set(handles.perfilt_panel, 'Visible', 'on');
        set(handles.flujos_bi_di, 'Visible', 'on');
        set(handles.dampingtext, 'Visible', 'off');
        set(handles.dampingedit, 'Visible', 'off');
    elseif metodonumero == 3
        Torre1 = NaphtaliSandholm();
        set(handles.newtonraphsondata, 'Visible', 'off')
        set(handles.bubblepointdata, 'Visible', 'off');
        set(handles.panelespecif, 'Visible', 'off');
        set(handles.uipanel3, 'Visible', 'on');
        set(handles.perfilt_panel, 'Visible', 'on');
        set(handles.perfil_v_l, 'Visible', 'on');
        set(handles.flujos_bi_di, 'Visible', 'on');
        set(handles.newtonraphsondata, 'Visible','on') 
        set(handles.dampingtext, 'Visible', 'on');
        damping = 0.6;
        set(handles.dampingedit, 'String', num2str(0.6));
        set(handles.dampingedit, 'Visible', 'on');
        set(handles.destiladofase,'Visible', 'on');
        set(handles.condiciondob2, 'ForegroundColor', [1,0,0])
        set(handles.reflujotext2, 'ForegroundColor', [1,0,0]);
        set(handles.panelespecif, 'ForegroundColor', [1,0,0]);
    end
    
else
    set(handles.cual_metodo, 'ForegroundColor', [1,0,0]);    
end



% --- Executes during object creation, after setting all properties.
function cual_metodo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cual_metodo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on key press with focus on cual_metodo and none of its controls.
function cual_metodo_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to cual_metodo (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
global prev_pc; global prev_tc; global prev_wa; global condicion1;
global condicion2;
global MulticompEdE; 
set(handles.popupmenu1, 'ForegroundColor', [0,0,0]);
condicion1 = 'T';
condicion2 = 'P';
correr=get(handles.popupmenu1, 'Value');
if correr == 2 || correr == 3 || correr == 4 || correr == 5
    set(handles.pc, 'ForegroundColor', [1,0,0]);
    set(handles.tc, 'ForegroundColor', [1,0,0]);
    set(handles.wa, 'ForegroundColor', [1,0,0]);
    set(handles.tc, 'Visible', 'on');
    set(handles.pc, 'Visible', 'on');
    set(handles.wa, 'Visible', 'on');    
    if correr == 2
        MulticompEdE = IdealEdE();
        set(handles.wa  , 'Visible', 'off');
    elseif correr == 3
        MulticompEdE = SRKEdE();
        MulticompEdE = RMVdW(MulticompEdE);
    elseif correr == 4
        MulticompEdE = PREdE();
        MulticompEdE = RMVdW(MulticompEdE);
    elseif correr == 5
        MulticompEdE = PRGEdE();
        MulticompEdE = RMVdW(MulticompEdE);
    end
elseif correr == 1
    set(handles.popupmenu1, 'ForegroundColor', [1,0,0]);
end


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in sustancias_lista.
function sustancias_lista_Callback(hObject, eventdata, handles)
% hObject    handle to sustancias_lista (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sustancias_lista contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sustancias_lista
global ident_sust; global formula; global masamolar;
global Tc; global Pc; global wa; global valor;
    comp2=get(handles.sustancias_lista,'Value');
    
    if comp2 ~= 2 && comp2 ~= 1
        set(handles.wa, 'Visible', 'on');
        load('indices1.mat'); %Carga la base de datos de sustancias
        ident_sust = database1(comp2-2);
        propiedades = get_properties(ident_sust);
        formula = propiedades{1};
        masamolar = propiedades{2};
        comp_repr = {strcat(ident_sust, ' ', '(', formula, ')')};
        Pc = propiedades{4};
        set(handles.pc, 'ForegroundColor', [0,0,0]); % Pc de rojo a negro  
        pc_str = strcat('Pc = ', num2str(Pc));
        set(handles.pc,'String', pc_str)
        Tc = propiedades{3};
        tc_str = strcat('Tc = ', num2str(Tc));
        set(handles.tc, 'ForegroundColor', [0,0,0]); % Tc de rojo a negro
        set(handles.tc,'String', tc_str)
        wa = propiedades{5};
        wa_str = strcat('w = ', num2str(wa));
        set(handles.wa, 'ForegroundColor', [0,0,0]); % wa de rojo a negro
        set(handles.wa,'String', wa_str)
        %%% Pasan a mostrarse los datos críticos de la sustancia
        %%% seleccionada
    elseif comp2 == 1
        set(handles.pc,'ForegroundColor',[1 0 0])
        pc_str = 'Pc (borrar y entrar)';
        set(handles.pc, 'String', pc_str);
        set(handles.tc, 'ForegroundColor', [1,0,0]);
        tc_str = 'Tc (borrar y entrar)';
        set(handles.tc, 'String', tc_str);
        set(handles.wa, 'ForegroundColor', [1,0,0]);
        wa_str = 'w (borrar y entrar)';
        set(handles.wa, 'String', wa_str);
        clearvars ident_sust; clearvars formula; clearvars masamolar; clearvars Tc; clearvars Pc; clearvars wa
        
    else
        set(handles.pc,'ForegroundColor',[1 0 0])
        pc_str = 'Pc (borrar y entrar)';
        set(handles.pc, 'String', pc_str);
        set(handles.tc, 'ForegroundColor', [1,0,0]);
        tc_str = 'Tc (borrar y entrar)';
        set(handles.tc, 'String', tc_str);
        set(handles.wa, 'ForegroundColor', [1,0,0]);
        wa_str = 'w (borrar y entrar)';
        set(handles.wa, 'String', wa_str);
        clearvars ident_sust; clearvars formula; clearvars masamolar; clearvars Tc; clearvars Pc; clearvars wa
        
    end

% --- Executes during object creation, after setting all properties.
function sustancias_lista_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sustancias_lista (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in componentes_visible.
function componentes_visible_Callback(hObject, eventdata, handles)
% hObject    handle to componentes_visible (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns componentes_visible contents as cell array
%        contents{get(hObject,'Value')} returns selected item from componentes_visible


% --- Executes during object creation, after setting all properties.
function componentes_visible_CreateFcn(hObject, eventdata, handles)
% hObject    handle to componentes_visible (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in anadir_sustancia.
function anadir_sustancia_Callback(hObject, eventdata, handles)
% hObject    handle to anadir_sustancia (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Tc; global Pc; global ident_sust; global masamolar; global formula
global wa; global Sust; global valor; global userdef;
if get(handles.sustancias_lista, 'Value') ~= 1 
        
    if ~isempty(valor)
        valor = valor + 1;
    else
        valor = 1;
    end
    if isa(Sust, 'double')
        Sust = Sustancia.empty(0);
    end

    if get(handles.sustancias_lista, 'Value') ~= 2
        comp_repr = {strcat(num2str(valor), '.', ' ', ident_sust, ' ', '(', formula, ')')};
    else 
        if isempty(userdef)
            userdef = 1;
        else userdef = userdef + 1;
        end
        ident_sust = strcat('UserDef', '.',num2str(userdef));
        formula = '';
        comp_repr = {strcat(num2str(valor), '.', ' ', ident_sust, ' ', '(', formula, ')')};
        set(handles.tc, 'String', 'Tc (borrar y entrar)');
        set(handles.pc, 'String', 'Pc (borrar y entrar)');
        set(handles.wa, 'String', 'w (borrar y entrar)');       
        set(handles.tc, 'ForegroundColor', [1,0,0]);
        set(handles.pc, 'ForegroundColor', [1,0,0]);
        set(handles.wa, 'ForegroundColor', [1,0,0]);
        Tc = [];
        Pc = [];
        wa = [];
    end
   
    
    Sust(end + 1) = Sustancia(ident_sust, masamolar, Tc, Pc, wa, [], [], [], [], formula);
    if length(Sust) >= 2
        set(handles.componentes_visible, 'ForegroundColor', [0,0,0]);
        set(handles.anadir_sustancia, 'ForegroundColor', [0,0,0]);
        set(handles.sustancias_lista, 'ForegroundColor', [0,0,0]);
    end
    
    prev_items = get(handles.componentes_visible, 'String'); 

    table_repr = strcat(ident_sust, ' ', '(', formula, ')');

    handles.uitable1.Data{valor, 1} = table_repr;

    if valor == 1
        prev_items = comp_repr;
    else
        prev_items(valor) = comp_repr;
    end
    set(handles.componentes_visible, 'Max', valor)
    set(handles.componentes_visible, 'String', prev_items, 'Value', valor)
end

% --- Executes on button press in remover_sust.
function remover_sust_Callback(hObject, eventdata, handles)
% hObject    handle to remover_sust (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Sust; global valor

removervalor = get(handles.componentes_visible, 'Value');
prev_items = get(handles.componentes_visible, 'String');
Sust(removervalor) = []; 
handles.uitable1.Data{valor,2} = '';
handles.uitable1.Data{valor,1} = '';
valor = valor - 1;
if valor < 2
    set(handles.popupmenu1, 'ForegroundColor', [1,0,0]);
    set(handles.anadir_sustancia, 'ForegroundColor', [1,0,0]);
    set(handles.sustancias_lista, 'ForegroundColor', [1,0,0]);
end
prev_items(removervalor) = [];
set(handles.componentes_visible, 'Value', 1);
set(handles.componentes_visible, 'String', prev_items);



function tc_Callback(hObject, eventdata, handles)
% hObject    handle to tc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Tc; 
Tc = str2double(get(handles.tc, 'String'));
tc_str = strcat('Tc = ',get(handles.tc, 'String'));
set(handles.tc, 'String', tc_str);
set(handles.tc, 'ForegroundColor', [0,0,0]);





% --- Executes during object creation, after setting all properties.
function tc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pc_Callback(hObject, eventdata, handles)
% hObject    handle to pc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Pc; 
Pc = str2double(get(handles.pc, 'String'));
pc_str = strcat('Pc = ',get(handles.pc, 'String'));
set(handles.pc, 'String', pc_str);
set(handles.pc, 'ForegroundColor', [0,0,0]);


% --- Executes during object creation, after setting all properties.
function pc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wa_Callback(hObject, eventdata, handles)
% hObject    handle to wa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global wa; 
wa = str2double(get(handles.wa, 'String'));
wa_str = strcat('wa = ',get(handles.wa, 'String'));
set(handles.wa, 'String', wa_str);
set(handles.wa, 'ForegroundColor', [0,0,0]);



% --- Executes during object creation, after setting all properties.
function wa_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in avanzadas.
function avanzadas_Callback(hObject, eventdata, handles)
% hObject    handle to avanzadas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Sust;    global kij;

if ~isempty(Sust) && length(Sust) >= 2
    set(handles.componentes_visible, 'Visible', 'off')
    set(handles.atras, 'Visible', 'on')
    set(handles.uitable2, 'Visible', 'on')
    mezclaparcial = Mezcla(Sust, ones(1, length(Sust)));
    if isempty(kij)
        kij = zeros(1, length(mezclaparcial.kbinario)); 
    end
    for ii = 1:length(mezclaparcial.kbinario)
        handles.uitable2.Data{ii,1} = mezclaparcial.kbinario{ii};
    end
end


% --- Executes on selection change in primera_condicion.
function primera_condicion_Callback(hObject, eventdata, handles)
% hObject    handle to primera_condicion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns primera_condicion contents as cell array
%        contents{get(hObject,'Value')} returns selected item from primera_condicion
global condicion1
whichisselected = get(handles.primera_condicion, 'Value');
if whichisselected == 1
    condicion1 = 'T';
elseif whichisselected == 2
    condicion1 = 'P';
elseif whichisselected == 3
    condicion1 = 'x';
elseif whichisselected == 4
    condicion1 = 'x';
end


% --- Executes during object creation, after setting all properties.
function primera_condicion_CreateFcn(hObject, eventdata, handles)
% hObject    handle to primera_condicion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in segunda_condicion.
function segunda_condicion_Callback(hObject, eventdata, handles)
% hObject    handle to segunda_condicion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global condicion2;
whichisselected = get(handles.segunda_condicion, 'Value');
if whichisselected == 1
    condicion2 = 'P';
elseif whichisselected == 2
    condicion2 = 'x';
elseif whichisselected == 3
    condicion2 = 'x';
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



function topoq1_Callback(hObject, eventdata, handles)
% hObject    handle to topoq1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of topoq1 as text
%        str2double(get(hObject,'String')) returns contents of topoq1 as a double
 global T_o_P_o_beta;
 
 if get(handles.primera_condicion, 'Value') ~= 3
    T_o_P_o_beta = str2double(get(handles.topoq1, 'String')); 
 else
    T_o_P_o_beta = 1 - str2double(get(handles.topoq1, 'String'));
 end
 
if str2double(get(handles.topoq1, 'String')) >=0
    set(handles.primera_condicion, 'ForegroundColor', [0,0,0])
end

% --- Executes during object creation, after setting all properties.
function topoq1_CreateFcn(hObject, eventdata, handles)
% hObject    handle toanatopoq1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function poqox2_Callback(hObject, eventdata, handles)
% hObject    handle to poqox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of poqox2 as text
%        str2double(get(hObject,'String')) returns contents of poqox2 as a double
global P_o_beta_o_q

if get(handles.segunda_condicion, 'Value') ~= 2
    P_o_beta_o_q = str2double(get(handles.poqox2, 'String')); 
else
    P_o_beta_o_q = 1 - str2double(get(handles.poqox2, 'String'));
end
if str2double(get(handles.poqox2, 'String')) >=0
    set(handles.segunda_condicion, 'ForegroundColor', [0,0,0])
end
 

% --- Executes during object creation, after setting all properties.
function poqox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to poqox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
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
global F

F = str2double(get(handles.F, 'String'));

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


% --- Executes on button press in calcular_corriente.
function calcular_corriente_Callback(hObject, eventdata, handles)
% hObject    handle to calcular_corriente (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global valor; global Sust; global MulticompEdE; global T_o_P_o_beta; 
global P_o_beta_o_q; global condicion1;
global condicion2; global kij; global F; global corriente1; global mezcla1; 
global conci; global numer_corrientes; global ident_corriente;
global alim; global etapasalim;

etapacorriente = str2double(get(handles.etapaplato, 'String'));
if ~isnan(etapacorriente)
    if get(handles.normalizar, 'Value') == 1 
        Ftotal = 0;
        conci = zeros(1, valor);

        for ittrety = 1:valor        
            flujoi = str2double(handles.uitable1.Data{ittrety, 2});
            conci(ittrety) = flujoi;
            handles.uitable1.Data{ittrety, 2} = num2str(conci(ittrety));
            Ftotal = Ftotal + flujoi;
        end
        conci = conci ./ Ftotal;
        ceroletras = '';
        if strcmpi(get(handles.F, 'String'), ceroletras)
            set(handles.F, 'String', num2str(Ftotal))
        else 
            Ftotal = str2double(get(handles.F, 'String'));
        end    
        if isempty(F)
            F = Ftotal;
        end
    else
        conci = zeros(1, valor);
        for ittrety = 1:valor
            conci(ittrety) = str2double(handles.uitable1.Data{ittrety, 2});
        end
        Ftotal = str2double(get(handles.F, 'String'));

        if isempty(F)
            F = Ftotal;
        end
        
    end
    try     
        colormodelotermo = get(handles.popupmenu1, 'ForegroundColor');
        colorsustancias_lista = get(handles.sustancias_lista, 'ForegroundColor');
        colorprimera_condicion = get(handles.primera_condicion, 'ForegroundColor');
        colorsegunda_condicion = get(handles.segunda_condicion, 'ForegroundColor');
        coloranadir_sustancia = get(handles.anadir_sustancia, 'ForegroundColor');
        negro = [0,0,0];
        if all(colormodelotermo == negro) && all(colorsustancias_lista == negro) && all(colorprimera_condicion == negro) && all(colorsegunda_condicion == negro) && all(coloranadir_sustancia == negro)
            if isempty(mezcla1)
                mezcla1 = Mezcla.empty(0,1);
            end
            if isempty(numer_corrientes)
                numer_corrientes = 1;
            else
                numer_corrientes = numer_corrientes + 1;
            end
            if isempty(corriente1)
                corriente1 = Corriente.empty(0,1);
            end 
            if isempty(ident_corriente) || ~isa(ident_corriente, 'cell')
                ident_corriente = cell.empty(0,1);
            end
            mezcla1(numer_corrientes) = Mezcla(Sust, conci, kij);
            display(mezcla1);
            ident_corriente{numer_corrientes} = strcat('S-', '0', num2str(numer_corrientes));
            corriente1(numer_corrientes) = Corriente(mezcla1(numer_corrientes), T_o_P_o_beta, condicion1, P_o_beta_o_q, condicion2, F, 'm', MulticompEdE, ident_corriente{numer_corrientes});
            display(corriente1);
            if isempty(alim)
                alim = cell.empty(1,0);
            end
            alim{end+1} = etapacorriente;
            alim{end+1} = corriente1(numer_corrientes);
            handles.uitable5.Data{numer_corrientes, 1} = ident_corriente{numer_corrientes};
            handles.uitable5.Data{numer_corrientes, 2} = num2str(etapacorriente);
            handles.uitable5.Data{numer_corrientes, 3} = 'Alimentación';
            handles.uitable5.Data{numer_corrientes, 4} = num2str(corriente1(numer_corrientes).T);
            handles.uitable5.Data{numer_corrientes, 5} = num2str(corriente1(numer_corrientes).P);
            handles.uitable5.Data{numer_corrientes, 6} = num2str(corriente1(numer_corrientes).q);
            handles.uitable5.Data{numer_corrientes, 7} = num2str(corriente1(numer_corrientes).H);
            handles.presionycalor.Data{2, 2} = num2str(0);
            handles.presionycalor.Data{3, 2} = num2str(0);
            handles.presionycalor.Data{1, 1} = num2str( num2str(corriente1(numer_corrientes).P));
            handles.presionycalor.Data{2, 1}= num2str( num2str(corriente1(numer_corrientes).P));
            handles.presionycalor.Data{3, 1} = num2str( num2str(corriente1(numer_corrientes).P));
            handles.presionycalor.Data{4, 1} = num2str( num2str(corriente1(numer_corrientes).P));
            for itterew = 1:corriente1.num_sust
                handles.flujosproductocomp.Data{1, itterew} = 0;
                handles.flujosproductocomp.Data{2, itterew} = 0;
            end
        else
            errordlg('Faltan algunos datos por rellenar, por favor revise áreas en rojo')
        end
    catch
        errordlg('Se ha producido un error inesperado. Intente rellenar áreas en rojo o reintroducir datos.')
        mezcla1 = Mezcla.empty(0,1);
        ident_corriente = cell.empty(0,1);
        corriente1 = Corriente.empty(0,1);
        alim = cell.empty(0,1);
        numer_corrientes = numer_corrientes - 1;
    end
    if length(corriente1) == 2
        I = imread('torre2.jpg');
        imshow(I);
    elseif length(corriente1) == 3
        I = imread('torre3.jpg');
        imshow(I);
    elseif length(corriente1) == 4
        I = imread('torre4.jpg');
        imshow(I);
    end
end
% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
global F; global Fi; global Sust;
if isempty(F)
    F = 0;    
    Fi = zeros(1, length(Sust));
end
ceroletras = '';
indices = eventdata.Indices;
if indices(2) == 2
    if strcmpi(eventdata.PreviousData, ceroletras)
        nuevoFi = str2double(eventdata.NewData);        
            F = F + nuevoFi;  
            Fi(indices(1)) = nuevoFi;        
        set(handles.F, 'String', num2str(F));
    else
        viejoFi = str2double(eventdata.PreviousData);
        nuevoFi = str2double(eventdata.NewData);
        F = F - viejoFi + nuevoFi;
        set(handles.F, 'String', num2str(F));
        Fi(indices(1)) = nuevoFi ;
    end
end


% --- Executes on button press in normalizar.
function normalizar_Callback(hObject, eventdata, handles)
% hObject    handle to normalizar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and useranadir data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of normalizar
if get(handles.normalizar, 'Value') == 1
    set(handles.normalizar, 'Value',0);
end


% --- Executes when entered data in editable cell(s) in uitable2.
function uitable2_CellEditCallback(hObject, eventdata, handles)
%uitable2 es la tabla popup que aparece para modificar parámetros de
%interacción binaria entre componentes
% hObject    handle to uitable2 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
global kij; global Sust;

mezclaparcial = Mezcla(Sust, ones(1,length(Sust)));
indices = eventdata.Indices;
indice = indices(1);
if isempty(kij)
    kij = zeros(1, length(mezclaparcial.kbinario));
end
if str2double(eventdata.NewData) < 0.5 && str2double(eventdata.NewData) > -0.5
    kij(indice) = str2double(eventdata.NewData);
else
    errordlg('El rango de valores para el parámetro de interacción binaria es -0.5<kij<0.5')
end    

% --- Executes on button press in atras.
function atras_Callback(hObject, eventdata, handles)
% hObject    handle to atras (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uitable2, 'Visible', 'off')
set(handles.atras, 'Visible', 'off')
set(handles.componentes_visible, 'Visible', 'on')
set(handles.avanzadas, 'Visible', 'on');


% --- Executes when entered data in editable cell(s) in presionycalor.
function presionycalor_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to presionycalor (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
global perfil_P; global perfil_q; global N

indice =  eventdata.Indices;
if indice(2) == 1
    perfil_P(indice(1)) = str2double(eventdata.EditData);
elseif indice(2) == 2 && indice(1) > 1 && indice(1) < N
    perfil_q(indice(1)-1) = str2double(eventdata.EditData);
elseif indice(2) == 2 && indice(1) == 1
    handles.presionycalor.Data{indice(1), indice(2)} = '';
    errordlg('En el método de punto de burbuja no se puede proveer calor estimado');
elseif indice(2) == 2 && indice(1) == N
    handles.presionycalor.Data{indice(1), indice(2)} = '';
    errordlg('En el método de punto de burbuja no se puede proveer calor estimado');
end


function etapaplato_Callback(hObject, eventdata, handles)
% hObject    handle to etapaplato (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etapaplato as text
%        str2double(get(hObject,'String')) returns contents of etapaplato as a double


% --- Executes during object creation, after setting all properties.
function etapaplato_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etapaplato (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in bubblepointdata.
function bubblepointdata_Callback(hObject, eventdata, handles)
% hObject    handle to bubblepointdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns bubblepointdata contents as cell array
%        contents{get(hObject,'Value')} returns selected item from bubblepointdata


% --- Executes during object creation, after setting all properties.
function bubblepointdata_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bubblepointdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in newtonraphsondata.
function newtonraphsondata_Callback(hObject, eventdata, handles)
% hObject    handle to newtonraphsondata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns newtonraphsondata contents as cell array
%        contents{get(hObject,'Value')} returns selected item from newtonraphsondata
condiciones = get(handles.newtonraphsondata, 'Value');

if condiciones == 1
    set(handles.panelespecif, 'Visible', 'off');
    set(handles.panelespecif2, 'Visible', 'off');
    set(handles.panelespecif3, 'Visible', 'off');
    set(handles.panelespecif4, 'Visible', 'off');
elseif condiciones == 2
    set(handles.panelespecif2, 'Visible', 'on');
    set(handles.panelespecif3, 'Visible', 'off');
    set(handles.panelespecif4, 'Visible', 'off');
    set(handles.panelespecif, 'Visible', 'on');
elseif condiciones == 3
    set(handles.panelespecif2, 'Visible', 'on');
    set(handles.panelespecif3, 'Visible', 'on');
    set(handles.panelespecif4, 'Visible', 'off');
    set(handles.panelespecif, 'Visible', 'on');
elseif condiciones == 4
    set(handles.panelespecif2, 'Visible', 'on');
    set(handles.panelespecif3, 'Visible', 'on');
    set(handles.panelespecif4, 'Visible', 'on');
    set(handles.panelespecif, 'Visible', 'on')
end

% --- Executes during object creation, after setting all properties.
function newtonraphsondata_CreateFcn(hObject, eventdata, handles)
% hObject    handle to newtonraphsondata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in destiladofase.
function destiladofase_Callback(hObject, eventdata, handles)
% hObject    handle to destiladofase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns destiladofase contents as cell array
%        contents{get(hObject,'Value')} returns selected item from destiladofase
global salidas; global Torre1;

if isa(Torre1, 'BubblePoint')
    if get(handles.destiladofase, 'Value') == 1
        salidas = cell.empty(1,0);
        set(handles.liquidodestilado, 'Visible', 'off')
    elseif get(handles.destiladofase, 'Value') == 3    
        set(handles.liquidodestilado, 'Visible', 'on')
    else
        if get(handles.condiciondbubble, 'Value') == 2
            set(handles.liquidodestilado, 'Visible', 'on')
        else
            set(handles.liquidodestilado, 'Visible', 'off')
        end
    end
elseif isa(Torre1, 'NaphtaliSandholm') && (get(handles.newtonraphsondata, 'Value') == 2)
    if get(handles.destiladofase, 'Value') == 1
        salidas = cell.empty(1,0);
        set(handles.liquidodestilado, 'Visible', 'off')
    elseif get(handles.destiladofase, 'Value') == 3    
        set(handles.liquidodestilado, 'Visible', 'on')
    else
        if get(handles.condiciondob2, 'Value') == 2
            set(handles.liquidodestilado, 'Visible', 'on')
        else
            set(handles.liquidodestilado, 'Visible', 'off')
        end
    end
else
    if get(handles.destiladofase, 'Value') ~= 1
        salidas = cell.empty(1,0);
        set(handles.liquidodestilado, 'Visible', 'on')
    end
end

% --- Executes during object creation, after setting all properties.
function destiladofase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to destiladofase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% Hints: get(hObject,'String') returns contents of nteoricas as text
%        str2double(get(hObject,'String')) returns contents of nteoricas as a double
function nteoricas_Callback(hObject, eventdata, handles)
% hObject    handle to nteoricas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global N; global perfil_P; global perfil_q; global corriente1; global numer_corrientes
global perfil_v; global perfil_vi; global perfil_l; global perfil_li;
global xsubj; global ysubj; global perfil_k; global bsubi; global dsubi;
global salidas;

if ~isempty(corriente1)
    N = str2double(get(handles.nteoricas, 'String'));
    if isempty(perfil_P)
        perfil_P = ones(1, N).*corriente1(1).P;
    end
    if isempty(perfil_q)
        perfil_q = zeros(1, N);
    end
    for ittere = 1:N
        handles.presionycalor.Data{ittere, 1} = num2str(corriente1(1).P);
        if ittere > 1 && ittere < N
            handles.presionycalor.Data{ittere, 2} = num2str(0);
        end
    end
    set(handles.panelnteoricas, 'ForegroundColor', [0,0,0]);
end


% --- Executes during object creation, after setting all properties.
function nteoricas_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nteoricas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function reflujoBP_Callback(hObject, eventdata, handles)
% hObject    handle to reflujoBP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of reflujoBP as text
%        str2double(get(hObject,'String')) returns contents of reflujoBP as a double
global R;

reflujo = str2double(get(handles.reflujoBP,'String'));
if ~isempty(reflujo) && isa(reflujo,'double') && (reflujo > 0)
    set(handles.reflujotext, 'ForegroundColor', [0,0,0])
    R = reflujo;
    colorcondiciondbubble = get(handles.condiciondbubble, 'ForegroundColor');
    colorpanelespecif = get(handles.condiciondbubble, 'ForegroundColor');
    negro = [0,0,0];
    if all(colorcondiciondbubble == negro) && all(colorpanelespecif == negro)
        set(handles.bubblepointdata, 'ForegroundColor', [0,0,0])
    end
else 
    set(handles.reflujotext, 'ForegroundColor', [1,0,0]);
    errordlg('Reflujo debe valor numérico ser mayor a 0');
end

% --- Executes during object creation, after setting all properties.
function reflujoBP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reflujoBP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in condiciondbubble.
function condiciondbubble_Callback(hObject, eventdata, handles)
% hObject    handle to condiciondbubble (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns condiciondbubble contents as cell array
%        contents{get(hObject,'Value')} returns selected item from condiciondbubble
global salidas

if get(handles.condiciondbubble, 'Value') == 2
    if get(handles.destiladofase, 'Value') == 2 || get(handles.destiladofase, 'Value') == 3

        set(handles.liquidodestilado, 'Visible', 'on')

    end
else
      set(handles.liquidodestilado, 'Visible', 'off')
end


% --- Executes during object creation, after setting all properties.
function condiciondbubble_CreateFcn(hObject, eventdata, handles)
% hObject    handle to condiciondbubble (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function flujod_o_b_Callback(hObject, eventdata, handles)
% hObject    handle to flujod_o_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of flujod_o_b as text
%        str2double(get(hObject,'String')) returns contents of flujod_o_b as a double
global D; global salidas; global dob

tipo_de_destilado = get(handles.destiladofase, 'Value'); 
if tipo_de_destilado == 1
    D = str2double(get(handles.flujod_o_b, 'String'));
    esdob = get(handles.condiciondbubble, 'Value');
    if esdob == 1
        dob = 'D';
    else
        dob = 'B';
    end
    set(handles.condiciondbubble,'ForegroundColor', [0,0,0]);
    set(handles.panelespecif, 'ForegroundColor', [0,0,0]);
    colorpanelespecif = get(handles.condiciondbubble, 'ForegroundColor');
    colorreflujotext = get(handles.reflujotext, 'ForegroundColor');
    negro = [0,0,0];
    if all(colorreflujotext == negro) && all(colorpanelespecif == negro)
        set(handles.bubblepointdata, 'ForegroundColor', [0,0,0])
    end
elseif tipo_de_destilado == 2
    D = str2double(get(handles.flujod_o_b, 'String'));
    esdob = get(handles.condiciondbubble, 'Value');
    if esdob == 1
        dob = 'D';
    else
        dob = 'B';
    end
    if esdob == 1
        salidas = {1,1,D};
        set(handles.panelespecif, 'ForegroundColor',[0,0,0])
    end
    set(handles.condiciondbubble, 'ForegroundColor',[0,0,0])
    colorpanelespecif = get(handles.condiciondbubble, 'ForegroundColor');
    colorreflujotext = get(handles.reflujotext, 'ForegroundColor');
    negro = [0,0,0];
    if all(colorreflujotext == negro) && all(colorpanelespecif == negro)
        set(handles.bubblepointdata, 'ForegroundColor', [0,0,0])
    end
else
    set(handles.liquidodestilado, 'Visible', 'on');
    D = str2double(get(handles.flujod_o_b, 'String'));
    esdob = get(handles.condiciondbubble, 'Value');
    if esdob == 1
        dob = 'D';
    else
        dob = 'B';
    end
    set(handles.condiciondbubble, 'ForegroundColor',[0,0,0])
    colorpanelespecif = get(handles.condiciondbubble, 'ForegroundColor');
    colorreflujotext = get(handles.reflujotext, 'ForegroundColor');
    negro = [0,0,0];
    if all(colorreflujotext == negro) && all(colorpanelespecif == negro)
        set(handles.bubblepointdata, 'ForegroundColor', [0,0,0])
    end
end


% --- Executes during object creation, after setting all properties.
function flujod_o_b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to flujod_o_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function liquidodestilado_Callback(hObject, eventdata, handles)
% hObject    handle to liquidodestilado (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of liquidodestilado as text
%        str2double(get(hObject,'String')) returns contents of liquidodestilado as a double
global D; global salidas

flujoliquido = str2double(get(handles.liquidodestilado, 'String'));
salidas = {1,1,flujoliquido};
set(handles.panelespecif, 'ForegroundColor', [0,0,0])

colorcondiciondbubble = get(handles.condiciondbubble, 'ForegroundColor');
colorreflujotext = get(handles.reflujotext, 'ForegroundColor');
negro = [0,0,0];
if all(colorcondiciondbubble == negro) && all(colorreflujotext == negro)
    set(handles.bubblepointdata, 'ForegroundColor', [0,0,0])
end



% --- Executes during object creation, after setting all properties.
function liquidodestilado_CreateFcn(hObject, eventdata, handles)
% hObject    handle to liquidodestilado (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in balanmasa.
function balanmasa_Callback(hObject, eventdata, handles)
% hObject    handle to balanmasa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global D; global dob; global salidas; global R; global corriente1; global Torre1;
global perfil_P; global perfil_q; global N; global alim; global iteraciones_max;
global supoc_di; global supoc_bi; global otradob; global damping;

if isempty(iteraciones_max)
    iteraciones_max = 100;
end

if isa(Torre1, 'BubblePoint')
    colorcondiciondbubble = get(handles.condiciondbubble, 'ForegroundColor');
    colorpanelespecif = get(handles.panelespecif, 'ForegroundColor');
    colorreflujotext = get(handles.reflujotext, 'ForegroundColor');
    negro = [0,0,0];
    
    if all(colorcondiciondbubble == negro) && all(colorpanelespecif == negro) && all(colorreflujotext == negro)
        try
            if all(perfil_P == corriente1(1).P) && all(perfil_q == 0)
                Torre1 = BubblePoint(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, 0, R, D, dob, iteraciones_max, 'T-001');

            elseif all(perfil_P == corriente1(1).P)

                Torre1 = BubblePoint(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, perfil_q, R, D, dob, iteraciones_max, 'T-001');
            elseif all(perfil_q == 0)

                Torre1 = BubblePoint(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, 0, R, D, dob, iteraciones_max, 'T-001');
            else

                Torre1 = BubblePoint(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, perfil_q, R, D, dob, iteraciones_max, 'T-001');
            end
        catch
            Torre1 = BubblePoint(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, perfil_q, R, D, dob, iteraciones_max, 'T-001');
        end
        Torre1.balanmasa();
        
        for ittreww = 1:N
            handles.perfil_v_l.Data{ittreww, 1} = num2str(Torre1.perfil_v(ittreww));
            handles.perfil_v_l.Data{ittreww, 2} = num2str(Torre1.perfil_l(ittreww));            
            handles.perfil_t.Data{ittreww, 1} = num2str(Torre1.perfil_t(ittreww));
        end
        for itterew = 1:length(alim{2}.mezcla.comp)
            handles.flujosproductocomp.Data{1, itterew} = num2str(Torre1.dsubi(itterew));
            handles.flujosproductocomp.Data{2, itterew} = num2str(Torre1.bsubi(itterew));
        end
    else
        errordlg('Por favor revise los datos faltantes en color rojo')
    end
    
elseif isa(Torre1, 'NaphtaliSandholm') && get(handles.newtonraphsondata, 'Value') == 2
    
    colorcondiciondob2 = get(handles.condiciondob2, 'ForegroundColor');
    colorpanelespecif = get(handles.panelespecif, 'ForegroundColor');
    colorreflujotext = get(handles.reflujotext2, 'ForegroundColor');
    negro = [0,0,0];
    
    if all(colorcondiciondob2 == negro) && all(colorpanelespecif == negro) && all(colorreflujotext == negro)
        try
            if all(perfil_P == corriente1(1).P) && all(perfil_q == 0)
                if (isempty(supoc_di) || all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, 0, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob);
                elseif (~isempty(supoc_di) && ~all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, 0, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob, supoc_di, 'di');
                elseif (isempty(supoc_di) || all(supoc_di == 0)) && (~isempty(supoc_bi) && ~all(supoc_bi == 0))   
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, 0, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob, supoc_bi, 'bi');
                else
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, 0, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob, supoc_di, 'di', supoc_bi, 'bi');
                end
            elseif all(perfil_P == corriente1(1).P)
                if (isempty(supoc_di) || all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, perfil_q, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob);
                elseif (~isempty(supoc_di) && ~all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, perfil_q, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob, supoc_di, 'di');
                elseif (isempty(supoc_di) || all(supoc_di == 0)) && (~isempty(supoc_bi) && ~all(supoc_bi == 0))   
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, perfil_q, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob, supoc_bi, 'bi');
                else
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, perfil_q, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob, supoc_di, 'di', supoc_bi, 'bi');
                end
            elseif all(perfil_q == 0)
                if (isempty(supoc_di) || all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, 0, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob);
                elseif (~isempty(supoc_di) && ~all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, 0, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob, supoc_di, 'di');
                elseif (isempty(supoc_di) || all(supoc_di == 0)) && (~isempty(supoc_bi) && ~all(supoc_bi == 0))   
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, 0, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob, supoc_bi, 'bi');
                else
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, 0, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob, supoc_di, 'di', supoc_bi, 'bi');
                end
            else
                if (isempty(supoc_di) || all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, perfil_q, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob);
                elseif (~isempty(supoc_di) && ~all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, perfil_q, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob, supoc_di, 'di');
                elseif (isempty(supoc_di) || all(supoc_di == 0)) && (~isempty(supoc_bi) && ~all(supoc_bi == 0))   
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, perfil_q, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob, supoc_bi, 'bi');
                else
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, perfil_q, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob, supoc_di, 'di', supoc_bi, 'bi');
                end
            end
        catch
            Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, perfil_q, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob, supoc_di, 'di', supoc_bi, 'bi');
        end
        Torre1.balanmasa();
        
        for ittreww = 1:N
            handles.perfil_v_l.Data{ittreww, 1} = num2str(Torre1.perfil_v(ittreww));
            handles.perfil_v_l.Data{ittreww, 2} = num2str(Torre1.perfil_l(ittreww));            
            handles.perfil_t.Data{ittreww, 1} = num2str(Torre1.perfil_t(ittreww));
        end
        for itterew = 1:length(alim{2}.mezcla.comp)
            handles.flujosproductocomp.Data{1, itterew} = num2str(Torre1.dsubi(itterew));
            handles.flujosproductocomp.Data{2, itterew} = num2str(Torre1.bsubi(itterew));
        end
    else
        errordlg('Por favor revise los datos faltantes en color rojo')
    end
    
elseif isa(Torre1, 'NaphtaliSandholm') && get(handles.newtonraphsondata, 'Value') == 3
    
    colorcondiciondob3 = get(handles.condiciondob3, 'ForegroundColor');
    colorcondiciondob4 = get(handles.condiciondob4, 'ForegroundColor');
    colorpanelespecif = get(handles.panelespecif, 'ForegroundColor');
    negro = [0,0,0];
    
    if all(colorcondiciondob3 == negro) && all(colorcondiciondob4 == negro) && all(colorpanelespecif == negro)
        try
            if all(perfil_P == corriente1(1).P) && all(perfil_q == 0)
                if (isempty(supoc_di) || all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, 0, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob);
                elseif (~isempty(supoc_di) && ~all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, 0, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob, supoc_di, 'di');
                elseif (isempty(supoc_di) || all(supoc_di == 0)) && (~isempty(supoc_bi) && ~all(supoc_bi == 0))   
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, 0, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob, supoc_bi, 'bi');
                else
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, 0, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob, supoc_di, 'di', supoc_bi, 'bi');
                end
            elseif all(perfil_P == corriente1(1).P)
                if (isempty(supoc_di) || all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, perfil_q, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob);
                elseif (~isempty(supoc_di) && ~all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, perfil_q, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob, supoc_di, 'di');
                elseif (isempty(supoc_di) || all(supoc_di == 0)) && (~isempty(supoc_bi) && ~all(supoc_bi == 0))   
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, perfil_q, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob, supoc_bi, 'bi');
                else
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, perfil_q, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob, supoc_di, 'di', supoc_bi, 'bi');
                end
            elseif all(perfil_q == 0)
                if (isempty(supoc_di) || all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, 0, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob);
                elseif (~isempty(supoc_di) && ~all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, 0, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob, supoc_di, 'di');
                elseif (isempty(supoc_di) || all(supoc_di == 0)) && (~isempty(supoc_bi) && ~all(supoc_bi == 0))   
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, 0, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob, supoc_bi, 'bi');
                else
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, 0, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob, supoc_di, 'di', supoc_bi, 'bi');
                end
            else
                if (isempty(supoc_di) || all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, perfil_q, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob);
                elseif (~isempty(supoc_di) && ~all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, perfil_q, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob, supoc_di, 'di');
                elseif (isempty(supoc_di) || all(supoc_di == 0)) && (~isempty(supoc_bi) && ~all(supoc_bi == 0))   
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, perfil_q, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob, supoc_bi, 'bi');
                else
                    Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, perfil_q, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob, supoc_di, 'di', supoc_bi, 'bi');
                end
            end
        catch
            Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, perfil_q, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob, supoc_di, 'di', supoc_bi, 'bi');
        end
        Torre1.balanmasa();
        
        for ittreww = 1:N
            handles.perfil_v_l.Data{ittreww, 1} = num2str(Torre1.perfil_v(ittreww));
            handles.perfil_v_l.Data{ittreww, 2} = num2str(Torre1.perfil_l(ittreww));            
            handles.perfil_t.Data{ittreww, 1} = num2str(Torre1.perfil_t(ittreww));
        end
        for itterew = 1:length(alim{2}.mezcla.comp)
            handles.flujosproductocomp.Data{1, itterew} = num2str(Torre1.dsubi(itterew));
            handles.flujosproductocomp.Data{2, itterew} = num2str(Torre1.bsubi(itterew));
        end
    else
        errordlg('Por favor revise los datos faltantes en color rojo')
    end
    
        
else
    errordlg('Seleccione un método de resolución riguroso ')
end



function iteraciones_Callback(hObject, eventdata, handles)
% hObject    handle to iteraciones (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global iteraciones_max

iteraciones_max = str2double(get(handles.iteraciones, 'String'));



% --- Executes during object creation, after setting all properties.
function iteraciones_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iteraciones (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calculartorre1.
function calculartorre1_Callback(hObject, eventdata, handles)
% hObject    handle to calculartorre1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global D; global dob; global salidas; global R; global corriente1; global Torre1;
global perfil_P; global perfil_q; global N; global alim; global iteraciones_max
global supoc_di; global supoc_bi; global otradob; global damping

if isempty(iteraciones_max)
    iteraciones_max = 100;
end

if isa(Torre1, 'BubblePoint')
    if isempty(Torre1.perfil_v)
        colorcondiciondbubble = get(handles.condiciondbubble, 'ForegroundColor');
        colorpanelespecif = get(handles.condiciondbubble, 'ForegroundColor');
        colorreflujotext = get(handles.reflujotext, 'ForegroundColor');
        negro = [0,0,0];

        if all(colorcondiciondbubble == negro) && all(colorpanelespecif == negro) && all(colorreflujotext == negro)
            if all(perfil_P == corriente1.P) && all(perfil_q == 0)
                Torre1 = BubblePoint(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1.P, 0, R, D, dob, iteraciones_max, 'T-001');

            elseif all(perfil_P == corriente1.P)

                Torre1 = BubblePoint(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1.P, perfil_q, R, D, dob, iteraciones_max, 'T-001');
            elseif all(perfil_q == 0)

                Torre1 = BubblePoint(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, 0, R, D, dob, iteraciones_max, 'T-001');
            else

                Torre1 = BubblePoint(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, perfil_q, R, D, dob, iteraciones_max, 'T-001');
            end

            Torre1.balanmasa();

            for ittreww = 1:N
                handles.perfil_v_l.Data{ittreww, 1} = num2str(Torre1.perfil_v(ittreww));
                handles.perfil_v_l.Data{ittreww, 2} = num2str(Torre1.perfil_l(ittreww));            
                handles.perfil_t.Data{ittreww, 1} = num2str(Torre1.perfil_t(ittreww));
            end
            for itterew = 1:length(alim{2}.mezcla.comp)
                handles.flujosproductocomp.Data{1, itterew} = num2str(Torre1.dsubi(itterew));
                handles.flujosproductocomp.Data{2, itterew} = num2str(Torre1.bsubi(itterew));
            end
            Torre1.wanghenke();
            for ittreww = 1:N
                handles.perfil_v_l.Data{ittreww, 3} = num2str(Torre1.perfil_v(ittreww));
                handles.perfil_v_l.Data{ittreww, 4} = num2str(Torre1.perfil_l(ittreww));            
                handles.perfil_t.Data{ittreww, 2} = num2str(Torre1.perfil_t(ittreww));
            end
            for itterew = 1:length(alim{2}.mezcla.comp)
                handles.flujosproductocomp.Data{1, itterew} = num2str(Torre1.dsubi(itterew));
                handles.flujosproductocomp.Data{2, itterew} = num2str(Torre1.bsubi(itterew));
            end                 
            handles.presionycalor.Data{1,2} = Torre1.qc;
            handles.presionycalor.Data{N,2} = Torre1.qb;
            set(handles.calorqc, 'String', num2str(Torre1.qc));
            set(handles.calorqb, 'String', num2str(Torre1.qb));
            set(handles.erroree, 'String', num2str(Torre1.actualvalI));
            set(handles.iteree, 'String', num2str(Torre1.actualiter));
            
        else
            errordlg('Por favor revise los datos faltantes en color rojo')
        end
    else
        Torre1.wanghenke();
        for ittreww = 1:N
            handles.perfil_v_l.Data{ittreww, 3} = num2str(Torre1.perfil_v(ittreww));
            handles.perfil_v_l.Data{ittreww, 4} = num2str(Torre1.perfil_l(ittreww));            
            handles.perfil_t.Data{ittreww, 2} = num2str(Torre1.perfil_t(ittreww));
        end
        for itterew = 1:length(alim{2}.mezcla.comp)
            handles.flujosproductocomp.Data{1, itterew} = num2str(Torre1.dsubi(itterew));
            handles.flujosproductocomp.Data{2, itterew} = num2str(Torre1.bsubi(itterew));
        end     
        handles.presionycalor.Data{1,2} = Torre1.qc;
        handles.presionycalor.Data{N,2} = Torre1.qb;
        set(handles.calorqc, 'String', num2str(Torre1.qc));
        set(handles.calorqb, 'String', num2str(Torre1.qb));
        set(handles.erroree, 'String', num2str(Torre1.actualvalI));
        set(handles.iteree, 'String', num2str(Torre1.actualiter));
    end
elseif isa(Torre1, 'NaphtaliSandholm') && get(handles.newtonraphsondata, 'Value') == 2
    if isempty(Torre1.perfil_v)
        colorcondiciondob2 = get(handles.condiciondob2, 'ForegroundColor');
        colorpanelespecif = get(handles.panelespecif, 'ForegroundColor');
        colorreflujotext = get(handles.reflujotext2, 'ForegroundColor');
        negro = [0,0,0];

        if all(colorcondiciondob2 == negro) && all(colorpanelespecif == negro) && all(colorreflujotext == negro)
            try
                if all(perfil_P == corriente1(1).P) && all(perfil_q == 0)
                    if (isempty(supoc_di) || all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, 0, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob);
                    elseif (~isempty(supoc_di) && ~all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, 0, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob, supoc_di, 'di');
                    elseif (isempty(supoc_di) || all(supoc_di == 0)) && (~isempty(supoc_bi) && ~all(supoc_bi == 0))   
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, 0, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob, supoc_bi, 'bi');
                    else
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, 0, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob, supoc_di, 'di', supoc_bi, 'bi');
                    end
                elseif all(perfil_P == corriente1(1).P)
                    if (isempty(supoc_di) || all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, perfil_q, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob);
                    elseif (~isempty(supoc_di) && ~all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, perfil_q, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob, supoc_di, 'di');
                    elseif (isempty(supoc_di) || all(supoc_di == 0)) && (~isempty(supoc_bi) && ~all(supoc_bi == 0))   
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, perfil_q, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob, supoc_bi, 'bi');
                    else
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, perfil_q, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob, supoc_di, 'di', supoc_bi, 'bi');
                    end
                elseif all(perfil_q == 0)
                    if (isempty(supoc_di) || all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, 0, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob);
                    elseif (~isempty(supoc_di) && ~all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, 0, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob, supoc_di, 'di');
                    elseif (isempty(supoc_di) || all(supoc_di == 0)) && (~isempty(supoc_bi) && ~all(supoc_bi == 0))   
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, 0, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob, supoc_bi, 'bi');
                    else
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, 0, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob, supoc_di, 'di', supoc_bi, 'bi');
                    end
                else
                    if (isempty(supoc_di) || all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, perfil_q, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob);
                    elseif (~isempty(supoc_di) && ~all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, perfil_q, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob, supoc_di, 'di');
                    elseif (isempty(supoc_di) || all(supoc_di == 0)) && (~isempty(supoc_bi) && ~all(supoc_bi == 0))   
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, perfil_q, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob, supoc_bi, 'bi');
                    else
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, perfil_q, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob, supoc_di, 'di', supoc_bi, 'bi');
                    end
                end
            catch
                Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, perfil_q, iteraciones_max, D, dob, [], [], damping, R, otradob, D, dob, supoc_di, 'di', supoc_bi, 'bi');
            end
            Torre1.balanmasa();

            for ittreww = 1:N
                handles.perfil_v_l.Data{ittreww, 1} = num2str(Torre1.perfil_v(ittreww));
                handles.perfil_v_l.Data{ittreww, 2} = num2str(Torre1.perfil_l(ittreww));            
                handles.perfil_t.Data{ittreww, 1} = num2str(Torre1.perfil_t(ittreww));
            end
            for itterew = 1:length(alim{2}.mezcla.comp)
                handles.flujosproductocomp.Data{1, itterew} = num2str(Torre1.dsubi(itterew));
                handles.flujosproductocomp.Data{2, itterew} = num2str(Torre1.bsubi(itterew));
            end
            Torre1.newtonraphson2();
            for ittreww = 1:N
                handles.perfil_v_l.Data{ittreww, 3} = num2str(Torre1.perfil_v(ittreww));
                handles.perfil_v_l.Data{ittreww, 4} = num2str(Torre1.perfil_l(ittreww));            
                handles.perfil_t.Data{ittreww, 2} = num2str(Torre1.perfil_t(ittreww));
            end
            for itterew = 1:length(alim{2}.mezcla.comp)
                handles.flujosproductocomp.Data{1, itterew} = num2str(Torre1.dsubi(itterew));
                handles.flujosproductocomp.Data{2, itterew} = num2str(Torre1.bsubi(itterew));
            end                 
            handles.presionycalor.Data{1,2} = Torre1.qc;
            handles.presionycalor.Data{N,2} = Torre1.qb;
            set(handles.calorqc, 'String', num2str(Torre1.qc));
            set(handles.calorqb, 'String', num2str(Torre1.qb));
            set(handles.erroree, 'String', num2str(Torre1.actualvalI));
            set(handles.iteree, 'String', num2str(Torre1.actualiter));
        else
            errordlg('Por favor revise los datos faltantes en color rojo')
        end
    else
        Torre1.newtonraphson2();
        for ittreww = 1:N
            handles.perfil_v_l.Data{ittreww, 3} = num2str(Torre1.perfil_v(ittreww));
            handles.perfil_v_l.Data{ittreww, 4} = num2str(Torre1.perfil_l(ittreww));            
            handles.perfil_t.Data{ittreww, 2} = num2str(Torre1.perfil_t(ittreww));
        end
        for itterew = 1:length(alim{2}.mezcla.comp)
            handles.flujosproductocomp.Data{1, itterew} = num2str(Torre1.dsubi(itterew));
            handles.flujosproductocomp.Data{2, itterew} = num2str(Torre1.bsubi(itterew));
        end                 
        handles.presionycalor.Data{1,2} = Torre1.qc;
        handles.presionycalor.Data{N,2} = Torre1.qb;
        set(handles.calorqc, 'String', num2str(Torre1.qc));
        set(handles.calorqb, 'String', num2str(Torre1.qb));
        set(handles.erroree, 'String', num2str(Torre1.actualvalI));
        set(handles.iteree, 'String', num2str(Torre1.actualiter));
    end
elseif isa(Torre1, 'NaphtaliSandholm') && get(handles.newtonraphsondata, 'Value') == 3
    if isempty(Torre1.perfil_v)
        colorcondiciondob3 = get(handles.condiciondob3, 'ForegroundColor');
        colorcondiciondob4 = get(handles.condiciondob4, 'ForegroundColor');
        colorpanelespecif = get(handles.panelespecif, 'ForegroundColor');
        negro = [0,0,0];

        if all(colorcondiciondob3 == negro) && all(colorcondiciondob4 == negro) && all(colorpanelespecif == negro)
            try
                if all(perfil_P == corriente1(1).P) && all(perfil_q == 0)
                    if (isempty(supoc_di) || all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, 0, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob);
                    elseif (~isempty(supoc_di) && ~all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, 0, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob, supoc_di, 'di');
                    elseif (isempty(supoc_di) || all(supoc_di == 0)) && (~isempty(supoc_bi) && ~all(supoc_bi == 0))   
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, 0, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob, supoc_bi, 'bi');
                    else
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, 0, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob, supoc_di, 'di', supoc_bi, 'bi');
                    end
                elseif all(perfil_P == corriente1(1).P)
                    if (isempty(supoc_di) || all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, perfil_q, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob);
                    elseif (~isempty(supoc_di) && ~all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, perfil_q, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob, supoc_di, 'di');
                    elseif (isempty(supoc_di) || all(supoc_di == 0)) && (~isempty(supoc_bi) && ~all(supoc_bi == 0))   
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, perfil_q, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob, supoc_bi, 'bi');
                    else
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1(1).P, perfil_q, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob, supoc_di, 'di', supoc_bi, 'bi');
                    end
                elseif all(perfil_q == 0)
                    if (isempty(supoc_di) || all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, 0, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob);
                    elseif (~isempty(supoc_di) && ~all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, 0, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob, supoc_di, 'di');
                    elseif (isempty(supoc_di) || all(supoc_di == 0)) && (~isempty(supoc_bi) && ~all(supoc_bi == 0))   
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, 0, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob, supoc_bi, 'bi');
                    else
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, 0, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob, supoc_di, 'di', supoc_bi, 'bi');
                    end
                else
                    if (isempty(supoc_di) || all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, perfil_q, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob);
                    elseif (~isempty(supoc_di) && ~all(supoc_di == 0)) && (isempty(supoc_bi) || all(supoc_bi == 0)) 
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, perfil_q, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob, supoc_di, 'di');
                    elseif (isempty(supoc_di) || all(supoc_di == 0)) && (~isempty(supoc_bi) && ~all(supoc_bi == 0))   
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, perfil_q, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob, supoc_bi, 'bi');
                    else
                        Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, perfil_q, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob, supoc_di, 'di', supoc_bi, 'bi');
                    end
                end
            catch
                Torre1 = NaphtaliSandholm(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, perfil_q, iteraciones_max, [], [], [], [], damping, R, otradob, D, dob, supoc_di, 'di', supoc_bi, 'bi');
            end
            Torre1.balanmasa();

            for ittreww = 1:N
                handles.perfil_v_l.Data{ittreww, 1} = num2str(Torre1.perfil_v(ittreww));
                handles.perfil_v_l.Data{ittreww, 2} = num2str(Torre1.perfil_l(ittreww));            
                handles.perfil_t.Data{ittreww, 1} = num2str(Torre1.perfil_t(ittreww));
            end
            for itterew = 1:length(alim{2}.mezcla.comp)
                handles.flujosproductocomp.Data{1, itterew} = num2str(Torre1.dsubi(itterew));
                handles.flujosproductocomp.Data{2, itterew} = num2str(Torre1.bsubi(itterew));
            end
            Torre1.newtonraphson3();
            for ittreww = 1:N
                handles.perfil_v_l.Data{ittreww, 3} = num2str(Torre1.perfil_v(ittreww));
                handles.perfil_v_l.Data{ittreww, 4} = num2str(Torre1.perfil_l(ittreww));            
                handles.perfil_t.Data{ittreww, 2} = num2str(Torre1.perfil_t(ittreww));
            end
            for itterew = 1:length(alim{2}.mezcla.comp)
                handles.flujosproductocomp.Data{1, itterew} = num2str(Torre1.dsubi(itterew));
                handles.flujosproductocomp.Data{2, itterew} = num2str(Torre1.bsubi(itterew));
            end                 
            handles.presionycalor.Data{1,2} = Torre1.qc;
            handles.presionycalor.Data{N,2} = Torre1.qb;
            set(handles.calorqc, 'String', num2str(Torre1.qc));
            set(handles.calorqb, 'String', num2str(Torre1.qb));
            set(handles.erroree, 'String', num2str(Torre1.actualvalI));
            set(handles.iteree, 'String', num2str(Torre1.actualiter));
        else
            errordlg('Por favor revise los datos faltantes en color rojo')
        end
    else
         Torre1.newtonraphson3();
        for ittreww = 1:N
            handles.perfil_v_l.Data{ittreww, 3} = num2str(Torre1.perfil_v(ittreww));
            handles.perfil_v_l.Data{ittreww, 4} = num2str(Torre1.perfil_l(ittreww));            
            handles.perfil_t.Data{ittreww, 2} = num2str(Torre1.perfil_t(ittreww));
        end
        for itterew = 1:length(alim{2}.mezcla.comp)
            handles.flujosproductocomp.Data{1, itterew} = num2str(Torre1.dsubi(itterew));
            handles.flujosproductocomp.Data{2, itterew} = num2str(Torre1.bsubi(itterew));
        end                 
        handles.presionycalor.Data{1,2} = Torre1.qc;
        handles.presionycalor.Data{N,2} = Torre1.qb;
        set(handles.calorqc, 'String', num2str(Torre1.qc));
        set(handles.calorqb, 'String', num2str(Torre1.qb));
        set(handles.erroree, 'String', num2str(Torre1.actualvalI));
        set(handles.iteree, 'String', num2str(Torre1.actualiter));
    end
else
    errordlg('Seleccione un método de resolución riguroso ')
end


% --- Executes when entered data in editable cell(s) in perfil_v_l.
function perfil_v_l_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to perfil_v_l (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)



function calorqc_Callback(hObject, eventdata, handles)
% hObject    handle to calorqc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of calorqc as text
%        str2double(get(hObject,'String')) returns contents of calorqc as a double


% --- Executes during object creation, after setting all properties.
function calorqc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calorqc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function calorqb_Callback(hObject, eventdata, handles)
% hObject    handle to calorqb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of calorqb as text
%        str2double(get(hObject,'String')) returns contents of calorqb as a double


% --- Executes during object creation, after setting all properties.
function calorqb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calorqb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function erroree_Callback(hObject, eventdata, handles)
% hObject    handle to erroree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of erroree as text
%        str2double(get(hObject,'String')) returns contents of erroree as a double


% --- Executes during object creation, after setting all properties.
function erroree_CreateFcn(hObject, eventdata, handles)
% hObject    handle to erroree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function iteree_Callback(hObject, eventdata, handles)
% hObject    handle to iteree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iteree as text
%        str2double(get(hObject,'String')) returns contents of iteree as a double


% --- Executes during object creation, after setting all properties.
function iteree_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iteree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pasoapaso.
function pasoapaso_Callback(hObject, eventdata, handles)
% hObject    handle to pasoapaso (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Torre1; global salidas; global damping; global prev_pc; 
global prev_tc; global prev_wa; global condicion1; global condicion2;
global MulticompEdE; global ident_sust; global formula; global masamolar;
global Tc; global Pc; global wa; global valor; global Sust; global kij;
global T_o_P_o_beta; global P_o_beta_o_q; global F; global corriente1;
global mezcla1; global numer_corrientes; global alim; global etapasalim;
global conci; global Fi; global perfil_P; global perfil_q; global N;
global perfil_v; global perfil_vi; global perfil_l; global perfil_li; 
global xsubj; global ysubj; global perfil_k; global bsubi; global dsubi; 
global R; global dob; global supoc_di; global supoc_bi; global otradob;

alim = cell.empty(0,1);
damping = 0.6;
condicion1 = double.empty(0,1);
condicion2 = double.empty(0,1);
perfil_q = double.empty(0,1);
perfil_P = double.empty(0,1);
Torre1 = double.empty(0,1);
conci = double.empty(0,1);
mezcla1 = Mezcla.empty(0,1);
P_o_beta_o_q = double.empty(0,1);
T_o_P_o_beta = double.empty(0,1);
R = double.empty(0,0);
dob = char.empty(0);
otradob = char.empty(0);
Tc = double.empty(0,1);
Pc = double.empty(0,1);
corriente1 = Corriente.empty(0,1);
MulticompEdE = RMVdW.empty(0,1);
wa = double.empty(0,1);
numer_corrientes = double.empty(0,1);
etapasalim = double.empty(0,1);
supoc_di = double.empty(0,1);
supoc_bi = double.empty(0,1);
bsubi = double.empty(0,1);
dsubi = double.empty(0,1);
ysubj = double.empty(0,1);
xsubj = double.empty(0,1);
Fi = double.empty(0,1);
F = double.empty(0,1);
kij = double.empty(0,1);
Sust = Sustancia.empty(0,1);
ident_sust = char.empty(0);
salidas = cell.empty(0,1);
formula = char.empty(0);
masamolar = double.empty(0,1);
valor = double.empty(0,1);
prev_tc = double.empty(0,1);
prev_pc = double.empty(0,1);
prev_wa = double.empty(0,1);
N = double.empty(0,1);
perfil_v = double.empty(0,1);
perfil_vi = double.empty(0,1);
perfil_l = double.empty(0,1);
perfil_k = double.empty(0,1);
perfil_li = double.empty(0,1);





% --- Executes on button press in avanza1paso.
function avanza1paso_Callback(hObject, eventdata, handles)
% hObject    handle to avanza1paso (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global D; global dob; global salidas; global R; global corriente1; global Torre1;
global perfil_P; global perfil_q; global N; global alim; global iteraciones_max

if isempty(iteraciones_max)
    iteraciones_max = 100;
end

if isa(Torre1, 'BubblePoint')
    if isempty(Torre1.perfil_v)
        colorcondiciondbubble = get(handles.condiciondbubble, 'ForegroundColor');
        colorpanelespecif = get(handles.condiciondbubble, 'ForegroundColor');
        colorreflujotext = get(handles.reflujotext, 'ForegroundColor');
        negro = [0,0,0];

        if all(colorcondiciondbubble == negro) && all(colorpanelespecif == negro) && all(colorreflujotext == negro)
            if all(perfil_P == corriente1.P) && all(perfil_q == 0)
                Torre1 = BubblePoint(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1.P, 0, R, D, dob, iteraciones_max, 'T-001');

            elseif all(perfil_P == corriente1.P)

                Torre1 = BubblePoint(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, corriente1.P, perfil_q, R, D, dob, iteraciones_max, 'T-001');
            elseif all(perfil_q == 0)

                Torre1 = BubblePoint(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, 0, R, D, dob, iteraciones_max, 'T-001');
            else

                Torre1 = BubblePoint(N, alim, Corriente.empty(0,1), Corriente.empty(0,1), salidas, perfil_P, perfil_q, R, D, dob, iteraciones_max, 'T-001');
            end

            Torre1.balanmasa();

            for ittreww = 1:N
                handles.perfil_v_l.Data{ittreww, 1} = num2str(Torre1.perfil_v(ittreww));
                handles.perfil_v_l.Data{ittreww, 2} = num2str(Torre1.perfil_l(ittreww));            
                handles.perfil_t.Data{ittreww, 1} = num2str(Torre1.perfil_t(ittreww));
            end
            for itterew = 1:length(alim{2}.mezcla.comp)
                handles.flujosproductocomp.Data{1, itterew} = num2str(Torre1.dsubi(itterew));
                handles.flujosproductocomp.Data{2, itterew} = num2str(Torre1.bsubi(itterew));
            end
            Torre1.avanza1paso();
             
            for ittreww = 1:N
                handles.perfil_v_l.Data{ittreww, 3} = num2str(Torre1.perfil_v(ittreww));
                handles.perfil_v_l.Data{ittreww, 4} = num2str(Torre1.perfil_l(ittreww));            
                handles.perfil_t.Data{ittreww, 2} = num2str(Torre1.perfil_t(ittreww));
            end
            for itterew = 1:length(alim{2}.mezcla.comp)
                handles.flujosproductocomp.Data{1, itterew} = num2str(Torre1.dsubi(itterew));
                handles.flujosproductocomp.Data{2, itterew} = num2str(Torre1.bsubi(itterew));
            end                 
            handles.presionycalor.Data{1,2} = Torre1.qc;
            handles.presionycalor.Data{N,2} = Torre1.qb;
            set(handles.calorqc, 'String', num2str(Torre1.qc));
            set(handles.calorqb, 'String', num2str(Torre1.qb));
            set(handles.erroree, 'String', num2str(Torre1.actualvalI));
            set(handles.iteree, 'String', num2str(Torre1.actualiter));
            
        else
            errordlg('Por favor revise los datos faltantes en color rojo')
        end
    else
        Torre1.avanza1paso();
        for ittreww = 1:N
            handles.perfil_v_l.Data{ittreww, 3} = num2str(Torre1.perfil_v(ittreww));
            handles.perfil_v_l.Data{ittreww, 4} = num2str(Torre1.perfil_l(ittreww));            
            handles.perfil_t.Data{ittreww, 2} = num2str(Torre1.perfil_t(ittreww));
        end
        for itterew = 1:length(alim{2}.mezcla.comp)
            handles.flujosproductocomp.Data{1, itterew} = num2str(Torre1.dsubi(itterew));
            handles.flujosproductocomp.Data{2, itterew} = num2str(Torre1.bsubi(itterew));
        end     
        handles.presionycalor.Data{1,2} = Torre1.qc;
        handles.presionycalor.Data{N,2} = Torre1.qb;
        set(handles.calorqc, 'String', num2str(Torre1.qc));
        set(handles.calorqb, 'String', num2str(Torre1.qb));
        set(handles.erroree, 'String', num2str(Torre1.actualvalI));
        set(handles.iteree, 'String', num2str(Torre1.actualiter));
    end
elseif isa(Torre1, 'NaphtaliSandholm')
else
    errordlg('Seleccione un método de resolución riguroso ')
end

function condiciondob2_Callback(hObject, eventdata, handles)
% hObject    handle to condiciondob2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns condiciondbubble contents as cell array
%        contents{get(hObject,'Value')} returns selected item from condiciondbubble
global salidas
if get(handles.condiciondob2, 'Value') == 2
    if get(handles.destiladofase, 'Value') == 2 || get(handles.destiladofase, 'Value') == 3

        set(handles.liquidodestilado, 'Visible', 'on')

    end
else
      set(handles.liquidodestilado, 'Visible', 'off')
end

% --- Executes during object creation, after setting all properties.
function condiciondob2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to condiciondob2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function flujod_o_b2_Callback(hObject, eventdata, handles)
% hObject    handle to flujod_o_b2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of flujod_o_b as text
%        str2double(get(hObject,'String')) returns contents of flujod_o_b as a double
global D; global salidas; global dob

tipo_de_destilado = get(handles.destiladofase, 'Value'); 
if tipo_de_destilado == 1
    D = str2double(get(handles.flujod_o_b2, 'String'));
    esdob = get(handles.condiciondob2, 'Value');
    if esdob == 1
        dob = 'D';
    else
        dob = 'B';
    end
    set(handles.condiciondob2,'ForegroundColor', [0,0,0]);
    colorreflujotext2 = get(handles.reflujotext2, 'ForegroundColor');
    if all(colorreflujotext2 == [0,0,0])
        set(handles.panelespecif, 'ForegroundColor', [0,0,0]);
    end
    colorpanelespecif = get(handles.condiciondob2, 'ForegroundColor');
    colorreflujotext = get(handles.reflujotext2, 'ForegroundColor');
    negro = [0,0,0];
    if all(colorreflujotext == negro) && all(colorpanelespecif == negro)
        set(handles.newtonraphsondata, 'ForegroundColor', [0,0,0])
    end
elseif tipo_de_destilado == 2
    D = str2double(get(handles.flujod_o_b2, 'String'));
    esdob = get(handles.condiciondob2, 'Value');
    if esdob == 1
        dob = 'D';
    else
        dob = 'B';
    end
    if esdob == 1
        salidas = {1,1,D};
        colorreflujotext2 = get(handles.panelespecif, 'ForegroundColor');
        if all(colorreflujotext2 == [0,0,0])
            set(handles.panelespecif, 'ForegroundColor', [0,0,0]);
        end
    end
    set(handles.condiciondob2, 'ForegroundColor',[0,0,0])
    colorpanelespecif = get(handles.condiciondob2, 'ForegroundColor');
    colorreflujotext = get(handles.reflujotext2, 'ForegroundColor');
    negro = [0,0,0];
    if all(colorreflujotext == negro) && all(colorpanelespecif == negro)
        set(handles.newtonraphsondata, 'ForegroundColor', [0,0,0])
    end
else
    set(handles.liquidodestilado, 'Visible', 'on');
    D = str2double(get(handles.flujod_o_b2, 'String'));
    esdob = get(handles.condiciondob2, 'Value');
    if esdob == 1
        dob = 'D';
    else
        dob = 'B';
    end
    set(handles.condiciondob2, 'ForegroundColor',[0,0,0])
    colorpanelespecif = get(handles.condiciondob2, 'ForegroundColor');
    colorreflujotext = get(handles.reflujotext2, 'ForegroundColor');
    negro = [0,0,0];
    if all(colorreflujotext == negro) && all(colorpanelespecif == negro)
        set(handles.newtonraphsondata, 'ForegroundColor', [0,0,0])
    end
end


% --- Executes during object creation, after setting all properties.
function flujod_o_b2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to flujod_o_b2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function flujod_o_b3_Callback(hObject, eventdata, handles)
% hObject    handle to flujod_o_b3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of flujod_o_b as text
%        str2double(get(hObject,'String')) returns contents of flujod_o_b as a double
global D; global salidas; global dob;

    if ~isa(D, 'cell') || isempty(D)
        D = cell.empty(1,0);
        D{1} = str2double(get(handles.flujod_o_b3, 'String'));
    else
        D{1} = str2double(get(handles.flujod_o_b3, 'String'));
    end
    dob = 'di';

    set(handles.condiciondob3,'ForegroundColor', [0,0,0]);
    set(handles.panelespecif, 'ForegroundColor', [0,0,0]);
    colorpanelespecif = get(handles.condiciondob3, 'ForegroundColor');
    colorreflujotext = get(handles.condiciondob4, 'ForegroundColor');
    negro = [0,0,0];
    if all(colorreflujotext == negro) && all(colorpanelespecif == negro)
        set(handles.newtonraphsondata, 'ForegroundColor', [0,0,0])
    end

% --- Executes during object creation, after setting all properties.
function flujod_o_b3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to flujod_o_b3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function flujod_o_b4_Callback(hObject, eventdata, handles)
% hObject    handle to flujod_o_b3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of flujod_o_b as text
%        str2double(get(hObject,'String')) returns contents of flujod_o_b as a double
global D; global salidas;

    if ~isa(D, 'cell') || isempty(D)
        D = cell.empty(1,0);
        D{2} = str2double(get(handles.flujod_o_b4, 'String'));
    else
        D{2} = str2double(get(handles.flujod_o_b4, 'String'));
    end
    
    set(handles.condiciondob3,'ForegroundColor', [0,0,0]);
    set(handles.panelespecif, 'ForegroundColor', [0,0,0]);
    colorpanelespecif = get(handles.condiciondob3, 'ForegroundColor');
    colorreflujotext = get(handles.condiciondob4, 'ForegroundColor');
    negro = [0,0,0];
    if all(colorreflujotext == negro) && all(colorpanelespecif == negro)
        set(handles.newtonraphsondata, 'ForegroundColor', [0,0,0])
    end

% --- Executes during object creation, after setting all properties.
function flujod_o_b4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to flujod_o_b3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function flujod_o_b5_Callback(hObject, eventdata, handles)
% hObject    handle to flujod_o_b3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of flujod_o_b as text
%        str2double(get(hObject,'String')) returns contents of flujod_o_b as a double
global R; global salidas; global dob; global otradob

    if ~isa(R, 'cell') || isempty(R)
        R = cell.empty(1,0);
        R{1} = str2double(get(handles.flujod_o_b5, 'String'));
    else
        R{1} = str2double(get(handles.flujod_o_b5, 'String'));
    end
    otradob = 'bi';
    
    set(handles.condiciondob4,'ForegroundColor', [0,0,0]);
    set(handles.panelespecif, 'ForegroundColor', [0,0,0]);
    colorpanelespecif = get(handles.condiciondob3, 'ForegroundColor');
    colorreflujotext = get(handles.condiciondob4, 'ForegroundColor');
    negro = [0,0,0];
    if all(colorreflujotext == negro) && all(colorpanelespecif == negro)
        set(handles.newtonraphsondata, 'ForegroundColor', [0,0,0])
    end

% --- Executes during object creation, after setting all properties.
function flujod_o_b5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to flujod_o_b3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function reflujons_Callback(hObject, eventdata, handles)
% hObject    handle to reflujoBP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of reflujoBP as text
%        str2double(get(hObject,'String')) returns contents of reflujoBP as a double
global R; global otradob;

reflujo = str2double(get(handles.reflujons,'String'));
if ~isempty(reflujo) && isa(reflujo,'double') && (reflujo > 0)
    set(handles.reflujotext2, 'ForegroundColor', [0,0,0])
    colorcondiciondob2 = get(handles.condiciondob2, 'ForegroundColor');
    if all(colorcondiciondob2 == [0,0,0])
        set(handles.panelespecif, 'ForegroundColor', [0,0,0])
    end
    R = reflujo;
    otradob = 'R';
    colorcondiciondob2 = get(handles.condiciondob2, 'ForegroundColor');
    colorpanelespecif = get(handles.panelespecif, 'ForegroundColor');
    negro = [0,0,0];
    if all(colorcondiciondob2 == negro) && all(colorpanelespecif == negro)
        set(handles.newtonraphsondata, 'ForegroundColor', [0,0,0])
    end
else 
    set(handles.reflujotext2, 'ForegroundColor', [1,0,0]);
    set(handles.panelespecif, 'ForegroundColor', [1,0,0]);    
    errordlg('Reflujo debe valor numérico ser mayor a 0');
end

% --- Executes during object creation, after setting all properties.
function reflujons_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reflujoBP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function flujod_o_b6_Callback(hObject, eventdata, handles)
% hObject    handle to flujod_o_b6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of flujod_o_b as text
%        str2double(get(hObject,'String')) returns contents of flujod_o_b as a double
global R; global salidas; global dob; global otradob

    if ~isa(R, 'cell') || isempty(R)
        R = cell.empty(1,0);
        R{2} = str2double(get(handles.flujod_o_b6, 'String'));
    else
        R{2} = str2double(get(handles.flujod_o_b6, 'String'));
    end
    otradob = 'bi';
    
    set(handles.condiciondob4,'ForegroundColor', [0,0,0]);
    set(handles.panelespecif, 'ForegroundColor', [0,0,0]);
    colorpanelespecif = get(handles.condiciondob3, 'ForegroundColor');
    colorreflujotext = get(handles.condiciondob4, 'ForegroundColor');
    negro = [0,0,0];
    if all(colorreflujotext == negro) && all(colorpanelespecif == negro)
        set(handles.newtonraphsondata, 'ForegroundColor', [0,0,0])
    end

% --- Executes during object creation, after setting all properties.
function flujod_o_b6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to flujod_o_b6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function flujod_o_b7_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function flujod_o_b7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to flujod_o_b3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function flujod_o_b8_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function flujod_o_b8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to flujod_o_b8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function flujosproductocomp_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
global supoc_di; global supoc_bi; global corriente1


if ~isempty(corriente1)
    indices = eventdata.Indices;
    if indices(1) == 1
        if isempty(supoc_di) || any(isnan(supoc_bi) == 1)
            supoc_di = zeros(1, corriente1(1).num_sust);
        end
        supoc_di(indices(2)) =  eventdata.NewData;        
    elseif indices(1) == 2  || any(isnan(supoc_bi) == 1)
        if isempty(supoc_bi)
            supoc_bi = zeros(1, corriente1(1).num_sust);
        end
        supoc_bi(indices(2)) = eventdata.NewData;      
    end
    
end


