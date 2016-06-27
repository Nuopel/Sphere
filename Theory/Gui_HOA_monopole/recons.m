function varargout = recons(varargin)
% RECONS MATLAB code for recons.fig
%      RECONS, by itself, creates a new RECONS or raises the existing
%      singleton*.
%
%      H = RECONS returns the handle to a new RECONS or the handle to
%      the existing singleton*.
%
%      RECONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RECONS.M with the given input arguments.
%
%      RECONS('Property','Value',...) creates a new RECONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before recons_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to recons_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help recons

% Last Modified by GUIDE v2.5 25-Feb-2016 09:16:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @recons_OpeningFcn, ...
                   'gui_OutputFcn',  @recons_OutputFcn, ...
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


% --- Executes just before recons is made visible.
function recons_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to recons (see VARARGIN)

% Choose default command line output for recons
handles.output = hObject;

handles.N = 1001;

handles.M = 1;
handles.rs = 5;
handles.phi = 1;
handles.theta= 1;
handles.f = 500;

handles.r_hp = 1;
handles.N_hp = 50;

handles.Mslider= linspace(0,25,handles.N).*handles.M;
handles.rsslider= linspace(0,5,handles.N).*handles.rs;
handles.phislider= linspace(0,180,handles.N).*handles.phi;
handles.thetaslider= linspace(0,360,handles.N).*handles.theta;
handles.fslider= logspace(log10(0.1),1,handles.N).*handles.f;

handles.r_hpslider= logspace(log10(0.1),1,handles.N).*handles.r_hp;

[m p] = min(abs(handles.Mslider-handles.M));
set(handles.M_slider, 'value', (p-1)/(length(handles.Mslider)-1));
set(handles.M_text, 'string', ['M order = ' num2str(round(handles.M)) ]);

[m p] = min(abs(handles.rsslider-handles.M));
set(handles.rs_slider, 'value', (p-1)/(length(handles.rsslider)-1));
set(handles.rs_text, 'string', ['Rs = ' num2str(round(handles.rs)) ]);

[m p] = min(abs(handles.phislider-handles.phi));
set(handles.phi_slider, 'value', (p-1)/(length(handles.phislider)-1));
set(handles.phi_text, 'string', ['Phi = ' num2str(round(handles.phi)) ]);

[m p] = min(abs(handles.thetaslider-handles.theta));
set(handles.theta_slider, 'value', (p-1)/(length(handles.thetaslider)-1));
set(handles.theta_text, 'string', ['Theta = ' num2str(round(handles.theta)) ]);

[m p] = min(abs(handles.fslider-handles.f));
set(handles.f_slider, 'value', (p-1)/(length(handles.fslider)-1));
set(handles.f_text, 'string', ['Frequency = ' num2str(round(handles.f)) ]);

[m p] = min(abs(handles.r_hpslider-handles.r_hp));
set(handles.r_hp_slider, 'value', (p-1)/(length(handles.r_hpslider)-1));
set(handles.r_hp_text, 'string', ['R hp = ' num2str(round(handles.r_hp)) ]);

% [m p] = min(abs(handles.N_hpslider-handles.N_hp));
% set(handles.N_hp_list, 'value', (p-1)/(length(handles.N_hpslider)-1));
% set(handles.N_hp_text, 'string', ['Nbr hp = ' num2str(round(handles.N_hp)) ]);

[handles.x1_vec,handles.Pmes_mat,handles.Precons_mat2,erreur_mat]=call_HOA(handles.M,handles.rs,handles.theta,handles.phi,handles.f,handles.r_hp,handles.N_hp);

pcolor(handles.Pmes_graph,handles.x1_vec,handles.x1_vec,real(handles.Pmes_mat));
cax = caxis(handles.Pmes_graph);
shading(handles.Pmes_graph, 'interp');
colorbar(handles.Pmes_graph,'location','EastOutside');

% shading interp
xlabel('Position x (m)')
ylabel('Position y (m)')

% hcb = colorbar('location','EastOutside');

pcolor(handles.Precons_graph,handles.x1_vec,handles.x1_vec,real(handles.Precons_mat2));
colorbar(handles.Precons_graph,'location','EastOutside');
shading(handles.Precons_graph, 'interp')
xlabel('Position x (m)')
ylabel('Position y (m)')
hold(handles.Precons_graph,'on')
contour(handles.Precons_graph,handles.x1_vec,handles.x1_vec,erreur_mat,[0 10])
caxis(handles.Precons_graph,cax)
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes recons wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = recons_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function M_slider_Callback(hObject, eventdata, handles)
% hObject    handle to M_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
D = get(hObject,'value');
p = round(D*(handles.N-1)+1);
handles.M = round(handles.Mslider(p));
str=sprintf('M order = %i ',handles.M);
set(handles.M_text, 'string', str);

[handles.x1_vec,handles.Pmes_mat,handles.Precons_mat2,erreur_mat]=call_HOA(handles.M,handles.rs,handles.theta,handles.phi,handles.f,handles.r_hp,handles.N_hp);


pcolor(handles.Pmes_graph,handles.x1_vec,handles.x1_vec,real(handles.Pmes_mat));
cax = caxis(handles.Pmes_graph);
shading(handles.Pmes_graph, 'interp');
colorbar(handles.Pmes_graph,'location','EastOutside');

% shading interp
xlabel('Position x (m)')
ylabel('Position y (m)')


pcolor(handles.Precons_graph,handles.x1_vec,handles.x1_vec,real(handles.Precons_mat2));
colorbar(handles.Precons_graph,'location','EastOutside');
shading(handles.Precons_graph, 'interp')
xlabel('Position x (m)')
ylabel('Position y (m)')
hold(handles.Precons_graph,'on')
contour(handles.Precons_graph,handles.x1_vec,handles.x1_vec,erreur_mat,[0 10])
caxis(handles.Precons_graph,cax)

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function M_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to M_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in N_hp_list.
function N_hp_list_Callback(hObject, eventdata, handles)
% hObject    handle to N_hp_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns N_hp_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from N_hp_list
contents=cellstr(get(hObject,'String'));
choix=cell2mat({get(hObject,'Value')});

handles.N_hp=str2double(contents(choix))
[handles.x1_vec,handles.Pmes_mat,handles.Precons_mat2,erreur_mat]=call_HOA(handles.M,handles.rs,handles.theta,handles.phi,handles.f,handles.r_hp,handles.N_hp);

pcolor(handles.Pmes_graph,handles.x1_vec,handles.x1_vec,real(handles.Pmes_mat));
shading(handles.Pmes_graph, 'interp')
colorbar(handles.Pmes_graph,'location','EastOutside');
cax = caxis(handles.Pmes_graph);
xlabel('Position x (m)')
ylabel('Position y (m)')
% hcb = colorbar('location','EastOutside');
pcolor(handles.Precons_graph,handles.x1_vec,handles.x1_vec,real(handles.Precons_mat2));
colorbar(handles.Precons_graph,'location','EastOutside');
shading(handles.Precons_graph, 'interp')
xlabel('Position x (m)')
ylabel('Position y (m)')
hold(handles.Precons_graph,'on')
contour(handles.Precons_graph,handles.x1_vec,handles.x1_vec,erreur_mat,[0 10])
caxis(handles.Precons_graph,cax)

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function N_hp_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N_hp_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function f_slider_Callback(hObject, eventdata, handles)
% hObject    handle to f_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
D = get(hObject,'value');
p = round(D*(handles.N-1)+1);
handles.f = round(handles.fslider(p));
str=sprintf('Frequency = %i ',handles.f);
set(handles.f_text, 'string', str);
[handles.x1_vec,handles.Pmes_mat,handles.Precons_mat2,erreur_mat]=call_HOA(handles.M,handles.rs,handles.theta,handles.phi,handles.f,handles.r_hp,handles.N_hp);

pcolor(handles.Pmes_graph,handles.x1_vec,handles.x1_vec,real(handles.Pmes_mat));
shading(handles.Pmes_graph, 'interp')
colorbar(handles.Pmes_graph,'location','EastOutside');
cax = caxis(handles.Pmes_graph);
xlabel('Position x (m)')
ylabel('Position y (m)')
% hcb = colorbar('location','EastOutside');
pcolor(handles.Precons_graph,handles.x1_vec,handles.x1_vec,real(handles.Precons_mat2));
colorbar(handles.Precons_graph,'location','EastOutside');
shading(handles.Precons_graph, 'interp')
xlabel('Position x (m)')
ylabel('Position y (m)')
hold(handles.Precons_graph,'on')
contour(handles.Precons_graph,handles.x1_vec,handles.x1_vec,erreur_mat,[0 10])
caxis(handles.Precons_graph,cax)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function f_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function rs_slider_Callback(hObject, eventdata, handles)
% hObject    handle to rs_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
D = get(hObject,'value');
p = round(D*(handles.N-1)+1);
handles.rs = round(handles.rsslider(p)*100)/100;
str=sprintf('R source = %i ',handles.rs);
set(handles.rs_text, 'string', str);
[handles.x1_vec,handles.Pmes_mat,handles.Precons_mat2,erreur_mat]=call_HOA(handles.M,handles.rs,handles.theta,handles.phi,handles.f,handles.r_hp,handles.N_hp);
cax = caxis(handles.Pmes_graph);
pcolor(handles.Pmes_graph,handles.x1_vec,handles.x1_vec,real(handles.Pmes_mat));
colorbar(handles.Pmes_graph,'location','EastOutside');
shading(handles.Pmes_graph, 'interp')
xlabel('Position x (m)')
ylabel('Position y (m)')
% hcb = colorbar('location','EastOutside');
pcolor(handles.Precons_graph,handles.x1_vec,handles.x1_vec,real(handles.Precons_mat2));
colorbar(handles.Precons_graph,'location','EastOutside');
shading(handles.Precons_graph, 'interp')
xlabel('Position x (m)')
ylabel('Position y (m)')
hold(handles.Precons_graph,'on')
contour(handles.Precons_graph,handles.x1_vec,handles.x1_vec,erreur_mat,[0 10])
caxis(handles.Precons_graph,cax)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function rs_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rs_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function r_hp_slider_Callback(hObject, eventdata, handles)
% hObject    handle to r_hp_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
D = get(hObject,'value');
p = round(D*(handles.N-1)+1);
handles.r_hp = round(handles.r_hpslider(p));
str=sprintf('R hp = %i ',handles.r_hp);
set(handles.r_hp_text, 'string', str);

[handles.x1_vec,handles.Pmes_mat,handles.Precons_mat2,erreur_mat]=call_HOA(handles.M,handles.rs,handles.theta,handles.phi,handles.f,handles.r_hp,handles.N_hp);

pcolor(handles.Pmes_graph,handles.x1_vec,handles.x1_vec,real(handles.Pmes_mat));
shading(handles.Pmes_graph, 'interp')
cax = caxis(handles.Pmes_graph);
colorbar(handles.Pmes_graph,'location','EastOutside');

xlabel('Position x (m)')
ylabel('Position y (m)')
% hcb = colorbar('location','EastOutside');
pcolor(handles.Precons_graph,handles.x1_vec,handles.x1_vec,real(handles.Precons_mat2));
colorbar(handles.Precons_graph,'location','EastOutside');
shading(handles.Precons_graph, 'interp')
xlabel('Position x (m)')
ylabel('Position y (m)')
hold(handles.Precons_graph,'on')
contour(handles.Precons_graph,handles.x1_vec,handles.x1_vec,erreur_mat,[0 10])
caxis(handles.Precons_graph,cax)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function r_hp_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r_hp_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function theta_slider_Callback(hObject, eventdata, handles)
% hObject    handle to theta_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
D = get(hObject,'value');
p = round(D*(handles.N-1)+1);
handles.theta=round(handles.thetaslider(p));
str=sprintf('Theta = %i ',handles.theta);
set(handles.theta_text, 'string', str);

[handles.x1_vec,handles.Pmes_mat,handles.Precons_mat2,erreur_mat]=call_HOA(handles.M,handles.rs,handles.theta,handles.phi,handles.f,handles.r_hp,handles.N_hp);

pcolor(handles.Pmes_graph,handles.x1_vec,handles.x1_vec,real(handles.Pmes_mat));
shading(handles.Pmes_graph, 'interp')
colorbar(handles.Pmes_graph,'location','EastOutside');
cax = caxis(handles.Pmes_graph);
xlabel('Position x (m)')
ylabel('Position y (m)')
% hcb = colorbar('location','EastOutside');
pcolor(handles.Precons_graph,handles.x1_vec,handles.x1_vec,real(handles.Precons_mat2));
colorbar(handles.Precons_graph,'location','EastOutside');
shading(handles.Precons_graph, 'interp')
xlabel('Position x (m)')
ylabel('Position y (m)')
hold(handles.Precons_graph,'on')
contour(handles.Precons_graph,handles.x1_vec,handles.x1_vec,erreur_mat,[0 10])
caxis(handles.Precons_graph,cax)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function theta_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theta_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function phi_slider_Callback(hObject, eventdata, handles)
% hObject    handle to phi_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
D = get(hObject,'value');
p = round(D*(handles.N-1)+1);
handles.phi = round(handles.phislider(p));
str=sprintf('Phi = %i ',handles.phi);
set(handles.phi_text, 'string', str);
[handles.x1_vec,handles.Pmes_mat,handles.Precons_mat2,erreur_mat]=call_HOA(handles.M,handles.rs,handles.theta,handles.phi,handles.f,handles.r_hp,handles.N_hp);

pcolor(handles.Pmes_graph,handles.x1_vec,handles.x1_vec,real(handles.Pmes_mat));
shading(handles.Pmes_graph, 'interp')
colorbar(handles.Pmes_graph,'location','EastOutside');
cax = caxis(handles.Pmes_graph);
xlabel('Position x (m)')
ylabel('Position y (m)')
% hcb = colorbar('location','EastOutside');
pcolor(handles.Precons_graph,handles.x1_vec,handles.x1_vec,real(handles.Precons_mat2));
colorbar(handles.Precons_graph,'location','EastOutside');
shading(handles.Precons_graph, 'interp')
xlabel('Position x (m)')
ylabel('Position y (m)')
hold(handles.Precons_graph,'on')
contour(handles.Precons_graph,handles.x1_vec,handles.x1_vec,erreur_mat,[0 10])
caxis(handles.Precons_graph,cax)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function phi_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phi_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
