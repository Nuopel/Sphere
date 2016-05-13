function varargout = Guifig(varargin)
% GUIFIG MATLAB code for Guifig.fig
%      GUIFIG, by itself, creates a new GUIFIG or raises the existing
%      singleton*.
%
%      H = GUIFIG returns the handle to a new GUIFIG or the handle to
%      the existing singleton*.
%
%      GUIFIG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIFIG.M with the given input arguments.
%
%      GUIFIG('Property','Value',...) creates a new GUIFIG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Guifig_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Guifig_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Guifig

% Last Modified by GUIDE v2.5 10-May-2016 14:57:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Guifig_OpeningFcn, ...
    'gui_OutputFcn',  @Guifig_OutputFcn, ...
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

% --- Executes just before Guifig is made visible.
function Guifig_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Guifig (see VARARGIN)

% Choose default command line output for Guifig
handles.output = hObject;
handles.ct.M=5;
handles.N = 1001;
handles.sliderfreq = logspace(log10(100),log10(5000),handles.N);
handles.ct.r_micsph=0.07;
handles.ct.hankel_order=2;

[~, p] = min(abs(handles.sliderfreq-(5000/10^0.5)));
set(handles.slider_freq, 'value', (p-1)/(length(handles.sliderfreq)-1));
set(handles.TextFreq, 'string', ['Frequency = ' num2str(round(5000/10^0.5)) ]);
set(handles.R, 'string', 'R = '  );
set(handles.Theta, 'string', 'Theta = ' );
set(handles.Phi, 'string', 'Phi = '  );

%% Create antenna 
ct.pas_m = 2.54*2e-2; % pas de l'antenne
N.nbrx_sca =20; % nombre de micros par ligne
N.nbry_sca =20; % nombre de micros par ligne
ct.N_mic=N.nbrx_sca*N.nbry_sca;
[ handles.Antenna ] = AntennArray( ct.pas_m,N.nbrx_sca,N.nbry_sca) ;

%% Define Spherical microphone set up (same as speaker but smaller)
[ handles.Sphmic, ct.N_mic ] = CreateSpeakerSystem(handles.ct.r_micsph);% create the sphere set up, sort in struc Array

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using Guifig.
if strcmp(get(hObject,'Visible'),'off')
    plot(rand(5));
end

% UIWAIT makes Guifig wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Guifig_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.Target);
cla;

popup_sel_index = get(handles.popupmenu1, 'Value');
switch popup_sel_index
    case 1
        plot(rand(5));
    case 2
        plot(sin(1:0.01:25.99));
    case 3
        bar(1:.5:10);
    case 4
        plot(membrane);
    case 5
        surf(peaks);
end


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
    ['Close ' get(handles.figure1,'Name') '...'],...
    'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


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

set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});


% --- Executes on slider movement.
function slider_freq_Callback(hObject, eventdata, handles)
% hObject    handle to slider_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
D = get(hObject,'value');
p = round(D*(handles.N-1)+1);
handles.slider = round(handles.sliderfreq(p));
str=sprintf('Frequency = %i ',handles.slider);
set(handles.TextFreq, 'string', str);


if isfield(handles,'ct')
if  isfield(handles.ct,'R')==1 && isfield(handles.ct,'Phi')==1 && isfield(handles.ct,'Theta')==1
        handles.ct.k = 2*pi*round(handles.sliderfreq(p))/340;
        % calculate targets
        [ Pressure ] = CallProcessingMicSphTh( handles.Antenna,handles.ct  );
        %Plot 
        var.empty=[];
        var=Pressure_map_gui(handles.Target,real(Pressure.monopole),handles.ct,handles.Antenna,var );
       
        Pressure_map_SphMic_gui(handles.AmbiTarget,handles.ct.M,Pressure.TargetAmbisonics,handles.ct,var,handles.Antenna);
        

end
end

if isfield(handles,'data') 

    if  isfield(handles.data,'OutSig')==1  && isfield(handles.data,'EntrySig')==1 && isfield(handles,'Bmn')==0 %%&& isfield(handles.data,'CalibMic')==1
        [handles.H_Data,handles.t,handles.var  ]= CallProcessingMicSphMeas( handles.data,handles.ct  );
         handles.H_Data.h_sig_fft=fft(handles.data.OutSig(1:240002,:));
    end
    ct.pos = closest( handles.ct.k,handles.t.Fsweep_avg*2*pi/340 );
    handles.ct.N_mic=50;
    handles.Bmn.recons = Bmn_encoding_sph( handles.H_Data.h_sig_fft(ct.pos,:),handles.Sphmic,handles.ct,handles.var );
%     handles.Bmn.recons = Bmn_encoding_sph( Pressure.difract,handles.Sphmic,handles.ct,handles.var );

    N.N_sweep=1;handles.ct.N_mic=handles.Antenna.N_mic;
    Pressure.decodemeas = Decoding_pressure_field(handles.ct.M,handles.Bmn.recons,handles.Antenna,handles.ct,handles.var,N ) ;
    
    Pressure_map_SphMic_gui(handles.Meas,handles.ct.M, real(Pressure.decodemeas),handles.ct,handles.var,handles.Antenna);
else
    disp('Load all data before plotting results')
    
    

end
 guidata(hObject, handles);   


% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in ListOrder.
function ListOrder_Callback(hObject, eventdata, handles)
% hObject    handle to ListOrder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents=cellstr(get(hObject,'String'));
choix=cell2mat({get(hObject,'Value')});
handles.M=str2double(contents(choix));

% Hints: contents = cellstr(get(hObject,'String')) returns ListOrder contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ListOrder
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ListOrder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ListOrder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in EntryFile.
function EntryFile_Callback(hObject, eventdata, handles)
% hObject    handle to EntryFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile('*.wav','File selector');
text.fullpath_entry = strcat( PathName, FileName );
[ handles.data.EntrySig, handles.ct.Fs ] = audioread(text.fullpath_entry);
PathWriter('Default\in.txt',text.fullpath_entry)
disp('In Data loaded')

guidata(hObject, handles);

% --- Executes on button press in OutFile.
function OutFile_Callback(hObject, eventdata, handles)
% hObject    handle to OutFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile('*.w64','File selector');
text.fullpath_out = strcat( PathName, FileName );
[ handles.data.OutSig, ct.Fs ] = audioread(text.fullpath_out);
if handles.ct.Fs~=ct.Fs
    disp('Error the sampling frequency is not the same for input and output');
    return;
end
if isfield(handles,'Bmn')
handles=rmfield(handles,'Bmn');
end
disp('Out Data loaded')
PathWriter('Default\out.txt',text.fullpath_out)

guidata(hObject, handles);

% --- Executes on button press in CalibMic.
function CalibMic_Callback(hObject, eventdata, handles)
% hObject    handle to CalibMic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile('*.w64','File selector');
text.fullpath_CalibMic = strcat( PathName, FileName );
[ handles.data.CalibMic, ct.Fs ] = audioread(text.fullpath_CalibMic);
if handles.ct.Fs~=ct.Fs
    disp('Error the sampling frequency is not the same for input and output');
    return;
end

PathWriter('Default\calibmic.txt',text.fullpath_CalibMic)
guidata(hObject, handles)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile('*.w64','File selector');
text.fullpathCalibFRF = strcat( PathName, FileName );
[ handles.data.CalibFrf, ct.Fs ] = audioread(text.fullpathCalibFRF);
if handles.ct.Fs~=ct.Fs
    disp('Error the sampling frequency is not the same for input and output');
    return;
end
PathWriter(filename,path)
guidata(hObject, handles);


% --- Executes on button press in Default.
function Default_Callback(hObject, eventdata, handles)
% hObject    handle to Default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles.data.OutSig = load('Default\out.txt');
% % handles.data.CalibMic = load('Default/calib.txt');
% handles.data.EntrySig = load('Default\in.txt');
handles.ct.R=load('Default/R.txt');
handles.ct.Theta=load('Default/Theta.txt');
handles.ct.Phi=load('Default/Phi.txt');


guidata(hObject, handles);



function PhiBox_Callback(hObject, eventdata, handles)
% hObject    handle to PhiBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ct.Phi = get(hObject,'String');
handles.ct.Phi = str2double(handles.ct.Phi);
PathWriter('Default/Phi.txt',sprintf('%i',handles.ct.Phi));


% Hints: get(hObject,'String') returns contents of PhiBox as text
%        str2double(get(hObject,'String')) returns contents of PhiBox as a double
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function PhiBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PhiBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ThetaBox_Callback(hObject, eventdata, handles)
% hObject    handle to ThetaBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ct.Theta = get(hObject,'String');
handles.ct.Theta = str2double(handles.ct.Theta);
PathWriter('Default/Theta.txt',sprintf('%i',handles.ct.Theta));

% Hints: get(hObject,'String') returns contents of ThetaBox as text
%        str2double(get(hObject,'String')) returns contents of ThetaBox as a double
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ThetaBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ThetaBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RBox_Callback(hObject, eventdata, handles)
% hObject    handle to RBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ct.R = get(hObject,'String');

handles.ct.R = str2double(handles.ct.R);
PathWriter('Default/R.txt',sprintf('%i',handles.ct.R));

% Hints: get(hObject,'String') returns contents of RBox as text
%        str2double(get(hObject,'String')) returns contents of RBox as a double
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function RBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
    
