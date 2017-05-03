function varargout = leadField_gui(varargin)
% LEADFIELD_GUI MATLAB code for leadField_gui.fig
%      LEADFIELD_GUI, by itself, creates a new LEADFIELD_GUI or raises the existing
%      singleton*.
%
%      H = LEADFIELD_GUI returns the handle to a new LEADFIELD_GUI or the handle to
%      the existing singleton*.
%
%      LEADFIELD_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LEADFIELD_GUI.M with the given input arguments.
%
%      LEADFIELD_GUI('Property','Value',...) creates a new LEADFIELD_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before leadField_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to leadField_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help leadField_gui

% Last Modified by GUIDE v2.5 12-Apr-2015 11:50:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @leadField_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @leadField_gui_OutputFcn, ...
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


% --- Executes just before leadField_gui is made visible.
function leadField_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to leadField_gui (see VARARGIN)

% Choose default command line output for leadField_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes leadField_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);
global R sigmas

R = [8.0, 8.1, 8.6, 9.2]; % cm
sigmas = [0.33, 1, 0.0042, 0.33];
% boundary_cols = [0.5,0.5,0;...
%     0,0.5,0.5];

set(handles.Rbrain, 'String', num2str(R(1)));
set(handles.RCSF, 'String', num2str(R(2)));
set(handles.Rskull, 'String', num2str(R(3)));
set(handles.Rscalp, 'String', num2str(R(4)));
set(handles.Sbrain, 'String', num2str(sigmas(1)));
set(handles.SCSF, 'String', num2str(sigmas(2)));
set(handles.Sskull, 'String', num2str(sigmas(3)));
set(handles.Sscalp, 'String', num2str(sigmas(4)));

set(handles.sensorRscalp, 'String', num2str(R(4)));
%set(handles.sensorRplus, 'String', num2str(0));
Rsensor = str2double(get(handles.sensorRscalp, 'String')) + ...
    str2double(get(handles.sensorRplus, 'String'));
set(handles.sensorR, 'String', num2str(Rsensor));

% thetaEl = [-pi/8; pi/8];
% rSens = R(end) * [zeros(size(thetaEl)), sin(thetaEl), cos(thetaEl)];
% wSens = [-1, 1]; erSens = [1, 1];
% senstype = 'EEG';
%plot_lead_ND(handles, rSens/100, wSens, erSens, R/100, sigmas, senstype)


set( gcf, 'toolbar', 'figure' )

clear_sensors(handles)

function clear_sensors(handles)

global rSens wSens erSens

axes(handles.axes2); cla
axes(handles.axes1); cla

change_circle_radii(handles)

rSens = []; wSens = []; erSens=[];

function change_circle_radii(handles)

global R

boundary_cols = [0.5,0.5,0;...
                 0,0.3,0.5*sqrt(2);...
                 0.7,0,0.7;...
                 0,0.5,0.5];
lw = [2,0.5,4,2];
ls = {'-','--','-','-'};

children = get(handles.axes1,'Children');
for jj = 1:length(children)
    if strcmp(get(children(jj), 'Type'), 'line')
        delete(children(jj))
    end
end

Rsensor = str2double(get(handles.sensorR, 'String'));

axes(handles.axes1)
for ii = 1:length(R)
    add_circle(R(ii), ls{ii}, lw(ii), boundary_cols(ii,:))
end
% add_circle(R(4), '-', 4, boundary_cols(2,:))
add_circle(Rsensor, '--', 2, 'c')
axis tight; axis xy


% --- Outputs from this function are returned to the command line.
function varargout = leadField_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Rbrain_Callback(hObject, eventdata, handles)
% hObject    handle to Rbrain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Rbrain as text
%        str2double(get(hObject,'String')) returns contents of Rbrain as a double
global R 
R(1) = str2double(get(handles.Rbrain, 'String'));
refresh_leads(handles)
change_circle_radii(handles)


% --- Executes during object creation, after setting all properties.
function Rbrain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rbrain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RCSF_Callback(hObject, eventdata, handles)
% hObject    handle to RCSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RCSF as text
%        str2double(get(hObject,'String')) returns contents of RCSF as a double
global R 
R(2) = str2double(get(handles.RCSF, 'String'));
refresh_leads(handles)
change_circle_radii(handles)


% --- Executes during object creation, after setting all properties.
function RCSF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RCSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Rskull_Callback(hObject, eventdata, handles)
% hObject    handle to Rskull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Rskull as text
%        str2double(get(hObject,'String')) returns contents of Rskull as a double
global R 
R(3) = str2double(get(handles.Rskull, 'String'));
refresh_leads(handles)
change_circle_radii(handles)


% --- Executes during object creation, after setting all properties.
function Rskull_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rskull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Rscalp_Callback(hObject, eventdata, handles)
% hObject    handle to Rscalp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Rscalp as text
%        str2double(get(hObject,'String')) returns contents of Rscalp as a double
global R 
R(4) = str2double(get(handles.Rscalp, 'String'));
refresh_leads(handles)
change_circle_radii(handles)


% --- Executes during object creation, after setting all properties.
function Rscalp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rscalp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Sbrain_Callback(hObject, eventdata, handles)
% hObject    handle to Sbrain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Sbrain as text
%        str2double(get(hObject,'String')) returns contents of Sbrain as a double
global sigmas
sigmas(1) = str2double(get(handles.Sbrain, 'String'));
refresh_leads(handles)


% --- Executes during object creation, after setting all properties.
function Sbrain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Sbrain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SCSF_Callback(hObject, eventdata, handles)
% hObject    handle to SCSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SCSF as text
%        str2double(get(hObject,'String')) returns contents of SCSF as a double
global sigmas
sigmas(2) = str2double(get(handles.SCSF, 'String'));
refresh_leads(handles)


% --- Executes during object creation, after setting all properties.
function SCSF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SCSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Sskull_Callback(hObject, eventdata, handles)
% hObject    handle to Sskull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Sskull as text
%        str2double(get(hObject,'String')) returns contents of Sskull as a double
global sigmas
sigmas(3) = str2double(get(handles.Sskull, 'String'));
refresh_leads(handles)


% --- Executes during object creation, after setting all properties.
function Sskull_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Sskull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Sscalp_Callback(hObject, eventdata, handles)
% hObject    handle to Sscalp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Sscalp as text
%        str2double(get(hObject,'String')) returns contents of Sscalp as a double
global sigmas
sigmas(4) = str2double(get(handles.Sscalp, 'String'));
refresh_leads(handles)


% --- Executes during object creation, after setting all properties.
function Sscalp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Sscalp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in sensorPopup.
function sensorPopup_Callback(hObject, eventdata, handles)
% hObject    handle to sensorPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sensorPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sensorPopup
global senstype
contents = cellstr(get(hObject,'String'));

senstype = contents{get(hObject,'Value')};
if strcmp(senstype(1:3), 'MEG')
    set(handles.sensorRplus, 'String', num2str(3));
    leadField_gui_OpeningFcn(hObject, eventdata, handles);
    sensorRplus_Callback(hObject, eventdata, handles);
elseif strcmp(senstype(1:3), 'EEG')
    set(handles.sensorRplus, 'String', num2str(0));
    leadField_gui_OpeningFcn(hObject, eventdata, handles);
    sensorRplus_Callback(hObject, eventdata, handles);
end
% --- Executes during object creation, after setting all properties.
function sensorPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensorPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
global senstype

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
items = {'EEG', 'MEG-mag', 'MEG-xgrad', 'MEG-ygrad'};
set(hObject,'String', items)

senstype = 'EEG';


function sensorRplus_Callback(hObject, eventdata, handles)
% hObject    handle to sensorRplus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sensorRplus as text
%        str2double(get(hObject,'String')) returns contents of sensorRplus as a double
global rSens

Rsensor = str2double(get(handles.sensorRscalp, 'String')) + ...
    str2double(get(handles.sensorRplus, 'String'));
set(handles.sensorR, 'String', num2str(Rsensor));
if size(rSens,1) > 0
    cur_Rsensors = sqrt(dot(rSens,rSens,2));
    rSens = rSens * (Rsensor/cur_Rsensors(1));
    refresh_leads(handles)
end
change_circle_radii(handles)


% --- Executes during object creation, after setting all properties.
function sensorRplus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sensorRplus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in bAddSensor.
function bAddSensor_Callback(hObject, eventdata, handles)
% hObject    handle to bAddSensor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bAddSensor

% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global rSens erSens wSens senstype 

clear coordinates
coordinates = get(hObject,'CurrentPoint');

if get(handles.bAddSensor,'Value') > 0

    % fix to sphere!
    Rpoint = sqrt(dot(coordinates(1,1:2), coordinates(1,1:2)));
    Rsensor = str2double(get(handles.sensorR, 'String'));
    coordinates(1,1:2) = coordinates(1,1:2)/Rpoint * Rsensor;
    sensorPos = [0, coordinates(1,1:2)]; % x == 0
    thetaSens = atan2(sensorPos(2),sensorPos(3));
%     disp(thetaSens/(2*pi) * 360)    
    
    sensSize = 12;
    if strcmp(senstype(1:3), 'EEG')
        rSens = [rSens; sensorPos]; % x == 0
        erSens = [erSens; [0,0,0]]; % not used
        if size(rSens,1) < 2 % first pick is ref
            sensCol = 'b';
            wSens = [-1];
        else
            sensCol = 'r';
%             wSens = [-1, ones(1,size(rSens,1)-1) ./ (size(rSens,1) - 1)];
            wSens = [-1, ones(1,size(rSens,1)-1)];
            
        end
    elseif strcmp(senstype(1:3), 'MEG')        
        sensCol = 'c';
        thetaMagRot = -thetaSens; % zero-angle is on y-axis, opens counter-cl.
        
        coilname = senstype(5:end);
        [rC, w, er] = coildefs(coilname);

        % first rotate in place
        rotMat = [1,            0,             0; ...
            0, cos(thetaMagRot), -sin(thetaMagRot); ...
            0, sin(thetaMagRot), cos(thetaMagRot)];
        
        rC = (rotMat * rC')';
        er = (rotMat * er')';
        
        % move center to coil position
        Rm = Rsensor * [zeros(size(thetaSens)), sin(thetaSens), cos(thetaSens)];
        
        rC = rC + repmat(Rm, size(rC, 1), 1);
        
        rSens = [rSens; rC];
        wSens = [wSens, w];
        erSens = [erSens; er];
        %disp(wSens)
    end
    axes(handles.axes1)
    hold on
    plot(rSens(end,2), rSens(end,3), sensCol,'MarkerFaceColor', sensCol,...
        'Marker', 'o', 'MarkerSize', sensSize, 'Hittest','off')
end

if strcmp(senstype(1:3), 'MEG') || ...
    (strcmp(senstype(1:3), 'EEG') && size(rSens,1) > 1)

    refresh_leads(handles);
end

function refresh_leads(handles)
global rSens wSens erSens R sigmas senstype

% is zero for mag and ygrad, but numerical roundoff gives noisy
% mag-plot, so just skip the quiver for now (doesn't look good anyway)
if (strcmp(senstype(1:3), 'MEG') && strcmp(senstype(5:end), 'xgrad')) || ...
        strcmp(senstype(1:3), 'EEG')
    b2Dquiver = 1;
else
    b2Dquiver = 0;
end    

if (strcmp(senstype(1:3), 'MEG') && size(rSens,1) > 0) || ...
    (strcmp(senstype(1:3), 'EEG') && size(rSens,1) > 1)

    plot_lead_ND(handles, rSens, wSens, erSens, R, sigmas, senstype, b2Dquiver);
    if strcmp(senstype, 'EEG')
        axes(handles.axes2); 
        hold on
        sh=plot3(rSens(1,1),rSens(1,2), rSens(1,3), 'bo', ...
            'MarkerSize', 12, 'MarkerFaceColor','b', 'Hittest','off');
        sh=plot3(rSens(2:end,1),rSens(2:end,2), rSens(2:end,3), 'ro', ...
            'MarkerSize', 12, 'MarkerFaceColor','r', 'Hittest','off');
        hold off
    else
        nSensors = size(rSens,1)/4;
        for jj = 1:nSensors
            ran = (jj-1)*4 + 1:jj*4;
            Rm = mean(rSens(ran,:),1)/100; % all MEG integration points are symmetric in x!
            thetaMagRot = -atan2(erSens(ran(1),2),erSens(ran(1),3));
            rotMat = [1,            0,             0; ...
                0, cos(thetaMagRot), -sin(thetaMagRot); ...
                0, sin(thetaMagRot), cos(thetaMagRot)];
            
            % %         Rm = [0,0,0.080];
            %         rotMat = eye(3);
            coilframe = 100*coil_outline(Rm, rotMat, senstype(5:end));
            hold on
            sh=plot3(coilframe(:,1),coilframe(:,2),coilframe(:,3),'Linewidth',2);
            set(sh,'Color',[0.8500    0.3250    0.0980]);
            hold off
        end
        Rsensor = str2double(get(handles.sensorR, 'String'));

        set(gca,'ZLim', 1.01*[-Rsensor,Rsensor])
        % coilfram
    end

    set(handles.axes1, 'ButtonDownFcn', {@axes1_ButtonDownFcn,handles});
    set(handles.axes2, 'ButtonDownFcn', {@axes2_ButtonDownFcn,handles});
%     set(handles.axes1, 'PickableParts', 'all');
%     set(handles.axes2, 'PickableParts', 'all');
end

% --- Executes on mouse press over axes background.
function axes2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ClearBut.
function ClearBut_Callback(hObject, eventdata, handles)
% hObject    handle to ClearBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear_sensors(handles)
