function varargout = volCond_gui(varargin)
% VOLCOND_GUI MATLAB code for volCond_gui.fig
%      VOLCOND_GUI, by itself, creates a new VOLCOND_GUI or raises the existing
%      singleton*.
%
%      H = VOLCOND_GUI returns the handle to a new VOLCOND_GUI or the handle to
%      the existing singleton*.
%
%      VOLCOND_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VOLCOND_GUI.M with the given input arguments.
%
%      VOLCOND_GUI('Property','Value',...) creates a new VOLCOND_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before volCond_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to volCond_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help volCond_gui

% Last Modified by GUIDE v2.5 19-Apr-2015 07:57:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @volCond_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @volCond_gui_OutputFcn, ...
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


% --- Executes just before volCond_gui is made visible.
function volCond_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to volCond_gui (see VARARGIN)

% Choose default command line output for volCond_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes volCond_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);
global rQ Q thetaR NRanP sigma Id

% rQ = [0,0, get(handles.rQSlider, 'Value')] * 0.092;
rQ = [0,0, get(handles.rQSlider, 'Value')];
QAngle = get(handles.QAngleSlider, 'Value');
Q = [0,sind(QAngle), cosd(QAngle)];

sigma = 0.33; set(handles.sSigma, 'String', num2str(sigma)); 
Id    = 10.0e-9; set(handles.sId, 'String', num2str(Id/1e-9)); 

set(handles.sQAngle, 'String', get(handles.QAngleSlider, 'Value'));
set(handles.srQ, 'String', get(handles.rQSlider, 'Value'));

thetaR = 1; set(handles.sThetaRanRad, 'String', num2str(thetaR)); 
NRanP = 50; set(handles.sNRangePoints, 'String', num2str(NRanP));

set( gcf, 'toolbar', 'figure' )

bShowCurrents_Callback(hObject, eventdata, handles)
% This has to be run explicitly (despite it already being at the end of
% bShowCurrents_Callback). Seems like there's something fishy about ML's
% initialisation, helps re-running the first time? Works anyway...
drawFieldPlots(handles)


function drawFieldPlots(handles)
global sigma Id rQ Q bShowCurrents

M=Id/(4*pi*sigma);
% axes(handles.axes1); cla
plotSphField('infiniteSpace', rQ, Q, 1, M, bShowCurrents,handles.axes1)
% infiniteSpace_BCK(rQ, Q, plot_streams)
% axes(handles.axes2); cla
% homogSphere_BCK(rQ, Q, plot_streams)
plotSphField('homogSphere', rQ, Q, 1, M, bShowCurrents, handles.axes2)
set(handles.axes1, 'ButtonDownFcn', {@axes1_ButtonDownFcn,handles});
set(handles.axes2, 'ButtonDownFcn', {@axes2_ButtonDownFcn,handles});
% set(handles.axes1, 'PickableParts', 'all');
% set(handles.axes2, 'PickableParts', 'all');

% --- Outputs from this function are returned to the command line.
function varargout = volCond_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function rQSlider_Callback(hObject, eventdata, handles)
% hObject    handle to rQSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global rQ
rQ = [0,0, round(get(handles.rQSlider, 'Value')*100)/100];
set(handles.srQ, 'String', rQ(3));

drawFieldPlots(handles)

% axes(handles.axes1); cla
% infiniteSpace(rQ,Q,0)
% axes(handles.axes2); cla
% homogSphere(rQ,Q,0)
% set(handles.axes1, 'ButtonDownFcn', {@axes1_ButtonDownFcn,handles});
% set(handles.axes2, 'ButtonDownFcn', {@axes2_ButtonDownFcn,handles});



% --- Executes during object creation, after setting all properties.
function rQSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rQSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function QAngleSlider_Callback(hObject, eventdata, handles)
% hObject    handle to QAngleSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global rQ Q
QAngle = round(get(handles.QAngleSlider, 'Value'));
set(handles.sQAngle, 'String', num2str(QAngle));
Q = [0,sind(QAngle), cosd(QAngle)];

drawFieldPlots(handles)

% --- Executes during object creation, after setting all properties.
function QAngleSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to QAngleSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function sQAngle_Callback(hObject, eventdata, handles)
% hObject    handle to sQAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sQAngle as text
%        str2double(get(hObject,'String')) returns contents of sQAngle as a double
QAngle = round(str2double(get(handles.sQAngle, 'String')));
set(handles.QAngleSlider, 'Value', QAngle);
QAngleSlider_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function sQAngle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sQAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function srQ_Callback(hObject, eventdata, handles)
% hObject    handle to srQ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of srQ as text
%        str2double(get(hObject,'String')) returns contents of srQ as a double


% --- Executes during object creation, after setting all properties.
function srQ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to srQ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in bSelRange.
function bSelRange_Callback(hObject, eventdata, handles)
% hObject    handle to bSelRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bSelRange
global range_points 
range_points = [];

% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes2_ButtonDownFcn(hObject, eventdata, handles)

% --- Executes on mouse press over axes background.
function axes2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global range_points thetaR NRanP rQ Q Id sigma bShowCurrents

clear coordinates
coordinates = get(hObject,'CurrentPoint');

cols = ['b','r'];
linestyle = '-';

M=Id/(4*pi*sigma);

if get(handles.bSelLinRange,'Value') > 0 || ... 
     get(handles.bSelThetaRange,'Value') > 0 % select range ON

    % if theta-range, fix to sphere!
    if get(handles.bSelThetaRange,'Value') > 0
         R = sqrt(dot(coordinates(1,1:2), coordinates(1,1:2)));
         coordinates(1,1:2) = coordinates(1,1:2)/R * thetaR;
         linestyle = 'none';
    end
    range_points = [range_points; coordinates(1,1:2)];
    for ii = 1:2,
        curax = eval(['handles.axes' num2str(ii)]);
        if size(range_points,1) == 1
            if bShowCurrents
                Nstreams = 15; % hard-coded here from plotSphField.m!
                to_delete = (ii+2) + Nstreams ;% left are arrow, circle (line) & contour + currents
            else
                to_delete = (ii+2);% arrow, circle (line) & contour
            end                
            children = get(curax,'Children');
            for jj = 1:length(children)-to_delete
%             for jj = 1:length(children)-(2) % left are contour & circle!
                delete(children(jj))
            end
        end
        plot(curax,range_points(:,1), range_points(:,2),'md',...
            'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','m',...
            'LineStyle',linestyle);

    end
    
    if size(range_points,1) > 1,

        for ii = 1:2,
            curax = eval(['handles.axes' num2str(ii)]);

            % if theta-range, fix to sphere!
            if get(handles.bSelThetaRange,'Value') > 0
                % these are -pi -> pi
                Th = atan2(range_points(:,2), range_points(:,1));
                % convert -pi/2 -> 2/3 pi, then re-wrap
                Th(Th<0) = Th(Th<0) + 2*pi;
                Th = Th + pi/2;
                Th(Th>2*pi) = Th(Th>2*pi) - 2*pi;

                Th_start = min(Th) - pi/2;
                Th_stop = max(Th) - pi/2;
                Th_step = (Th_stop-Th_start) / NRanP;
                theta = [Th_start:Th_step:Th_stop];
                %theta = linspace(Th_start,Th_stop,nTheta);
                plot(curax, thetaR*cos(theta), thetaR*sin(theta), ...
                    'Color','m', ...
                    'LineWidth', 2, 'Hittest', 'off',...
                    'LineStyle','-.');
                Y = cos(fliplr(theta)); Z = sin(fliplr(theta));
%                 x = pi/2-Th_stop:Th_step:pi/2-Th_start;
                x = (pi/2-fliplr(theta))/(2*pi)*360;
                xlab = 'Theta (degrees)';
            elseif get(handles.bSelLinRange,'Value') > 0
                Y = linspace(range_points(1,1), range_points(2,1), NRanP);
                Z = linspace(range_points(1,2), range_points(2,2), NRanP);
                x = sqrt((Y-range_points(1,1)).^2 + (Z-range_points(1,2)).^2);
                xlab = 'Distance on line (a.u.)';
            end
            figure(101); clf; hold on; grid on
            V_is = infiniteSpace(rQ, Q, thetaR, M, Y, Z);
            V_hs = homogSphere(rQ, Q, thetaR, M, Y, Z);
            range_rad = sqrt(Y.^2 + Z.^2);
            V_hs(range_rad > 1) = 0;
            %V = calc_field(src_vec, sig, X, Y);
            %d = sqrt(X.^2 + Y.^2);
            plot(x, 1e9*V_is, 'Color', [0 0.4470 0.7410], 'Linewidth', 2);
            plot(x, 1e9*V_hs, 'Color', [0.8500 0.3250 0.0980], 'Linewidth', 2);
            hl = legend('infiniteSpace','homogSphere', 'Location', 'best');
            set(hl,'FontSize',14);
            xlabel(xlab)
            ylabel('Electric potential (nV)')
            set(gca, 'FontSize', 14)
        end
        range_points = [];
        set(handles.bSelLinRange,'Value',0)
        set(handles.bSelThetaRange,'Value',0)
    end
end

% --- Executes on button press in bSelLinRange.
function bSelLinRange_Callback(hObject, eventdata, handles)
% hObject    handle to bSelLinRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bSelLinRange
global range_points 
range_points = [];

% --- Executes on button press in bSelThetaRange.
function bSelThetaRange_Callback(hObject, eventdata, handles)
% hObject    handle to bSelThetaRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bSelThetaRange
global range_points 
range_points = [];



function sThetaRanRad_Callback(hObject, eventdata, handles)
% hObject    handle to sThetaRanRad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sThetaRanRad as text
%        str2double(get(hObject,'String')) returns contents of sThetaRanRad as a double
global thetaR
thetaR = str2num(get(handles.sThetaRanRad, 'String'));



function sNRangePoints_Callback(hObject, eventdata, handles)
% hObject    handle to sNRangePoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sNRangePoints as text
%        str2double(get(hObject,'String')) returns contents of sNRangePoints as a double
global NRanP
NRanP = str2num(get(handles.sNRangePoints, 'String'));



function sSigma_Callback(hObject, eventdata, handles)
% hObject    handle to sSigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sSigma as text
%        str2double(get(hObject,'String')) returns contents of sSigma as a double
global sigma
sigma = str2double(get(handles.sSigma, 'String'));
rQSlider_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function sSigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sSigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sId_Callback(hObject, eventdata, handles)
% hObject    handle to sId (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sId as text
%        str2double(get(hObject,'String')) returns contents of sId as a double
global Id
Id = 1e-9 * str2double(get(handles.sId, 'String'));
rQSlider_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function sId_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sId (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in bShowCurrents.
function bShowCurrents_Callback(hObject, eventdata, handles)
% hObject    handle to bShowCurrents (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bShowCurrents
global bShowCurrents

bShowCurrents = get(handles.bShowCurrents,'Value');
drawFieldPlots(handles)
