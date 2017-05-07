function varargout = ecd_gui(varargin)
% ECD_GUI MATLAB code for ecd_gui.fig
%      ECD_GUI, by itself, creates a new ECD_GUI or raises the existing
%      singleton*.
%
%      H = ECD_GUI returns the handle to a new ECD_GUI or the handle to
%      the existing singleton*.
%
%      ECD_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ECD_GUI.M with the given input arguments.
%
%      ECD_GUI('Property','Value',...) creates a new ECD_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ecd_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ecd_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ecd_gui

% Last Modified by GUIDE v2.5 27-Mar-2015 14:43:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ecd_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @ecd_gui_OutputFcn, ...
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


% --- Executes just before ecd_gui is made visible.
function ecd_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ecd_gui (see VARARGIN)

% Choose default command line output for ecd_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global current_sink_source totcurrent src_vec add_amps sig elec_coords
current_sink_source = -1;
elec_coords = -1;
totcurrent = 0;
src_vec = []; % [x, y, I]
set(handles.totcurrent,'String',num2str(totcurrent))
set(handles.axes1,'Xlim', [-1, 1])
set(handles.axes1,'Ylim', [-2, 2])
set(handles.axes2,'Xlim', [-2, 2])
set(handles.axes2,'Ylim', [-2, 2])
% reset handle values to default
set(handles.sinkbut, 'Value', 1.0)
set(handles.sourcebut, 'Value', 0.0)
add_amps = str2num(get(handles.sinkstren,'String'));

img=imread('pyrasketch.jpg');
xx=linspace(-1,1,size(img,2));
yy=linspace(-2,2,size(img,1));
% To prevent problems with older versions of flipud, we'll manually flip
% the first dimension
img = img(end:-1:1,:,:);
imh = imagesc(xx,yy,img, 'Parent', handles.axes1, 'Hittest','off');% hold on
axis(handles.axes1, 'xy')
hold on

% WOW: This is a BIG gotcha'! 
% Any plotting command will by default replace all the callbacks!
% To prevent calls to PLOT, BAR, etc. from replacing the buttondownfcn, set the 'NextPlot' property of the axes to 'replacechildren' instead of the default 'replace'. This will clear the axes contents and plot the new data without replacing the axes object itself.
% If you want to plot data over your existing data, then set the 'NextPlot' property of the axes to 'Add'.
% An alternative workaround is to assign 'ButtonDownFcn' for the axes again after calling the plot.
set(handles.axes1, 'ButtonDownFcn', {@axes1_ButtonDownFcn,handles});

sig = 0.33;
set(handles.sigmaStr, 'String', num2str(sig));

% elecCoordsStr = strsplit(get(handles.elecCoordsStr,'String'), ' ; ');
% elec_coords(1) = str2double(elecCoordsStr(1));
% elec_coords(2) = str2double(elecCoordsStr(2));
% set(handles.elecCoordsStr, 'String', '')
% set(handles.elecValStr, 'String', '')

cla(handles.axes2)
% set(handles.axes2, 'PickableParts', 'all');

% UIWAIT makes ecd_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ecd_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in calcV.
function calcV_Callback(hObject, eventdata, handles)
% hObject    handle to calcV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)    elseif strcmp(field_type{ff}, 'open')
global src_vec sig clim_locked elec_coords
obsDist = str2double(get(handles.obsDist, 'String'));
clim_auto = plot_field(src_vec,obsDist, sig, handles.axes2);
set(handles.axes2, 'ButtonDownFcn', {@axes2_ButtonDownFcn,handles});
% set(handles.axes2, 'PickableParts', 'all');
set(handles.elecCoordsStr, 'String', '')
set(handles.elecValStr, 'String', '')
if length(clim_locked) > 1
    set(handles.axes2, 'CLim', clim_locked)
else
    set(handles.axes2, 'CLim', clim_auto)
end

if length(elec_coords) > 1
    ph = plot(elec_coords(1), elec_coords(2), 'x', 'HitTest','off',...
        'MarkerSize',30,'MarkerEdgeColor','g',...
        'LineWidth', 3);
    
    V = calc_field(src_vec, sig, elec_coords(1), elec_coords(2));
    
    elecValStr = sprintf('%.2f', V);
    set(handles.elecValStr, 'String', elecValStr)
    elecCoordsStr = sprintf('%.1f ; %.1f', elec_coords(1), elec_coords(2));
    set(handles.elecCoordsStr, 'String', elecCoordsStr)

end

% --- Executes on button press in sinkbut.
function sinkbut_Callback(hObject, eventdata, handles)
% hObject    handle to sinkbut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sinkbut
global current_sink_source add_amps
current_sink_source = -1;
add_amps = str2num(get(handles.sinkstren,'String'));

% --- Executes on button press in sourcebut.
function sourcebut_Callback(hObject, eventdata, handles)
% hObject    handle to sourcebut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sourcebut
global current_sink_source add_amps
current_sink_source = 1;
add_amps = str2num(get(handles.sourcestren,'String'));


function sinkstren_Callback(hObject, eventdata, handles)
% hObject    handle to sinkstren (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sinkstren as text
%        str2double(get(hObject,'String')) returns contents of sinkstren as a double
global add_amps
add_amps = str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function sinkstren_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sinkstren (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sourcestren_Callback(hObject, eventdata, handles)
% hObject    handle to sourcestren (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sourcestren as text
%        str2double(get(hObject,'String')) returns contents of sourcestren as a double
global add_amps
add_amps = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function sourcestren_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sourcestren (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global current_sink_source add_amps totcurrent src_vec

clear coordinates
coordinates = get(hObject,'CurrentPoint');
clickPos = coordinates(1,1:2);

% place a circle marking the point
%cla
if current_sink_source < 0 % sink
    col = 'r';
    markerStyle = 'v';
else
    col = 'b';
    markerStyle = '^';
end    
cur_amp = current_sink_source*add_amps;
totcurrent = totcurrent + cur_amp;

src_vec = [src_vec; clickPos, cur_amp];

% added because on ML < 2014b it seems to have been forgotten from the
% openingFunction!
hold on
%text(clickPos(1), clickPos(2), 'O', 'HitTest','off','FontSize',10*abs(cur_amp),'Color',col,'HorizontalAlignment','Center')
ph=plot(clickPos(1), clickPos(2), markerStyle, 'HitTest','off',...
    'MarkerSize',10*abs(cur_amp),'MarkerEdgeColor',col, ...
    'LineWidth', 2, 'MarkerFaceColor',col);
% drawScenario();
%disp([num2str(clickPos(1)) ',' num2str(clickPos(2)), ':' num2str(cur_amp)]);


set(handles.totcurrent,'String',num2str(totcurrent))
% 
% 


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on button press in clearsrc.
function clearsrc_Callback(hObject, eventdata, handles)
% hObject    handle to clearsrc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ecd_gui_OpeningFcn(hObject, eventdata, handles)
% callbackCell = get(handles.sinkbut,'Callback');
% callbackCell{1}(handles.sinkbut,[],callbackCell{2:end});


% --- Executes on mouse press over axes background.
function axes2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global elec_coords src_vec sig range_points

clear coordinates
coordinates = get(hObject,'CurrentPoint');

children = get(handles.axes2,'Children');
for ii = 1:length(children)-1
    delete(children(ii))
end

if get(handles.bSelRange,'Value') < 1 % select range OFF
    elec_coords = coordinates(1,1:2);
    elecCoordsStr = sprintf('%.1f ; %.1f', elec_coords(1), elec_coords(2));
    set(handles.elecCoordsStr, 'String', elecCoordsStr)
    
    % text(elec_coords(1), elec_coords(2)-0.05, 'V', 'HitTest','off','FontSize',24,...
    %     'Color','g','HorizontalAlignment','Center','VerticalAlignment','Bottom')
    ph = plot(elec_coords(1), elec_coords(2), 'x', 'HitTest','off',...
        'MarkerSize',30,'MarkerEdgeColor','g',...
        'LineWidth', 3);
    
    V = calc_field(src_vec, sig, elec_coords(1), elec_coords(2));
    
    elecValStr = sprintf('%.2f', V);
    set(handles.elecValStr, 'String', elecValStr)
else
    % force selected points to be outside source/sinks!
    cR = sqrt(sum(coordinates(1,1:2).*coordinates(1,1:2), 2));
    if cR < 3.
        coordinates(1,1:2) = coordinates(1,1:2)/cR * 3;
    end
    range_points = [range_points; coordinates(1,1:2)];
    ph = plot(coordinates(1,1), coordinates(1,2), 'd', 'HitTest','off',...
        'MarkerSize',10,'MarkerEdgeColor','m',...
        'LineWidth', 2);
    if size(range_points,1) > 1,
        plot(range_points(:,1), range_points(:,2),'md',...
            'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','m',...
            'LineStyle','-');
        figure;
        X = linspace(range_points(1,1), range_points(2,1), 50);
        Y = linspace(range_points(1,2), range_points(2,2), 50);
        V = calc_field(src_vec, sig, X, Y);
        d = sqrt(X.^2 + Y.^2);
        plot(d,V); xlabel('Distance from origin (a.u.)'); ylabel('Potential (mV)')
        range_points = [];
        set(handles.bSelRange,'Value',0)
        % assign some values for plotting to workspace variable
        [foo, iV_ex] = max(abs(V));
        extrema = struct('x_ex', d(iV_ex), 'V_ex', V(iV_ex), 'x_max', max(d));
        assignin('base','extrema',extrema)
    end
end

function obsDist_Callback(hObject, eventdata, handles)
% hObject    handle to obsDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of obsDist as text
%        str2double(get(hObject,'String')) returns contents of obsDist as a double
calcV_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function obsDist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to obsDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigmaStr_Callback(hObject, eventdata, handles)
% hObject    handle to sigmaStr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigmaStr as text
%        str2double(get(hObject,'String')) returns contents of sigmaStr as a double
global sig

sig = str2double(get(handles.sigmaStr, 'String'));

% --- Executes during object creation, after setting all properties.
function sigmaStr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigmaStr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in bLockColour.
function bLockColour_Callback(hObject, eventdata, handles)
% hObject    handle to bLockColour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bLockColour
global clim_locked

if get(handles.bLockColour,'Value') > 0
    clim_locked = get(handles.axes2,'CLim');
else
    clim_locked = -1;
end


% --- Executes on button press in bSelRange.
function bSelRange_Callback(hObject, eventdata, handles)
% hObject    handle to bSelRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bSelRange
global range_points 
range_points = [];
