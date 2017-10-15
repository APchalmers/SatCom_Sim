%********************************************************************
% Wireless, Photonics and Space Engineering MSc. Programme
%
% Alvaro Perez Ortega
% Martin Anderberg 
% Han Zhou
% Chris J
%
% Chalmers University
%********************************************************************
% This work is licensed under a Creative Commons 
% Attribution-ShareAlike 4.0 International License.
% Info: https://creativecommons.org/licenses/by-sa/4.0/
%********************************************************************

function varargout = GUI(varargin)
% GUI MATLAB code for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 08-Oct-2017 10:24:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_OutputFcn, ...
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


% --- Executes just before GUI is made visible.
function GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI (see VARARGIN)

% Initialitation
%       Struct for GUI
set(handles.warning_text, 'String', 'LOADING...');
handles.nantennas = 1;
handles.theta = 1;
handles.azimuth = 0;
handles.roll = 0;
handles.sspointLon = 0;
handles.shape = 0;
handles.currentPopmenu = '1';
handles.error = 0;
set(handles.Theta_warning_text, 'String', '');
set(handles.Azimuth_warning_text, 'String', '');
set(handles.Roll_warning_text, 'String', '');
set(handles.sspointLon_warning_text, 'String', '');

% Default Run
%       Considered Spherical Map
[handles.Lat, handles.Lon] = intersecction(handles.theta,handles.azimuth,handles.roll);
handles.Lon = handles.Lon + handles.sspointLon*ones(1,size(handles.Lon,2));
handles.E = struct('Geometry','Polygon','Lon',handles.Lon','Lat',handles.Lat');
handles.map = worldmap([-90 90], [-180 180]);
handles.land = shaperead('landareas','UseGeoCoords',true);
geoshow(handles.land, 'FaceColor', [0.5 0.7 0.5]);
geoshow(handles.E,'FaceColor','none','DefaultEdgeColor','red','Linewidth',1.5);

set(handles.warning_text, 'String', 'READY');

% Choose default command line output for GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Final_pushbutton.
function Final_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Final_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Loading Frame
set(handles.warning_text, 'String', 'LOADING...');

% Conditions in case figure characteristics are changed
if handles.shape == 0
    cla reset;
    handles.map = worldmap([-90 90], [-180 180]);
    handles.land = shaperead('landareas','UseGeoCoords',true);
    geoshow(gca, handles.land, 'FaceColor', [0.5 0.7 0.5]);
else
    cla reset;
    geoshow('landareas.shp','FaceColor',[0.5 0.7 0.5]);   
    set(gca,'Xlim',[-180,180]);
    set(gca,'Ylim',[-90,90]);
    grid on;
    xlabel('Longitude [degrees]');
    ylabel('Lattitude [degrees]');
    set(gca,'XTick',[-180:20:180]);
    set(gca,'YTick',[-90:15:90]);
end

% Ploting Loop
for i = 1:handles.nantennas
    [handles.Lat, handles.Lon] = intersecction(handles.theta(i),handles.azimuth(i),handles.roll(i));
    handles.Lon = handles.Lon + handles.sspointLon*ones(1,size(handles.Lon,2));
    handles.E(i) = struct('Geometry','Polygon','Lon',handles.Lon','Lat',handles.Lat');
    geoshow(handles.E(i),'FaceColor','none','DefaultEdgeColor','red','Linewidth',1.5);
end

% Ready frame
set(handles.warning_text, 'String', 'READY');

guidata(hObject, handles);


% --- Executes on selection change in Figure_popmenu.
function Figure_popmenu_Callback(hObject, eventdata, handles)
% hObject    handle to Figure_popmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Case for Figure popmenu
str = get(hObject,'String');
val = get(hObject,'Value');
switch str{val}
    case 'Rectangular'
        handles.shape = 1;
    case 'Spherical'
        handles.shape = 0;
end
guidata(hObject,handles)
% Hints: contents = cellstr(get(hObject,'String')) returns Figure_popmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Figure_popmenu


% --- Executes during object creation, after setting all properties.
function Figure_popmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Figure_popmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function NumAntennas_edit_Callback(hObject, eventdata, handles)
% hObject    handle to NumAntennas_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumAntennas_edit as text
%        str2double(get(hObject,'String')) returns contents of NumAntennas_edit as a double


% --- Executes during object creation, after setting all properties.
function NumAntennas_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumAntennas_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Theta_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Theta_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Simply SAVING the edit value into the structure(NO BUTTON NEEDED)
val = get(handles.Antenna_popmenu,'Value');
handles.theta(val) = str2double(get(hObject,'String'));
guidata(hObject,handles);
% Hints: get(hObject,'String') returns contents of Theta_edit as text
%        str2double(get(hObject,'String')) returns contents of Theta_edit as a double


% --- Executes during object creation, after setting all properties.
function Theta_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Theta_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Antenna_popmenu.
function Antenna_popmenu_Callback(hObject, eventdata, handles)
% hObject    handle to Antenna_popmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Simply LOADING the values into the edit fields (NO BUTTON NEEDED)
val = get(hObject,'Value');
set(handles.Theta_edit,'String',num2str(handles.theta(val)));
set(handles.Azimuth_edit,'String',num2str(handles.azimuth(val)));
set(handles.Roll_edit,'String',num2str(handles.roll(val)));
guidata(hObject,handles);
% Hints: contents = cellstr(get(hObject,'String')) returns Antenna_popmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Antenna_popmenu


% --- Executes during object creation, after setting all properties.
function Antenna_popmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Antenna_popmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Azimuth_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Azimuth_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Simply SAVING the edit value into the structure(NO BUTTON NEEDED)
val = get(handles.Antenna_popmenu,'Value');
handles.azimuth(val) = str2double(get(hObject,'String'));
guidata(hObject,handles);
% Hints: get(hObject,'String') returns contents of Azimuth_edit as text
%        str2double(get(hObject,'String')) returns contents of Azimuth_edit as a double

% --- Executes during object creation, after setting all properties.
function Azimuth_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Azimuth_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Roll_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Roll_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Simply SAVING the edit value into the structure(NO BUTTON NEEDED)
val = get(handles.Antenna_popmenu,'Value');
handles.roll(val) = str2double(get(hObject,'String'));
guidata(hObject,handles);
% Hints: get(hObject,'String') returns contents of Roll_edit as text
%        str2double(get(hObject,'String')) returns contents of Roll_edit as a double

% --- Executes during object creation, after setting all properties.
function Roll_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Roll_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function sspointLon_edit_Callback(hObject, eventdata, handles)
% hObject    handle to sspointLon_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Simply SAVING the edit value into the structure(NO BUTTON NEEDED)
handles.sspointLon = str2double(get(hObject,'String'));
guidata(hObject,handles);
% Hints: get(hObject,'String') returns contents of sspointLon_edit as text
%        str2double(get(hObject,'String')) returns contents of sspointLon_edit as a double


% --- Executes during object creation, after setting all properties.
function sspointLon_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sspointLon_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Apply_button.
function Apply_button_Callback(hObject, eventdata, handles)
% hObject    handle to Apply_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Saves the number of antennas inputed
handles.nantennas = str2double(get(handles.NumAntennas_edit,'String'));

% In case that you want more antennas, it keeps the old values for the
% first ones
if (handles.nantennas < get(handles.Antenna_popmenu,'Value'))
    set(handles.Theta_edit,'String',num2str(handles.theta(handles.nantennas)));
    set(handles.Azimuth_edit,'String',num2str(handles.azimuth(handles.nantennas)));
    set(handles.Roll_edit,'String',num2str(handles.roll(handles.nantennas)));
    set(handles.Antenna_popmenu,'Value',handles.nantennas);
% Initializes the new antennas
elseif handles.nantennas > get(handles.Antenna_popmenu,'Value')
    for i=(get(handles.Antenna_popmenu,'Value')+1):handles.nantennas
        handles.theta(i) = 1;
        handles.azimuth(i) = 0;
        handles.roll(i) = 0;
    end
end

% Creates a dynamic popup menu for the number of antennas selected. This is
% done by creating a string with \n, which after used with sprintf creates
% a several line string required by the popup menu
handles.currentPopmenu = '';
for i=1:handles.nantennas
    handles.currentPopmenu = strcat(handles.currentPopmenu,num2str(i));
    if i ~= handles.nantennas
        handles.currentPopmenu = strcat(handles.currentPopmenu,' \n');
    end
end
handles.currentPopmenu = sprintf(handles.currentPopmenu);
set(handles.Antenna_popmenu,'String',handles.currentPopmenu);

% Loads the value of the last antenna directly, useful for checking
set(handles.Theta_edit,'String',num2str(handles.theta(handles.nantennas)));
set(handles.Azimuth_edit,'String',num2str(handles.azimuth(handles.nantennas)));
set(handles.Roll_edit,'String',num2str(handles.roll(handles.nantennas)));
set(handles.Antenna_popmenu,'Value',handles.nantennas);

guidata(hObject,handles);

function [lattitude, longitude] = intersecction(Theta_3dB,az,roll)
% Parameters (beware of the maybe weird cone steering)
R_E=earthRadius('kilometers');             % [km] Radius of Earth
H=35786;                % Altitude of satellite [km]
% Theta_3dB=4;            % 3 dB beamwidth in degrees
Cone_tilt=az;            % Cone tilt in degrees
Cone_rotation=roll;       % Cone rotation (clockwise) in degrees

% NOTE: Longitude reference is at satellite subpoint

% X-Y plane
rho=linspace(0,R_E,100);
Theta=linspace(0,2*pi);
[Theta_grid,rho_grid]=meshgrid(Theta,rho);
%[X_grid,Y_grid]=pol2cart(Theta_grid,rho_grid);
[X_grid,Y_grid]=meshgrid(linspace(-R_E,R_E));

% Move sphere with satellite as reference
rho_pos=R_E+H;
azimuth_pos=deg2rad(Cone_rotation);
elevation_pos=-deg2rad(90+Cone_tilt);
[x_pos,y_pos,z_pos]=sph2cart(azimuth_pos,elevation_pos,rho_pos);

% Sphere
Z_Sphere=-(z_pos+real(sqrt(R_E^2-(X_grid-x_pos).^2-(Y_grid-y_pos).^2)));

% Cone
r=tand(Theta_3dB/2)*(H+R_E);        % Radius of cone at centre of Earth
% Z_cone=sqrt(X_grid.^2+Y_grid.^2)/((r/(H+R_E))^2);
Z_cone=sqrt(X_grid.^2+Y_grid.^2)/tand(Theta_3dB/2);

% Find intersection
diffZ=Z_Sphere-Z_cone;
h = figure(4);
cut=contour(X_grid-x_pos,Y_grid+y_pos,diffZ,[0 0]);
close (h);

% Map the intersection onto sphere
Z_sphereIntersect=interp2(X_grid-x_pos,Y_grid+y_pos,Z_Sphere,cut(1,2:end),cut(2,2:end));

% Compute latitude and longitude

lattitude=rad2deg(asin(-cut(2,2:end)./(ones(1,length(cut(2,2:end)))*R_E)));
longitude=rad2deg(asin(cut(1,2:end)./sqrt(ones(1,length(cut(2,2:end)))*R_E.^2-cut(2,2:end).^2)));
