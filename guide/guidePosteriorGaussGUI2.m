function varargout = guidePosteriorGaussGUI(varargin)
% guidePosteriorGaussGUI MATLAB code for guidePosteriorGaussGUI.fig
%      guidePosteriorGaussGUI, by itself, creates a new guidePosteriorGaussGUI or raises the existing
%      singleton*.
%
%      H = guidePosteriorGaussGUI returns the handle to a new guidePosteriorGaussGUI or the handle to
%      the existing singleton*.
%
%      guidePosteriorGaussGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in guidePosteriorGaussGUI.M with the given input arguments.
%
%      guidePosteriorGaussGUI('Property','Value',...) creates a new guidePosteriorGaussGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before guidePosteriorGaussGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to guidePosteriorGaussGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guidePosteriorGaussGUI

% Last Modified by GUIDE v2.5 11-Oct-2016 09:34:14

% Begin initialization code - DO NOT EDIT
addpath(utils.system.matlabPath('dynamicalSystems','guide'));
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guidePosteriorGaussGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @guidePosteriorGaussGUI_OutputFcn, ...
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


% --- Executes just before guidePosteriorGaussGUI is made visible.
function guidePosteriorGaussGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to guidePosteriorGaussGUI (see VARARGIN)

% Choose default command line output for guidePosteriorGaussGUI
handles.output = hObject;
dsobj = varargin{1};
assert(isa(dsobj, 'dynamicalSystem'), 'input is not a dynamicalSystems object!');
if ~isfield(dsobj.posterior, 'smooth')
    fprintf('(%s) --- Not all posteriors have been run for input. Running...\n', datestr(now, 'HH:MM:SS'));
    dsobj = dsobj.posteriorSmooth;
    fprintf('(%s) --- Complete...\n\n', datestr(now, 'HH:MM:SS'));
end

handles.dsObj = dsobj;

% series:
%    0 - Ground truth
%    1 - Filter
%    2 - Smooth
handles.showEmission = true;
handles.pSeries1     = 0;
handles.pSeries2     = 1;
doPlot(hObject, eventdata, handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes guidePosteriorGaussGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = guidePosteriorGaussGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in series1.
function series1_Callback(hObject, eventdata, handles)
% hObject    handle to series1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(hObject,'String'));
handles.pSeries1 = contents{get(hObject,'Value')} - 1;
doPlot(hObject, eventdata, handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function series1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to series1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in series2.
function series2_Callback(hObject, eventdata, handles)
% hObject    handle to series2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(hObject,'String'));
handles.pSeries2 = contents{get(hObject,'Value')} - 1;
doPlot(hObject, eventdata, handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function series2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to series2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxEmission.
function checkboxEmission_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxEmission (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxEmission
handles.showEmission = get(hObject,'Value');
doPlot(hObject, eventdata, handles);
guidata(hObject, handles);



function doPlot(hObject, eventdata, handles)
% main plot

col1 = [109,166,84]/255;
col2 = [160,105,190]/255;

clf(handles.axes1);
hold(handles.axes1, 'on');

if handles.showEmission
    plot(handles.axes1, handles.dsObj.y(1,:), handles.dsObj.y(2,:), 'k*');
end

if handles.pSeries1 == 0 || handles.pSeries2 == 0
    col = col1;
    if handles.pSeries2 == 0; col = col2; end
    plot(handles.axes1, handles.dsObj.y(1,:), handles.dsObj.y(2,:), '*-', 'Color', col);
end

if handles.pSeries1 == 1 || handles.pSeries2 == 1
    col    = col1;
    if handles.pSeries2 == 1; col = col2; end
    fMu    = handles.dsObj.posterior.filter.mu;
    fSigma = handles.dsObj.posterior.filter.sigma;
    for tt=1:handles.dsObj.T
        plot(handles.axes1, fMu(1,tt), fMu(2,tt), '+', 'Color', col);
        lc  = utils.plot.gaussian2DLevelCurve(1, fMu(:,tt), fSigma{tt}, 100);
        plot(handles.axes1, lc(:,1), lc(:,2), '-', 'Color', col);
    end
end

if handles.pSeries1 == 2 || handles.pSeries2 == 2
    col    = col1;
    if handles.pSeries2 == 0; col = col2; end
    fMu    = handles.dsObj.posterior.smooth.mu;
    fSigma = handles.dsObj.posterior.smooth.sigma;
    for tt=1:handles.dsObj.T
        plot(handles.axes1, fMu(1,tt), fMu(2,tt), '+', 'Color', col);
        lc  = utils.plot.gaussian2DLevelCurve(1, fMu(:,tt), fSigma{tt}, 100);
        plot(handles.axes1, lc(:,1), lc(:,2), '-', 'Color', col);
    end
end
hold(handles.axes1, 'off');