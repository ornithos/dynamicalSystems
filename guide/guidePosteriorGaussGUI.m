function varargout = guidePosteriorGaussGUI(varargin)
% GUIDEPOSTERIORGAUSSGUI MATLAB code for guidePosteriorGaussGUI.fig
%      GUIDEPOSTERIORGAUSSGUI, by itself, creates a new GUIDEPOSTERIORGAUSSGUI or raises the existing
%      singleton*.
%
%      H = GUIDEPOSTERIORGAUSSGUI returns the handle to a new GUIDEPOSTERIORGAUSSGUI or the handle to
%      the existing singleton*.
%
%      GUIDEPOSTERIORGAUSSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIDEPOSTERIORGAUSSGUI.M with the given input arguments.
%
%      GUIDEPOSTERIORGAUSSGUI('Property','Value',...) creates a new GUIDEPOSTERIORGAUSSGUI or raises the
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

% Last Modified by GUIDE v2.5 11-Oct-2016 16:57:30

% Begin initialization code - DO NOT EDIT
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
if isempty(varargin)
    guidata(hObject, handles);
    return
end

dsobj = varargin{1};
assert(isa(dsobj, 'dynamicalSystem'), 'input is not a dynamicalSystems object!');
sp1   = varargin{2};
sp2   = varargin{3};

handles.dsObj = dsobj;
handles.sp    = {sp1, sp2};
% series:
%    0 - Ground truth
%    1 - Filter
%    2 - Smooth
handles.showEmission = true;
handles.pSeries1     = 0;
handles.sp2_1        = 0;
handles.pSeries2     = 1;
handles.sp2_2        = 0;

handles.doLC1        = 1;
handles.doLC2        = 1;
handles.lcAlpha      = 1;   % standard dev. of level curve plotted for Gaussian

handles.popupS1.String = {'Ground Truth', [sp1.descr, '::filter'], [sp1.descr, '::smooth'], ...
    [sp2.descr, '::filter'], [sp2.descr, '::smooth']};
handles.popupS2.String = {'Ground Truth', [sp1.descr, '::filter'], [sp1.descr, '::smooth'], ...
    [sp2.descr, '::filter'], [sp2.descr, '::smooth']};
handles.popupS2.Value = 1 + (handles.pSeries2 + 0) + handles.sp2_1*2;

%% TRANSFORM LATENT INTO OBSERVATION SPACE
handles.yhat = NaN(size(dsobj.y));
handles.sp{1}.infer.filter.yhat = NaN(size(dsobj.y));
handles.sp{2}.infer.filter.yhat = NaN(size(dsobj.y));
handles.sp{1}.infer.smooth.yhat = NaN(size(dsobj.y));
handles.sp{2}.infer.smooth.yhat = NaN(size(dsobj.y));

handles.sp{1}.infer.filter.yhatSigma = cell(size(dsobj.y, 2),1);
handles.sp{2}.infer.filter.yhatSigma = cell(size(dsobj.y, 2),1);
handles.sp{1}.infer.smooth.yhatSigma = cell(size(dsobj.y, 2),1);
handles.sp{2}.infer.smooth.yhatSigma = cell(size(dsobj.y, 2),1);

[~, ~, h, Dh] = dsobj.functionInterfaces;

for tt = 1:dsobj.d.T
    if dsobj.emiLinear
        % transform posterior mean into observation space
        handles.yhat(:,tt) = dsobj.par.H * dsobj.x(:,tt);
        handles.sp{1}.infer.filter.yhat(:,tt) = dsobj.par.H * handles.sp{1}.infer.filter.mu(:,tt);
        handles.sp{2}.infer.filter.yhat(:,tt) = dsobj.par.H * handles.sp{2}.infer.filter.mu(:,tt);
        handles.sp{1}.infer.smooth.yhat(:,tt) = dsobj.par.H * handles.sp{1}.infer.smooth.mu(:,tt);
        handles.sp{2}.infer.smooth.yhat(:,tt) = dsobj.par.H * handles.sp{2}.infer.smooth.mu(:,tt);
        
        % transformed covariances into observation space
        handles.sp{1}.infer.filter.yhatSigma{tt} = dsobj.par.H * handles.sp{1}.infer.filter.sigma{tt} * dsobj.par.H' + dsobj.par.R;
        handles.sp{2}.infer.filter.yhatSigma{tt} = dsobj.par.H * handles.sp{2}.infer.filter.sigma{tt} * dsobj.par.H' + dsobj.par.R;
        handles.sp{1}.infer.smooth.yhatSigma{tt} = dsobj.par.H * handles.sp{1}.infer.smooth.sigma{tt} * dsobj.par.H' + dsobj.par.R;
        handles.sp{2}.infer.smooth.yhatSigma{tt} = dsobj.par.H * handles.sp{2}.infer.smooth.sigma{tt} * dsobj.par.H' + dsobj.par.R;
    else
        % transform posterior mean into observation space
        handles.yhat(:,tt) = dsobj.par.h(dsobj.x(:,tt));
        handles.sp{1}.infer.filter.yhat(:,tt) = h(handles.sp{1}.infer.filter.mu(:,tt));
        handles.sp{2}.infer.filter.yhat(:,tt) = h(handles.sp{2}.infer.filter.mu(:,tt));
        handles.sp{1}.infer.smooth.yhat(:,tt) = h(handles.sp{1}.infer.smooth.mu(:,tt));
        handles.sp{2}.infer.smooth.yhat(:,tt) = h(handles.sp{2}.infer.smooth.mu(:,tt));
        
        % transformed covariances into observation space
        H               = Dh(dsobj.par.h(handles.sp{1}.infer.filter.mu(:,tt)));
        handles.sp{1}.infer.filter.yhatSigma{tt} = H * handles.sp{1}.infer.filter.sigma{tt} * H' + dsobj.par.R;
        H               = Dh(dsobj.par.h(handles.sp{2}.infer.filter.mu(:,tt)));
        handles.sp{2}.infer.filter.yhatSigma{tt} = H * handles.sp{2}.infer.filter.sigma{tt} * H' + dsobj.par.R;
        H               = Dh(dsobj.par.h(handles.sp{1}.infer.smooth.mu(:,tt)));
        handles.sp{1}.infer.smooth.yhatSigma{tt} = H * handles.sp{1}.infer.smooth.sigma{tt} * H' + dsobj.par.R;
        H               = Dh(dsobj.par.h(handles.sp{2}.infer.smooth.mu(:,tt)));
        handles.sp{2}.infer.smooth.yhatSigma{tt} = H * handles.sp{2}.infer.smooth.sigma{tt} * H' + dsobj.par.R;
    end
end

%% Cut down to 2D for visualisation
handles.y           = dsobj.y;

if dsobj.d.y > 2
    if dsobj.opts.warning; warning('Output space is not two dimensional. Only the first two will be used.'); end
    handles.y           = handles.y(1:2,:);
    handles.yhat        = handles.yhat(1:2,:);
    handles.sp{1}.infer.filter.yhat = handles.sp{1}.infer.filter.yhat(1:2,:);
    handles.sp{2}.infer.filter.yhat = handles.sp{2}.infer.filter.yhat(1:2,:);
    handles.sp{1}.infer.smooth.yhat = handles.sp{1}.infer.smooth.yhat(1:2,:);
    handles.sp{2}.infer.smooth.yhat = handles.sp{2}.infer.smooth.yhat(1:2,:);
    for tt = 1:dsobj.d.T
        handles.sp{1}.infer.filter.yhatSigma{tt} = handles.sp{1}.infer.filter.yhatSigma{tt}(1:2,1:2);
        handles.sp{2}.infer.filter.yhatSigma{tt} = handles.sp{2}.infer.filter.yhatSigma{tt}(1:2,1:2);
        handles.sp{1}.infer.smooth.yhatSigma{tt} = handles.sp{1}.infer.smooth.yhatSigma{tt}(1:2,1:2);
        handles.sp{2}.infer.smooth.yhatSigma{tt} = handles.sp{2}.infer.smooth.yhatSigma{tt}(1:2,1:2);
    end
    
elseif dsobj.d.y == 1
    if dsobj.opts.warnings; warning('Output space is not two dimensional. A first ''time'' dimension will be added.'); end
    handles.y          = [1:dsobj.d.T; handles.y];
    handles.yhat       = [1:dsobj.d.T; handles.yhat];
    handles.sp{1}.infer.filter.yhat = [1:dsobj.d.T; handles.sp{1}.infer.filter.yhat];
    handles.sp{2}.infer.filter.yhat = [1:dsobj.d.T; handles.sp{2}.infer.filter.yhat];
    handles.sp{1}.infer.smooth.yhat = [1:dsobj.d.T; handles.sp{1}.infer.smooth.yhat];
    handles.sp{2}.infer.smooth.yhat = [1:dsobj.d.T; handles.sp{2}.infer.smooth.yhat];
    for tt = 1:dsobj.d.T
        handles.sp{1}.infer.filter.yhatSigma{tt} = [1e-4, 0; 0, handles.sp{1}.infer.filter.yhatSigma{tt}];
        handles.sp{2}.infer.filter.yhatSigma{tt} = [1e-4, 0; 0, handles.sp{2}.infer.filter.yhatSigma{tt}];
        handles.sp{1}.infer.smooth.yhatSigma{tt} = [1e-4, 0; 0, handles.sp{1}.infer.smooth.yhatSigma{tt}];
        handles.sp{2}.infer.smooth.yhatSigma{tt} = [1e-4, 0; 0, handles.sp{2}.infer.smooth.yhatSigma{tt}];
    end
end

%% finit
doPlot(hObject, eventdata, handles);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes guidePosteriorGaussGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% --- Executes during object creation, after setting all properties.
function popupS1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupS2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
    set(hObject, 'Value', handles.pSeries2 + 1 + handles.sp2_2*2);
end

function editLvlSD_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Outputs from this function are returned to the command line.
function varargout = guidePosteriorGaussGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupS1.
function popupS1_Callback(hObject, eventdata, handles)
% hObject    handle to popupS1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(hObject,'String'));
inpVal   = get(hObject,'Value');
if inpVal > 3
    handles.sp2_1 = 1;
    inpVal = inpVal - 2;
else
    handles.sp2_1 = 0;
end
handles.pSeries1 = inpVal - 1;
doPlot(hObject, eventdata, handles);
guidata(hObject, handles);



% --- Executes on selection change in popupS2.
function popupS2_Callback(hObject, eventdata, handles)
% hObject    handle to popupS2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupS2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupS2
contents = cellstr(get(hObject,'String'));
inpVal   = get(hObject,'Value');
if inpVal > 3
    handles.sp2_2 = 1;
    inpVal = inpVal - 2;
else
    handles.sp2_2 = 0;
end
handles.pSeries2 = inpVal - 1;
doPlot(hObject, eventdata, handles);
guidata(hObject, handles);



% --- Executes on button press in checkboxObs.
function checkboxObs_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxObs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxObs
handles.showEmission = get(hObject,'Value');
doPlot(hObject, eventdata, handles);
guidata(hObject, handles);

function editLvlSD_Callback(hObject, eventdata, handles)
% hObject    handle to editLvlSD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editLvlSD as text
%        str2double(get(hObject,'String')) returns contents of editLvlSD as a double
editval         = get(hObject,'String');
testval         = regexp(editval, '^[0-9]+(\.[0-9]+)?$');
if(iscell(testval)); testval = testval{1}; end
assert(~isempty(1), 'Input sigma value must be a positive real (rational) number');
handles.lcAlpha = str2double(editval);
doPlot(hObject, eventdata, handles);
guidata(hObject, handles);


function doPlot(hObject, eventdata, handles)
% main plot

col1 = [109,166,84]/255;
col2 = [160,105,190]/255;

cla(handles.axes1);
hold(handles.axes1, 'on');

% observations
if handles.showEmission
    plot(handles.axes1, handles.y(1,:), handles.y(2,:), 'k*');
    plot(handles.axes1, handles.y(1,:), handles.y(2,:), 'c:');
end

% ground truth
%if handles.pSeries1 == 0 || handles.pSeries2 == 0
%    col = col1;
%    if handles.pSeries2 == 0; col = col2; end
%    plot(handles.axes1, handles.yhat(1,:), handles.yhat(2,:), '*-', 'Color', col);
%end

% other plots
doPlotInner(handles.pSeries1, 1-handles.sp2_1, handles, col1, handles.doLC1);
doPlotInner(handles.pSeries2, 1-handles.sp2_2, handles, col2, handles.doLC2);

hold(handles.axes1, 'off');



function doPlotInner(type, doSP1, handles, col, doLC)

    % ground truth
    if type == 0 
        plot(handles.axes1, handles.yhat(1,:), handles.yhat(2,:), '*-', 'Color', col);
        return
    end
    
    if type == 1 && doSP1
        mu      = handles.sp{1}.infer.filter.yhat;
        sigma   = handles.sp{1}.infer.filter.yhatSigma;
    elseif type == 1 && ~doSP1
        mu      = handles.sp{2}.infer.filter.yhat;
        sigma   = handles.sp{2}.infer.filter.yhatSigma;
    elseif type == 2 && doSP1
        mu      = handles.sp{1}.infer.smooth.yhat;
        sigma   = handles.sp{1}.infer.smooth.yhatSigma;
    elseif type == 2 && ~doSP1
        mu      = handles.sp{2}.infer.smooth.yhat;
        sigma   = handles.sp{2}.infer.smooth.yhatSigma;
    end
    
    if doLC
        for tt=1:handles.dsObj.d.T
            plot(handles.axes1, mu(1,tt), mu(2,tt), '+', 'Color', col);
            lc  = utils.plot.gaussian2DLevelCurve(handles.lcAlpha, mu(:,tt), sigma{tt}, 100);
            plot(handles.axes1, lc(:,1), lc(:,2), '-', 'Color', col);
        end
        plot(handles.axes1, mu(1,:), mu(2,:), '-', 'Color', col);
    else
        plot(handles.axes1, mu(1,:), mu(2,:), '+-', 'Color', col);
    end
    
% --- Executes on button press in checkLC1.
function checkLC1_Callback(hObject, eventdata, handles)
% hObject    handle to checkLC1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.doLC1   = get(hObject,'Value');
doPlot(hObject, eventdata, handles);
guidata(hObject, handles);

% --- Executes on button press in checkLC2.
function checkLC2_Callback(hObject, eventdata, handles)
% hObject    handle to checkLC2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.doLC2   = get(hObject,'Value');
doPlot(hObject, eventdata, handles);
guidata(hObject, handles);


% --- Executes on key press with focus on popupS2 and none of its controls.
function popupS2_KeyPressFcn(hObject, eventdata, handles)
%eventdata % Let's see the KeyPress event data
switch eventdata.Key
    case 'f'
        handles.pSeries2 = 1;
    case 'g'
        handles.pSeries2 = 0;
    case 's'
        handles.pSeries2 = 2;
    case 'e'
        handles.sp2_2 = 1;
    case 'd'
        handles.sp2_2 = 0;
end
set(hObject, 'Value', handles.pSeries2 + 1 + handles.sp2_2*2);
doPlot(hObject, eventdata, handles);
guidata(hObject, handles);
