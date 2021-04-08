% SASICA() - perform automatic detection of artifactual components
%
% Usage:
%   >> [EEG com] = SASICA();
%       Raises a GUI
%       Pressing "Compute" uses the selected methods to detect artifactual
%       components.
%       Available methods are:
%              Autocorrelation: detects noisy components with weak
%                               autocorrelation (muscle artifacts usually)
%              Focal components: detects components that are too focal and
%                               thus unlikely to correspond to neural
%                               activity (bad channel or muscle usually).
%              Focal trial activity: detects components with focal trial
%                               activity, with same algorhithm as focal
%                               components above. Results similar to trial
%                               variability.
%              Signal to noise ratio: detects components with weak signal
%                               to noise ratio between arbitrary baseline
%                               and interest time windows.
%              Dipole fit residual variance: detects components with high
%                               residual variance after subtraction of the
%                               forward dipole model. Note that the inverse
%                               dipole modeling using DIPFIT2 in EEGLAB
%                               must have been computed to use this
%                               measure.
%              EOG correlation: detects components whose time course
%                               correlates with EOG channels.
%              Bad channel correlation: detects components whose time course
%                               correlates with any channel(s).
%              ADJUST selection: use ADJUST routines to select components
%                               (see Mognon, A., Jovicich, J., Bruzzone,
%                               L., & Buiatti, M. (2011). ADJUST: An
%                               automatic EEG artifact detector based on
%                               the joint use of spatial and temporal
%                               features. Psychophysiology, 48(2), 229Ã¢ÂÂ240.
%                               doi:10.1111/j.1469-8986.2010.01061.x)
%              FASTER selection: use FASTER routines to select components
%                               (see Nolan, H., Whelan, R., & Reilly, R. B.
%                               (2010). FASTER: Fully Automated Statistical
%                               Thresholding for EEG artifact Rejection.
%                               Journal of Neuroscience Methods, 192(1),
%                               152Ã¢ÂÂ162. doi:16/j.jneumeth.2010.07.015)
%              MARA selection:  use MARA classification engine to select components
%                               (see Winkler I, Haufe S, Tangermann M.
%                               2011. Automatic Classification of
%                               Artifactual ICA-Components for Artifact
%                               Removal in EEG Signals. Behavioral and
%                               Brain Functions. 7:30.)
%
%   >> [EEG com] = SASICA( [], 'key', 'val');
%       Takes optional key val pairs to set specific options from the
%       command line.
%   >> SASICA('resetprefs')
%       Resets all options to default values prior to opening the gui.
%   >> cfg = SASICA('getdefs')
%       Returns default options in a structure cfg, which can be used as
%       input to eeg_SASICA
%
% Inputs:
%   EEG       - input EEG structure
%
% Optional inputs: are generated when you press OK in the GUI.
%                  Please explore from the LASTCOM variable in your workspace.
%
%
% Outputs:
%   EEG          - EEGLAB data structure
%   com          - [str] command generating the same selection of
%                   components
%
% If you use this program in your research, please cite the following
% article: 
%   Chaumon M, Bishop DV, Busch NA. A Practical Guide to the Selection of
%   Independent Components of the Electroencephalogram for Artifact
%   Correction. Journal of neuroscience methods. 2015  
%
%   SASICA is a software that helps select independent components of
%   the electroencephalogram based on various signal measures.
%     Copyright (C) 2014  Maximilien Chaumon
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


function varargout = SASICA(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SASICA_OpeningFcn, ...
                   'gui_OutputFcn',  @SASICA_OutputFcn, ...
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


% --- Executes just before SASICA is made visible.
function SASICA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SASICA (see VARARGIN)


% first deal with getting EEG
try
    if ~isempty(varargin) && isstruct(varargin{1}) && isfield(varargin{1},'setname')
        EEG = varargin{1};
        varargin(1) = [];
    else
        EEG = evalin('base','EEG;');
        if ~isempty(varargin) && isempty(varargin{1})
            varargin(1) = [];
        end
    end
catch
    EEG = pop_loadset;
end
if isempty(EEG)
    errordlg('No EEG loaded')
    return
end
if isempty([EEG.chanlocs.X])
    errordlg('No electrode locations provided.')
    error('No electrode locations provided.')
end

assignin('base','EEG',EEG);
EEG.times = 1000*(EEG.xmin:1/EEG.srate:EEG.xmax);
if isempty(EEG.icawinv)
    errordlg('No ica weights in the current EEG dataset! Compute ICA on your data first.')
    error('No ica weights! Compute ICA on your data first.')
end
try
    EEG.icaact = eeg_getica(EEG);
    % note here we overwrite any EEG structure in the base workspace.
    assignin('base','EEG',EEG);
end
handles.EEG = EEG;
% then with cfg
if numel(varargin) == 1 && isstruct(varargin{1})
    % assume we have only cfg as input
    cfg = varargin{1};
    cfg = setdef(cfg,getpref('SASICA','cfg',getdefs));
elseif numel(varargin) > 1 && ischar(varargin{1})
    % assume we have param, argument pairs
    cfg = vararg2struct(varargin,'_');
    cfg = setdef(cfg,getpref('SASICA','cfg',getdefs));
    % run eeg_SASICA straight
    EEG = eeg_SASICA(EEG,cfg);
    setpref('SASICA','cfg',cfg);
    handles.EEG = EEG;
    return;
else
    cfg = getpref('SASICA','cfg',getdefs);
end
% create command line for eegh
if not(isempty(cfg))
    handles.com = ['[EEG] = SASICA([],' vararg2str(struct2vararg(cfg)) ');'];
else
    handles.com = '[EEG] = SASICA();';
end
guidata(hObject,handles);

if not(isempty(cfg))
    fs = fieldnames(handles);
    fedits = fs(regexpcell(fs,'edit_'));
    fchecks = fs(regexpcell(fs,'check_'));
    f = regexp(strrep(fchecks,'check_',''),'[^_]*','match');
    for i_f = 1:numel(fchecks)
        if not(isfield(cfg,f{i_f}{1}))
            continue
        end
        set(handles.(fchecks{i_f}),'Value',cfg.(f{i_f}{1}).(f{i_f}{2}));
        check_enable(handles.(fchecks{i_f}),handles)
    end
    f = regexp(strrep(fedits,'edit_',''),'[^_]*','match');
    for i_f = 1:numel(fedits)
        if not(isempty(strfind(fedits{i_f},'channame')))
            chname = channamecheck(cfg.(f{i_f}{1}).(f{i_f}{2}),hObject);
            set(handles.(fedits{i_f}),'String',chname);
        else
            set(handles.(fedits{i_f}),'String',num2str(cfg.(f{i_f}{1}).(f{i_f}{2})));
        end
    end
end

if not(ndims(EEG.icaact) == 3)
    set(handles.check_trialfoc_enable,'value',0);
    check_enable(handles.check_trialfoc_enable,handles);
end
if not(isfield(EEG,'dipfit')) || isempty(EEG.dipfit)
    set(handles.check_resvar_enable,'value',0);
    check_enable(handles.check_resvar_enable,handles);
end


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SASICA wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SASICA_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if isfield(handles,'EEG')
    varargout = {handles.EEG};
    varargout{2} = handles.com;
else
    try
        varargout = {evalin('base','EEG') fastif(isfield(handles,'com'),handles.com,[])};
    catch
        varargout = {[] []};
    end
end



function edit_autocorr_dropautocorr_Callback(hObject, eventdata, handles)
% hObject    handle to edit_autocorr_dropautocorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_autocorr_dropautocorr as text
%        str2double(get(hObject,'String')) returns contents of edit_autocorr_dropautocorr as a double
setThreshString(hObject)

function setThreshString(hObject)
dat = regexp(get(hObject,'String'),'(auto)(.*)','tokens');
if not(isempty(dat)) && not(isempty(strfind(dat{1}{1},'auto')))
    if isempty(strtrim(dat{1}{2}))
        dat{1}{2} = '2';
    end
    dat = ['auto ' num2str(str2double(dat{1}{2}))];
else
    dat = num2str(str2double(get(hObject,'String')),'%g');
end
set(hObject,'String',dat)


% --- Executes during object creation, after setting all properties.
function edit_autocorr_dropautocorr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_autocorr_dropautocorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_autocorr_autocorrint_Callback(hObject, eventdata, handles)
% hObject    handle to edit_autocorr_autocorrint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_autocorr_autocorrint as text
%        str2double(get(hObject,'String')) returns contents of edit_autocorr_autocorrint as a double
setThreshString(hObject)


% --- Executes during object creation, after setting all properties.
function edit_autocorr_autocorrint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_autocorr_autocorrint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_focalcomp_focalICAout_Callback(hObject, eventdata, handles)
% hObject    handle to edit_focalcomp_focalICAout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_focalcomp_focalICAout as text
%        str2double(get(hObject,'String')) returns contents of edit_focalcomp_focalICAout as a double
setThreshString(hObject)


% --- Executes during object creation, after setting all properties.
function edit_focalcomp_focalICAout_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_focalcomp_focalICAout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_SNR_snrPOI_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SNR_snrPOI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SNR_snrPOI as text
%        str2double(get(hObject,'String')) returns contents of edit_SNR_snrPOI as a double
set(hObject,'String',num2str(str2num(get(hObject,'String')),'%g '))


% --- Executes during object creation, after setting all properties.
function edit_SNR_snrPOI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SNR_snrPOI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_SNR_snrBL_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SNR_snrBL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SNR_snrBL as text
%        str2double(get(hObject,'String')) returns contents of edit_SNR_snrBL as a double
set(hObject,'String',num2str(str2num(get(hObject,'String')),'%g '))


% --- Executes during object creation, after setting all properties.
function edit_SNR_snrBL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SNR_snrBL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_SNR_snrcut_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SNR_snrcut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SNR_snrcut as text
%        str2double(get(hObject,'String')) returns contents of edit_SNR_snrcut as a double
set(hObject,'String',num2str(str2num(get(hObject,'String')),'%g'))


% --- Executes during object creation, after setting all properties.
function edit_SNR_snrcut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SNR_snrcut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_EOGcorr_corthreshV_Callback(hObject, eventdata, handles)
% hObject    handle to edit_EOGcorr_corthreshV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_EOGcorr_corthreshV as text
%        str2double(get(hObject,'String')) returns contents of edit_EOGcorr_corthreshV as a double
setThreshString(hObject)


% --- Executes during object creation, after setting all properties.
function edit_EOGcorr_corthreshV_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_EOGcorr_corthreshV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_EOGcorr_Veogchannames_Callback(hObject, eventdata, handles)
% hObject    handle to edit_EOGcorr_Veogchannames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_EOGcorr_Veogchannames as text
%        str2double(get(hObject,'String')) returns contents of edit_EOGcorr_Veogchannames as a double
numchan = get(hObject,'String');
chname = channamecheck(numchan,hObject);
set(hObject,'String',chname);

% --- Executes during object creation, after setting all properties.
function edit_EOGcorr_Veogchannames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_EOGcorr_Veogchannames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_EOGcorr_corthreshH_Callback(hObject, eventdata, handles)
% hObject    handle to edit_EOGcorr_corthreshH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_EOGcorr_corthreshH as text
%        str2double(get(hObject,'String')) returns contents of edit_EOGcorr_corthreshH as a double
setThreshString(hObject)


% --- Executes during object creation, after setting all properties.
function edit_EOGcorr_corthreshH_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_EOGcorr_corthreshH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_EOGcorr_Heogchannames_Callback(hObject, eventdata, handles)
% hObject    handle to edit_EOGcorr_Heogchannames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_EOGcorr_Heogchannames as text
%        str2double(get(hObject,'String')) returns contents of edit_EOGcorr_Heogchannames as a double
numchan = get(hObject,'String');
chname = channamecheck(numchan,hObject);
set(hObject,'String',chname);


% --- Executes during object creation, after setting all properties.
function edit_EOGcorr_Heogchannames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_EOGcorr_Heogchannames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_chancorr_corthresh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_chancorr_corthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_chancorr_corthresh as text
%        str2double(get(hObject,'String')) returns contents of edit_chancorr_corthresh as a double
setThreshString(hObject)


% --- Executes during object creation, after setting all properties.
function edit_chancorr_corthresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_chancorr_corthresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_chancorr_channames_Callback(hObject, eventdata, handles)
% hObject    handle to edit_chancorr_channames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_chancorr_channames as text
%        str2double(get(hObject,'String')) returns contents of edit_chancorr_channames as a double
numchan = get(hObject,'String');
chname = channamecheck(numchan,hObject);
set(hObject,'String',chname);


% --- Executes during object creation, after setting all properties.
function edit_chancorr_channames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_chancorr_channames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_help.
function push_help_Callback(hObject, eventdata, handles)
% hObject    handle to push_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pophelp('SASICA');

% --- Executes on button press in push_close.
function push_close_Callback(hObject, eventdata, handles)
% hObject    handle to push_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(gcf);

% --- Executes on button press in push_ok.
function push_ok_Callback(hObject, eventdata, handles)
% hObject    handle to push_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    handles.EEG = evalin('base','EEG');
    if isempty(handles.EEG)
        error;
    end
catch
    errordlg('Please load a proper EEG structure first');
    error('Please load a proper EEG structure first');
end
EEG = handles.EEG;
fs = fieldnames(handles);
fedits = fs(regexpcell(fs,'edit_'));
fchecks = fs(regexpcell(fs,'check_'));
f = regexp(strrep(fchecks,'check_',''),'[^_]*','match');
for i_f = 1:numel(fchecks)
    cfg.(f{i_f}{1}).(f{i_f}{2}) = eval(['[' num2str(get(handles.(fchecks{i_f}),'Value')) ']']);
end
f = regexp(strrep(fedits,'edit_',''),'[^_]*','match');
for i_f = 1:numel(fedits)
    if not(isempty(strfind(fedits{i_f},'channame')))
        chan = get(handles.(fedits{i_f}),'String');
        cfg.(f{i_f}{1}).(f{i_f}{2}) = channamecheck(chan,hObject);
    else
        try % if it's a number this is no error
            cfg.(f{i_f}{1}).(f{i_f}{2}) = eval(['[' get(handles.(fedits{i_f}),'String') ']']);
        catch % if it's a string (auto #) then we can't evaluate
            cfg.(f{i_f}{1}).(f{i_f}{2}) = get(handles.(fedits{i_f}),'String');
        end

    end
end
cfg.opts.noplot = get(handles.checkNoPlot,'Value');
cfg.opts.nocompute = get(handles.checkReplot,'Value');
cfg = setdef(cfg,getpref('SASICA','cfg',getdefs));

try
    [EEG, cfg] = eeg_SASICA(EEG,cfg);
catch ME
    disp('================================');
    disp('================================');
    disp('ERROR. Please send the entire error log below to max.chaumon@gmail.com. Thanks for your help!');
    disp('================================');
    disp(['This is ' eegplugin_SASICA])
    disp(['This is MATLAB ' version])
    disp(['Running on ' computer])
    dispstruct(cfg)
    rethrow(ME)
end
setpref('SASICA','cfg',cfg);

com = ['[EEG] = SASICA(EEG,' vararg2str(struct2vararg(cfg)) ');'];
assignin('base','EEG',EEG);
assignin('base','LASTCOM',com);
% close(handles.figure1);

function dispstruct(S,space)
if nargin == 1
    space = 0;
end
f = fieldnames(S);
for i_f = 1:numel(f)
    if ~isstruct(S.(f{i_f}))
        if isnumeric(S.(f{i_f}))
            disp([repmat('  ',space), f{i_f} ': ', num2str(S.(f{i_f}))])
        else
            disp([repmat('  ',space), f{i_f} ': ', S.(f{i_f})])
        end
    else
        disp([repmat('  ',space), f{i_f} ':'])
        dispstruct(S.(f{i_f}),space+1)
    end
end

function edit_trialfoc_focaltrialout_Callback(hObject, eventdata, handles)
% hObject    handle to edit_trialfoc_focaltrialout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_trialfoc_focaltrialout as text
%        str2double(get(hObject,'String')) returns contents of edit_trialfoc_focaltrialout as a double
setThreshString(hObject)


% --- Executes during object creation, after setting all properties.
function edit_trialfoc_focaltrialout_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_trialfoc_focaltrialout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in check_autocorr_enable.
function check_autocorr_enable_Callback(hObject, eventdata, handles)
% hObject    handle to check_autocorr_enable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_autocorr_enable
check_enable(hObject,handles)


% --- Executes on button press in check_focalcomp_enable.
function check_focalcomp_enable_Callback(hObject, eventdata, handles)
% hObject    handle to check_focalcomp_enable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_focalcomp_enable
check_enable(hObject,handles)


% --- Executes on button press in check_trialfoc_enable.
function check_trialfoc_enable_Callback(hObject, eventdata, handles)
% hObject    handle to check_trialfoc_enable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_trialfoc_enable
try
    handles.EEG = evalin('base','EEG');
end
if ndims(handles.EEG.icaact) < 3
    errordlg({'EEG contains continuous data (1 trial). Focal trial activity does not make sense...'})
    set(hObject,'value',0)
end
check_enable(hObject,handles)


% --- Executes on button press in check_SNR_enable.
function check_SNR_enable_Callback(hObject, eventdata, handles)
% hObject    handle to check_SNR_enable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_SNR_enable
check_enable(hObject,handles)


% --- Executes on button press in check_EOGcorr_enable.
function check_EOGcorr_enable_Callback(hObject, eventdata, handles)
% hObject    handle to check_EOGcorr_enable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_EOGcorr_enable
check_enable(hObject,handles)

% --- Executes on button press in check_chancorr_enable.
function check_chancorr_enable_Callback(hObject, eventdata, handles)
% hObject    handle to check_chancorr_enable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_chancorr_enable
check_enable(hObject,handles)

% --- Executes on button press in check_resvar_enable.
function check_resvar_enable_Callback(hObject, eventdata, handles)
% hObject    handle to check_resvar_enable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_resvar_enable

try
    handles.EEG = evalin('base','EEG');
end
if ~(isfield(handles.EEG,'dipfit') && not(isempty(handles.EEG.dipfit)))
    errordlg({'You must have dipole fits computed in your dataset in order to use this feature.'
        'See tools > locate dipoles with dipfit in EEGLAB'})
    set(hObject,'value',0)
end
check_enable(hObject,handles)

function check_enable(hObject,handles)
fs = fieldnames(handles);
theones = regexp(get(hObject,'Tag'),'check_([^_]+)','tokens');
theones = theones{1}{1};
fedits = fs(regexpcell(fs,['edit_' theones]));
for i_f = 1:numel(fedits)
    if get(hObject,'value')
       set(handles.(fedits{i_f}),'enable','on');
    else
       set(handles.(fedits{i_f}),'enable','off');
    end
end




function v = struct2vararg(s,tag)

if not(exist('tag','var'))
    tag = '';
end
if numel(s) ~= 1
    error('input should be 1x1 structure array')
end

fn = fieldnames(s);
v = {};
for i_f = 1:numel(fn)
    if isstruct(s.(fn{i_f}));
        if not(isempty(tag))
            v = [v struct2vararg(s.(fn{i_f}),[tag '_' fn{i_f}])];
        else
            v = [v struct2vararg(s.(fn{i_f}),[fn{i_f}])];
        end
        continue
    end
    if not(isempty(tag))
        v{end+1} = [tag '_' fn{i_f}];
    else
        v{end+1} = [fn{i_f}];
    end
    v{end+1} = s.(fn{i_f});
end


function s = vararg2struct(v,tag)

if not(exist('tag','var'))
    tag = '_';
end

f = regexp(v(1:2:end),['[^'  regexptranslate('escape',tag) ']*'],'match');
for i_f = 1:numel(f)
    str = 's';
    for i_ff = 1:numel(f{i_f})
        str = [str '.' f{i_f}{i_ff}];
    end
    str = [str ' = v{i_f*2};'];
    eval(str);
end

% --------------------------------------------------------------------
function Opts_Callback(hObject, eventdata, handles)
% hObject    handle to Opts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function resetprefs_Callback(hObject, eventdata, handles)
% hObject    handle to resetprefs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

def = getdefs;
setpref('SASICA','cfg',def);
SASICA(def)

function resetprefs()
resetprefs_Callback()

function def = getdefs

def.autocorr.enable = false;
def.autocorr.dropautocorr = 'auto';
def.autocorr.autocorrint = 20;% will compute autocorrelation with this many milliseconds lag

def.focalcomp.enable = false;
def.focalcomp.focalICAout = 'auto';

def.trialfoc.enable = false;
def.trialfoc.focaltrialout = 'auto';

def.resvar.enable = false;
def.resvar.thresh = 15;% %residual variance allowed

def.SNR.enable = false;
def.SNR.snrPOI = [0 Inf];% period of interest (signal)
def.SNR.snrBL = [-Inf 0];% period of no interest (noise)
def.SNR.snrcut = 1;% SNR below this threshold will be dropped

def.EOGcorr.enable = false;
def.EOGcorr.corthreshV = 'auto 4';% threshold correlation with vertical EOG
def.EOGcorr.Veogchannames = [];% vertical channel(s)
def.EOGcorr.corthreshH = 'auto 4';% threshold correlation with horizontal EOG
def.EOGcorr.Heogchannames = [];% horizontal channel(s)

def.chancorr.enable = false;
def.chancorr.corthresh = 'auto 4';% threshold correlation
def.chancorr.channames = [];% channel(s)

def.FASTER.enable = false;
def.FASTER.blinkchanname = [];

def.ADJUST.enable = false;

def.MARA.enable = false;

def.opts.FontSize = 14;
def.opts.noplot = 0;
def.opts.nocompute = 0;



function edit_resvar_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_resvar_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_resvar_thresh as text
%        str2double(get(hObject,'String')) returns contents of edit_resvar_thresh as a double
setThreshString(hObject)


% --- Executes during object creation, after setting all properties.
function edit_resvar_thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_resvar_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in check_ADJUST_enable.
function check_ADJUST_enable_Callback(hObject, eventdata, handles)
% hObject    handle to check_ADJUST_enable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_ADJUST_enable
check_enable(hObject,handles)


% --- Executes on button press in check_FASTER_enable.
function check_FASTER_enable_Callback(hObject, eventdata, handles)
% hObject    handle to check_FASTER_enable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_FASTER_enable
check_enable(hObject,handles)



function edit_FASTER_blinkchanname_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FASTER_blinkchanname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_FASTER_blinkchanname as text
%        str2double(get(hObject,'String')) returns contents of edit_FASTER_blinkchanname as a double
numchan = get(hObject,'String');
chname = channamecheck(numchan,hObject);
set(hObject,'String',chname);

function chname = channamecheck(inchan,hObject)
handles = guidata(hObject);
% inchan can be numeric or not
if isempty(inchan)
    chname = '';
elseif isnumeric(inchan)
    [dum dum chname] = chnb(inchan,handles.EEG);
elseif not(isempty(str2num(inchan)))
    % if it's a string of numbers
    [dum dum chname] = chnb(str2num(inchan),handles.EEG);
else
    [dum dum chname] = chnb(inchan,handles.EEG);
end
if isempty(chname)
    chname = 'No channel';
end
    
% --- Executes during object creation, after setting all properties.
function edit_FASTER_blinkchanname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FASTER_blinkchanname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function s = setdef(s,d)
% s = setdef(s,d)
% Merges the two structures s and d recursively.
% Adding the default field values from d into s when not present or empty.

if isstruct(s) && not(isempty(s))
    fields = fieldnames(d);
    for i_f = 1:numel(fields)
        if isfield(s,fields{i_f})
            s.(fields{i_f}) = setdef(s.(fields{i_f}),d.(fields{i_f}));
        else
            s.(fields{i_f}) = d.(fields{i_f});
        end
    end
elseif not(isempty(s))
    s = s;
elseif isempty(s);
    s = d;
end


% --- Executes on button press in check_MARA_enable.
function check_MARA_enable_Callback(hObject, eventdata, handles)
% hObject    handle to check_MARA_enable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_MARA_enable
check_enable(hObject,handles)


% --- Executes on button press in checkNoPlot.
function checkNoPlot_Callback(hObject, eventdata, handles)
% hObject    handle to checkNoPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkNoPlot


% --- Executes on button press in checkReplot.
function checkReplot_Callback(hObject, eventdata, handles)
% hObject    handle to checkReplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkReplot
if get(hObject,'Value')
    set(handles.push_ok,'String','Replot')
else
    set(handles.push_ok,'String','Compute')
end
