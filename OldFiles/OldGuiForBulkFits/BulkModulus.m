function varargout = BulkModulus(varargin)
% BULKMODULUS MATLAB code for BulkModulus.fig
%      BULKMODULUS, by itself, creates a new BULKMODULUS or raises the existing
%      singleton*.
%
%      H = BULKMODULUS returns the handle to a new BULKMODULUS or the handle to
%      the existing singleton*.
%
%      BULKMODULUS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BULKMODULUS.M with the given input arguments.
%
%      BULKMODULUS('Property','Value',...) creates a new BULKMODULUS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BulkModulus_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BulkModulus_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BulkModulus

% Last Modified by GUIDE v2.5 25-Sep-2015 00:01:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BulkModulus_OpeningFcn, ...
                   'gui_OutputFcn',  @BulkModulus_OutputFcn, ...
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


% --- Executes just before BulkModulus is made visible.
function BulkModulus_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BulkModulus (see VARARGIN)

% Choose default command line output for BulkModulus
handles.output = hObject;


FilesInFolder=dir([pwd '\ExperimentalData\']);
CSV_files_indicator=cellfun(@(x)~isempty(x),strfind({FilesInFolder.name},'.csv'));
CSV_files={FilesInFolder(CSV_files_indicator>0).name};

% currentFolder = pwd;
set(handles.listbox1,'String',CSV_files)

% Update handles structure
guidata(hObject, handles);
% 1;
% get(handles.listbox1,'String')
main_bulk(handles)
% main_bulk(handles)

% UIWAIT makes BulkModulus wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BulkModulus_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
main_bulk(handles)

% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Command_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Command_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Command_edit as text
%        str2double(get(hObject,'String')) returns contents of Command_edit as a double
main_bulk(handles)


% --- Executes during object creation, after setting all properties.
function Command_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Command_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
