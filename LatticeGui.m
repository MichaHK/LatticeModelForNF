function varargout = LatticeGui(varargin)
% LATTICEGUI MATLAB code for LatticeGui.fig
%      LATTICEGUI, by itself, creates a new LATTICEGUI or raises the existing
%      singleton*.
%
%      H = LATTICEGUI returns the handle to a new LATTICEGUI or the handle to
%      the existing singleton*.
%
%      LATTICEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LATTICEGUI.M with the given input arguments.
%
%      LATTICEGUI('Property','Value',...) creates a new LATTICEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LatticeGui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LatticeGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LatticeGui

% Last Modified by GUIDE v2.5 08-Mar-2016 15:03:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LatticeGui_OpeningFcn, ...
                   'gui_OutputFcn',  @LatticeGui_OutputFcn, ...
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


% --- Executes just before LatticeGui is made visible.
function LatticeGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LatticeGui (see VARARGIN)

% Choose default command line output for LatticeGui
handles.output = hObject;
set(handles.P_total_checkbox,'Value',1)
set(handles.L_Stoich_edit,'string','1');
set(handles.A_Stoich_edit,'string','0');
set(handles.M_Stoich_edit,'string','0');
set(handles.H_Stoich_edit,'string','0');
set(handles.Salt_edit,'string','150');
set(handles.P_conf_edit,'String','1');
set(handles.P_ion_edit,'String','1');
set(handles.P_mono_edit,'String','1');
set(handles.L_truncate,'string','0');



set(handles.aa_edit,'String','0.34');
set(handles.cyl_radius_edit,'String','5');
set(handles.chi_flory_edit,'String','0.5');
set(handles.kuhn_length_edit,'String','4');

% 1; 
% h=get(handles.ModelPanel)
% h.SelectedObject=h.Children(1)
% 1;
% ModelPanel_SelectionChangeFcn
% ModelPanel_SelectionChangeFcn(hObject, eventdata, handles)
% handles.ButtonSelection=get(eventdata.NewValue,'Tag');



FilesInFolder=dir([pwd '\ExperimentalData\']);
CSV_files_indicator=cellfun(@(x)~isempty(x),strfind({FilesInFolder.name},'.csv'));
CSV_files={FilesInFolder(CSV_files_indicator>0).name};

currentFolder = pwd;
Zhu_FilesInFolder=dir([currentFolder '\ZhuData']);
Zhu_CSV_files_indicator=cellfun(@(x)~isempty(x),strfind({Zhu_FilesInFolder.name},'.csv'));
Zhu_CSV_files={Zhu_FilesInFolder(Zhu_CSV_files_indicator>0).name};
% for i=1:length(CSV_files)
%     Exp_data(i)=regexp(CSV_files{i},'(?<Types>\w+)_(?<Salt>\d+)mM.csv','names');  
% end
% handles.Exp_data
set(handles.listbox1,'String',CSV_files)
set(handles.listbox2,'String',Zhu_CSV_files)

% Update handles structure
guidata(hObject, handles);

main_lattice(handles);
% main_lattice(handles);

% UIWAIT makes LatticeGui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LatticeGui_OutputFcn(hObject, eventdata, handles) 
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
main_lattice(handles);

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
main_lattice(handles);


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



function A_Stoich_edit_Callback(hObject, eventdata, handles)
% hObject    handle to A_Stoich_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A_Stoich_edit as text
%        str2double(get(hObject,'String')) returns contents of A_Stoich_edit as a double
main_lattice(handles);

% --- Executes during object creation, after setting all properties.
function A_Stoich_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A_Stoich_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function M_Stoich_edit_Callback(hObject, eventdata, handles)
% hObject    handle to M_Stoich_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of M_Stoich_edit as text
%        str2double(get(hObject,'String')) returns contents of M_Stoich_edit as a double
main_lattice(handles);

% --- Executes during object creation, after setting all properties.
function M_Stoich_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to M_Stoich_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function H_Stoich_edit_Callback(hObject, eventdata, handles)
% hObject    handle to H_Stoich_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of H_Stoich_edit as text
%        str2double(get(hObject,'String')) returns contents of H_Stoich_edit as a double
main_lattice(handles);

% --- Executes during object creation, after setting all properties.
function H_Stoich_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to H_Stoich_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function L_Stoich_edit_Callback(hObject, eventdata, handles)
% hObject    handle to L_Stoich_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of L_Stoich_edit as text
%        str2double(get(hObject,'String')) returns contents of L_Stoich_edit as a double
main_lattice(handles);

% --- Executes during object creation, after setting all properties.
function L_Stoich_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to L_Stoich_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Salt_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Salt_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Salt_edit as text
%        str2double(get(hObject,'String')) returns contents of Salt_edit as a double
main_lattice(handles);

% --- Executes during object creation, after setting all properties.
function Salt_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Salt_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in P_total_checkbox.
function P_total_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to P_total_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of P_total_checkbox
main_lattice(handles);

% --- Executes on button press in P_conf_checkbox.
function P_conf_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to P_conf_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of P_conf_checkbox
main_lattice(handles);

% --- Executes on button press in P_ion_checkbox.
function P_ion_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to P_ion_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of P_ion_checkbox
main_lattice(handles);

% --- Executes on button press in P_mono_checkbox.
function P_mono_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to P_mono_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of P_mono_checkbox
main_lattice(handles);



function aa_edit_Callback(hObject, eventdata, handles)
% hObject    handle to aa_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of aa_edit as text
%        str2double(get(hObject,'String')) returns contents of aa_edit as a double
main_lattice(handles);


% --- Executes during object creation, after setting all properties.
function aa_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to aa_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function P_conf_edit_Callback(hObject, eventdata, handles)
% hObject    handle to P_conf_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of P_conf_edit as text
%        str2double(get(hObject,'String')) returns contents of P_conf_edit as a double
main_lattice(handles);


% --- Executes during object creation, after setting all properties.
function P_conf_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to P_conf_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function P_ion_edit_Callback(hObject, eventdata, handles)
% hObject    handle to P_ion_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of P_ion_edit as text
%        str2double(get(hObject,'String')) returns contents of P_ion_edit as a double
main_lattice(handles);


% --- Executes during object creation, after setting all properties.
function P_ion_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to P_ion_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function P_mono_edit_Callback(hObject, eventdata, handles)
% hObject    handle to P_mono_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of P_mono_edit as text
%        str2double(get(hObject,'String')) returns contents of P_mono_edit as a double
main_lattice(handles);


% --- Executes during object creation, after setting all properties.
function P_mono_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to P_mono_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cyl_radius_edit_Callback(hObject, eventdata, handles)
% hObject    handle to cyl_radius_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cyl_radius_edit as text
%        str2double(get(hObject,'String')) returns contents of cyl_radius_edit as a double
main_lattice(handles);


% --- Executes during object creation, after setting all properties.
function cyl_radius_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cyl_radius_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in H_phos.
function H_phos_Callback(hObject, eventdata, handles)
% hObject    handle to H_phos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of H_phos
main_lattice(handles);


% --- Executes on button press in M_phos.
function M_phos_Callback(hObject, eventdata, handles)
% hObject    handle to M_phos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of M_phos
main_lattice(handles);


% --- Executes on button press in L_phos.
function L_phos_Callback(hObject, eventdata, handles)
% hObject    handle to L_phos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of L_phos
main_lattice(handles);


% --- Executes on button press in A_phos.
function A_phos_Callback(hObject, eventdata, handles)
% hObject    handle to A_phos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of A_phos
main_lattice(handles);



function A_value_Callback(hObject, eventdata, handles)
% hObject    handle to A_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A_value as text
%        str2double(get(hObject,'String')) returns contents of A_value as a double
main_lattice(handles);


% --- Executes during object creation, after setting all properties.
function A_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function B_Value_Callback(hObject, eventdata, handles)
% hObject    handle to B_Value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of B_Value as text
%        str2double(get(hObject,'String')) returns contents of B_Value as a double
main_lattice(handles);


% --- Executes during object creation, after setting all properties.
function B_Value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B_Value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function C_Value_Callback(hObject, eventdata, handles)
% hObject    handle to C_Value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of C_Value as text
%        str2double(get(hObject,'String')) returns contents of C_Value as a double
main_lattice(handles);


% --- Executes during object creation, after setting all properties.
function C_Value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to C_Value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function ModelPanel_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to ModelPanel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected object is changed in ModelPanel.
function ModelPanel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in ModelPanel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
% switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
%     case 'radiobutton1'
%         display('Radio button 1');
%     case 'radiobutton2'
%         display('Radio button 2');
%     case 'togglebutton1'
%         display('Toggle button 1');
%     case 'togglebutton2'
%         display('Toggle button 2');
% end
% 1;
handles.ButtonSelection=get(eventdata.NewValue,'Tag');
% if get(eventdata.NewValue,'Tag')==get(eventdata.OldValue,'Tag')
%     return;
% else 
    main_lattice(handles);
% end


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2
main_lattice(handles);


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PlotZhu_checkbox.
function PlotZhu_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to PlotZhu_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PlotZhu_checkbox
main_lattice(handles);



function chi_flory_edit_Callback(hObject, eventdata, handles)
% hObject    handle to chi_flory_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of chi_flory_edit as text
%        str2double(get(hObject,'String')) returns contents of chi_flory_edit as a double
main_lattice(handles);


% --- Executes during object creation, after setting all properties.
function chi_flory_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chi_flory_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kuhn_length_edit_Callback(hObject, eventdata, handles)
% hObject    handle to kuhn_length_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kuhn_length_edit as text
%        str2double(get(hObject,'String')) returns contents of kuhn_length_edit as a double
main_lattice(handles);


% --- Executes during object creation, after setting all properties.
function kuhn_length_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kuhn_length_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ExportModel.
function ExportModel_Callback(hObject, eventdata, handles)
% hObject    handle to ExportModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_lattice(handles);



function L_truncate_Callback(hObject, eventdata, handles)
% hObject    handle to L_truncate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of L_truncate as text
%        str2double(get(hObject,'String')) returns contents of L_truncate as a double
main_lattice(handles);


% --- Executes during object creation, after setting all properties.
function L_truncate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to L_truncate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in p_minimize_checkbox.
function push_minimize_Callback(hObject, eventdata, handles)
% hObject    handle to p_minimize_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_lattice(handles);


% --- Executes on selection change in popupmenu_X.
function popupmenu_X_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_X contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_X
main_lattice(handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_X_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Auto_checkbox.
function Auto_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Auto_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Auto_checkbox
main_lattice(handles);


% --- Executes on button press in P_minimize_checkbox.
function P_minimize_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to P_minimize_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of P_minimize_checkbox
main_lattice(handles);


% --- Executes on selection change in Model_listbox.
function Model_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to Model_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Model_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Model_listbox
main_lattice(handles);


% --- Executes during object creation, after setting all properties.
function Model_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Model_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
main_lattice(handles);



function Ex_exp_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Ex_exp_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ex_exp_edit as text
%        str2double(get(hObject,'String')) returns contents of Ex_exp_edit as a double
main_lattice(handles);


% --- Executes during object creation, after setting all properties.
function Ex_exp_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ex_exp_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
main_lattice(handles);



function A_n_edit_Callback(hObject, eventdata, handles)
% hObject    handle to A_n_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A_n_edit as text
%        str2double(get(hObject,'String')) returns contents of A_n_edit as a double
main_lattice(handles);


% --- Executes during object creation, after setting all properties.
function A_n_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A_n_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in mini_A_n_checkbox.
function mini_A_n_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to mini_A_n_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mini_A_n_checkbox
main_lattice(handles);


% --- Executes on button press in ExportATable.
function ExportATable_Callback(hObject, eventdata, handles)
% hObject    handle to ExportATable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ExportCoeffTable=struct;
FileList=get(handles.listbox1,'String');
handles.CylRadius=[];
handles.EffAlpha=[];
handles.D=[];
handles.V_diff=[];
handles.sum_P_conf=[];
handles.ExportModelName=[];
% guidata(handles.output, handles);

for i=1:numel(FileList)
    filename=FileList{i};
    FileTokens=regexp(filename,'(?<CompositionFromFile>[\w-]*)_(?<SaltFromFile>\w*)mM(?<suffix>\w*)','names');
    ExportCoeffTable(i).TableComposition=FileTokens.CompositionFromFile;
    ExportCoeffTable(i).TableSalt=FileTokens.SaltFromFile;
    set(handles.listbox1,'Value',i);
    set(handles.ExportATable,'Value',0);
    1;
    handles=main_lattice(handles);
    D=handles.D;
    % better set the filename, send to the program main (if possible)
    % and then read charge fraction and A parameters.
    1;
    ExportCoeffTable(i).EffAlpha=handles.EffAlpha;
    A_P_conf_prefactor=get(handles.P_conf_edit,'String');
    ExportCoeffTable(i).A_Norm=A_P_conf_prefactor;
    R=handles.CylRadius;
    D_dependent_term=(D-2*R)./handles.V_diff(D);
    ExportCoeffTable(i).A_simple=max(-handles.sum_P_conf./D_dependent_term);
%     ExportCoeffTable(i).A_simple2=min(handles.sum_P_conf./D_dependent_term);
    1;
end
1;
ExportFileName=[handles.ExportModelName '_TableOfConfPrefactorA'];
% TempExportCell=struct2cell(ExportCoeffTable);
% numel(ExportCoeffTable)
% for j=1:numel(ExportCoeffTable)
%     ExportCoeffTable(1)
% end
TempExportTable=struct2table(ExportCoeffTable);
writetable(TempExportTable,ExportFileName);
disp(['File saved: ' ExportFileName])
% xlswrite('X.xls',TempExportCell)
1;
%     save(ExportFileName,ExportCoeffTable);


% --- Executes on button press in push_MegaPlot.
function push_MegaPlot_Callback(hObject, eventdata, handles)
% hObject    handle to push_MegaPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

FileList=get(handles.listbox1,'String');
handles.x_exp=[];
handles.y_exp=[];
handles.axes2=[];
handles.x=[]; 
handles.y=[];

HorSubPlot=4;
VertSubPlot=4;
j=1;
h=figure;

ind=0;
for i=1:numel(FileList)
    filename=FileList{i};
    FileTokens=regexp(filename,'(?<CompositionFromFile>[\w-]*)_(?<SaltFromFile>\w*)mM(?<suffix>\w*)','names');   
    PlotIt=1;
    set(handles.listbox1,'Value',i);
    set(handles.push_MegaPlot,'Value',0);
    temp=handles.edit25.String{1};
    if handles.pop_megaplot.Value==2
        SaltFilterTokens=regexp(temp,'(?<LowSalt>\w*),(?<HighSalt>\w*)','names');
        if (str2double(FileTokens.SaltFromFile)>str2double(SaltFilterTokens.HighSalt) || str2double(FileTokens.SaltFromFile)<str2double(SaltFilterTokens.LowSalt))
            PlotIt=0;
        end
    end
    if handles.pop_megaplot.Value==3
        if ~strcmp(handles.edit25.String,FileTokens.CompositionFromFile)
            PlotIt=0;
        end
    end
    if PlotIt==1
        ind=ind+1;
        axes(handles.axes1)
        handles=main_lattice(handles);
%         D=handles.D;
        x_exp=handles.x_exp;
        y_exp=handles.y_exp;
        figure(h)
        subplot(HorSubPlot,VertSubPlot,j)
        j=j+1;  
        plot(gca,handles.x,handles.y,'r','LineWidth',2); % the model
        hold on;
        scatter(x_exp,y_exp);                           % the experimental data
        
        xlim([10,100])
        ylim([10,10^7])
        set(gca,'Yscale','log'); % Otherwise the numbers on the axis come out funny. 
        set(gca,'XScale','linear');
        xlabel('\it{D} \rm{(nm)}','FontName','Times New Roman')
        ylabel('\Pi (Pa)','FontName','Times New Roman')
        A_n=get(handles.A_n_edit,'String');
        exponent_A=get(handles.Ex_exp_edit,'String');
        A_P_conf_prefactor=get(handles.P_conf_edit,'String');
        Salt_mM=get(handles.Salt_edit,'String');
        CompOut=RenameCompositionForExport(FileTokens.CompositionFromFile);
        %     legend(['NF:' CompositionFromFile ' , A_{' exponent_A '} = ' A_n],'Fontsize',20)
%         l=legend(['NF:' FileTokens.CompositionFromFile ', k=' A_P_conf_prefactor ', ' 10 FileTokens.SaltFromFile 'mM']);
%         title([ CompOut ', k=' A_P_conf_prefactor ', ' FileTokens.SaltFromFile 'mM']);
        set(gca,'FontSize',12,'FontName','Times New Roman');
        set(gca,'Ytick',[10 10^3 10^5 10^7])
        hold on;
        fakeplt=plot(-10:1:0,-10:1:0,'w.');
        suffix= [];
        if ~isempty(FileTokens.suffix)
            if strcmp('_Eti',FileTokens.suffix)
                suffix=' (*)'
            end
            if strcmp('_Roy',FileTokens.suffix)
                suffix=' (+)'
            end
            if strcmp('_Micha',FileTokens.suffix)
                suffix=' (-)'
            end
        end
        t=legend(fakeplt,{[ CompOut ', \it{k}\rm=' A_P_conf_prefactor ', ' 10 FileTokens.SaltFromFile 'mM' suffix]});
        t.Box='off';
        temp=t.Position;
        t.Location='none';
        t.Position=temp+[0.09 0.05 0.008 0.01];
%         text(0.02,0.98,char('a' + ind - 1),'Units', 'Normalized', 'VerticalAlignment', 'Top')
        text(-0.2,1.3,['(' char('a'+ind-1) ')'],'Units', 'Normalized', 'VerticalAlignment', 'Top')
        axes(handles.axes1)
        
    end
end
figure(h)
fig=gcf;
set(fig,'PaperUnits' , 'inches')
set(fig,'PaperPosition' , [0 0 8.27 10])
set(fig,'PaperPositionMode' , 'manual')
% print('try.eps','-dpng','-r0')
axes(handles.axes1)

save2eps('try2.eps',h)
% --- Executes on selection change in pop_megaplot.
function CompOut=RenameCompositionForExport(CompIn)
CompOut=CompIn;
switch CompIn
    case 'L'
        CompOut='NF-L';
    case 'L5'
        CompOut='NF-L5';
    case 'L11'
        CompOut='NF-L11';
    case 'LM'
        CompOut='NF-L:M';
    case 'L5M'
        CompOut='NF-L5:M';
    case 'L11M'
        CompOut='NF-L11:M';
    case 'LH'
        CompOut='NF-L:H';
    case 'L-deH'    
        CompOut='deNF-L:H';
    case 'L-deM'
        CompOut='deNF-L:M';
    case 'LMH'
        CompOut='NF-L:M:H';
    case 'de-LMH'
        CompOut='deNF-L:M:H';
end


function pop_megaplot_Callback(hObject, eventdata, handles)
% hObject    handle to pop_megaplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_megaplot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_megaplot


% --- Executes during object creation, after setting all properties.
function pop_megaplot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_megaplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double


% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in New_P_Conf_check.
function New_P_Conf_check_Callback(hObject, eventdata, handles)
% hObject    handle to New_P_Conf_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of New_P_Conf_check
main_lattice(handles);


% --- Executes on selection change in Geometry_popupmenu.
function Geometry_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to Geometry_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Geometry_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Geometry_popupmenu
main_lattice(handles);


% --- Executes during object creation, after setting all properties.
function Geometry_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Geometry_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
