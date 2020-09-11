function gui_start
% Opens a GUI where the source folder can be selected and the different
% steps of the analysis can be accessed.
%
% The file "gui_path.txt" has the list of folders that can be quickly
% accessed. New lines can be added directly to the text file or via the
% GUI.
%%

gui_info = [];

gui_info.S1_list = {};
gui_info.S3_list = {};
gui_info.S5_list = {};

gui_info.tiff_structure = cell(0,1);

% get source path
sourcepath = mfilename('fullpath');
idx = strfind(sourcepath,filesep);
gui_info.source_path = sourcepath(1:idx(end));

% read file "gui_path.txt" to get initial path list
gui_info.data_path = read_path(gui_info.source_path);

% create and modify the GUI
[~,h] = gui_create();


%%
% ----- GUI generation and update
    function [hfig,h] = gui_create()
        
        % make new figure
        hfig = open([gui_info.source_path,'gui_start.fig']);
        set(hfig,'Tag','inifig')
        h = guihandles(hfig);
        set(hfig,'Name','Start window')
        
        % ------ PATH
        set(h.dataPath,'String',gui_info.data_path)
        set(h.addNew,'callback',@addNew)
        
        % ------ STEP 1: Pre process
        set(h.S1_align,'Value',1)
        set(h.S1_range,'Value',0)
        set(h.S1_maxShift,'String','5')
        set(h.S1_concatenate,'Value',true)
        
        set(h.S1_nncorr,'Value',0)
        set(h.S1_chunkyMax,'Value',0)
        
        set(h.S1_list,'String',gui_info.S1_list)
        set(h.S1_add,'callback',{@list_manager,'S1_list','add'})
        set(h.S1_delete,'callback',{@list_manager,'S1_list','delete'})
        set(h.S1_run,'callback',@call_preprocess)
        
        % ------- STEP 2: Segmentation
        set(h.S2_loadExport,'Value',0)
        set(h.S2_run,'callback',@call_segmentation)
        
        % ------- STEP 3: Signal extraction
        set(h.S3_list,'String',gui_info.S3_list)
        set(h.S3_add,'callback',{@list_manager,'S3_list','add'})
        set(h.S3_delete,'callback',{@list_manager,'S3_list','delete'})
        set(h.S3_run,'callback',@call_extraction)
        
        % ------- STEP 4: Analysis
        set(h.S4_run,'callback',@call_process)
        
        % ------- STEP 5: Join result
        set(h.S5_list,'String',gui_info.S5_list)
        set(h.S5_add,'callback',{@list_manager,'S5_list','add'})
        set(h.S5_delete,'callback',{@list_manager,'S5_list','delete'})
        set(h.S5_run,'callback',@call_export)
    end

    function gui_update(updateTag,hObject,value)
        switch updateTag
            case 'update_color'
                if value
                    set(hObject,'BackgroundColor',[0.9333 0.8667 0.5098])
                    set(hObject,'Enable','off')
                else
                    set(hObject,'BackgroundColor',[1 1 1]*0.941)
                    set(hObject,'Enable','on')
                end
        end
        drawnow
    end

% ----- function to add a new path to the list
    function addNew(hObject,eventdata)
        new_folder_name = uigetdir(gui_info.data_path{end});
        gui_info.data_path{length(gui_info.data_path)+1} = [new_folder_name,filesep];
        
        set(h.dataPath,'String',gui_info.data_path)
        set(h.dataPath,'Value',length(gui_info.data_path))
        
        save_path(gui_info.source_path,new_folder_name)
    end

% ----- function caller
    function call_preprocess(hObject,eventdata)
        
        gui_update('update_color',hObject,true)
        
        try
            do_align = get(h.S1_align,'Value');
            do_clip = get(h.S1_range,'Value');
            max_shift = str2double(get(h.S1_maxShift,'String'));
            do_corr = get(h.S1_nncorr,'Value');
            do_chunkyMax = get(h.S1_chunkyMax,'Value');
            
            for set_num = 1:length(gui_info.tiff_structure)
                disp(['===== Opening data set ',gui_info.S1_list{set_num}])
                helper_preprocess(gui_info.tiff_structure{set_num},do_align,do_clip,max_shift,do_corr,do_chunkyMax)
            end
            disp('===== FINISHED!')
            gui_update('update_color',hObject,false)
        catch ierr
            gui_update('update_color',hObject,false)
            rethrow(ierr)
        end
    end

    function call_segmentation(hObject,eventdata)
        gui_update('update_color',hObject,true)
        try
            ini_path = gui_info.data_path{get(h.dataPath,'Value')};
            if get(h.S2_loadExport,'Value')
                gui_segmentation('exported',ini_path)
            elseif get(h.S2_loadSelected,'Value')
                if ~isempty(gui_info.tiff_structure)
                    gui_segmentation('file_structure',gui_info.tiff_structure{get(h.S1_list,'Value')})
                else
                    gui_segmentation('gui',ini_path)
                end
            else
                gui_segmentation('gui',ini_path)
            end
            gui_update('update_color',hObject,false)
        catch ierr
            gui_update('update_color',hObject,false)
            rethrow(ierr)
        end
    end

    function call_extraction(hObject,eventdata)
        gui_update('update_color',hObject,true)
        
        try
            for file_ii = 1:length(gui_info.S3_list)
                disp('===========================================')
                disp(['Reading data from: ',gui_info.S3_list{file_ii}])
                
                destination_file = signal_extract(gui_info.S3_list{file_ii});
                
                disp('Finished reading data')
                disp(['Extracted signal saved in: ',destination_file])
                disp('===========================================')
            end
            gui_update('update_color',hObject,false)
        catch ierr
            gui_update('update_color',hObject,false)
            rethrow(ierr)
        end
    end

    function call_process(hObject,eventdata)
        gui_update('update_color',hObject,true)
        
        try
            ini_path = gui_info.data_path{get(h.dataPath,'Value')};
            gui_analyze(ini_path)
            gui_update('update_color',hObject,false)
        catch ierr
            gui_update('update_color',hObject,false)
            rethrow(ierr)
        end
    end

    function call_export(hObject,eventdata)
        gui_update('update_color',hObject,true)
        try
            if length(gui_info.S5_list)>0
                signal_export(gui_info.S5_list)
            end
            gui_update('update_color',hObject,false)
        catch ierr
            gui_update('update_color',hObject,false)
            rethrow(ierr)
        end
    end

% ------ manage list of files
    function list_manager(hObject,eventdata,current_list,task)
        switch task
            case 'delete'
                if get(h.(current_list),'Value')<=length(gui_info.(current_list))
                    gui_info.(current_list)(get(h.(current_list),'Value')) = [];
                    % if pre processing, remove from tif structure
                    if strcmp(current_list,'S1_list')
                        gui_info.tiff_structure(get(h.(current_list),'Value')) = [];
                    end
                    set(h.(current_list),'Value',1)
                    set(h.(current_list),'String',gui_info.(current_list))
                end
            case 'add'
                
                ini_path = gui_info.data_path{get(h.dataPath,'Value')};
                
                switch current_list
                    
                    case {'S3_list','S5_list'}
                        % get ini path
                        if strcmp(current_list,'S3_list')
                            ini_path_text = [gui_info.data_path{get(h.dataPath,'Value')},'*.mat'];
                        elseif strcmp(current_list,'S5_list')
                            ini_path_text = [gui_info.data_path{get(h.dataPath,'Value')},'*_extracted.mat'];
                        end
                        % select the file
                        [FileName,PathName] = uigetfile(ini_path_text,'MultiSelect','on','Select ROI file to open');
                        % update path or leave
                        if ~ischar(PathName)
                            error('Error selecting files to open')
                        end
                        % if multiple selected, add one by one
                        file_ii = length(gui_info.(current_list));
                        if iscell(FileName)
                            for f_ii = 1:length(FileName)
                                file_ii = file_ii +1;
                                gui_info.(current_list){file_ii} = [PathName,FileName{f_ii}];
                            end
                        else
                            % add to structure
                            file_ii = file_ii +1;
                            gui_info.(current_list){file_ii} = [PathName,FileName];
                        end
                        % remove duplicates
                        gui_info.(current_list) = unique(gui_info.(current_list)');
                        % change gui
                        set(h.(current_list),'String',gui_info.(current_list))
                        
                    case 'S1_list'
                        
                        % make list of tiffs for each channel
                        num_channels = 0;
                        file_names = [];
                        while num_channels <= 2
                            
                            num_channels = num_channels + 1;
                            
                            % get file names for next channel
                            [FileName,PathName] = uigetfile({[ini_path,'/*.tif']},'MultiSelect','on',['Select tif files for channel ',num2str(num_channels)]);
                            if ischar(PathName)
                                ini_path = PathName;
                            end
                            
                            % stop adding channels if canceled
                            if isnumeric(FileName)
                                num_channels = num_channels - 1;
                                break
                            end
                            
                            % convert input into cell "file_names"
                            if iscell(FileName) %i.e multiple files
                                file_names_channel = cell(length(FileName),1);
                                for ii=1:length(FileName)
                                    file_names_channel{ii} = [PathName,FileName{ii}];
                                end
                            else %i.e only one file
                                file_names_channel = cell(1,1);
                                file_names_channel{1,1} = [PathName,FileName];
                            end
                            
                            % if it is the first channel, just put in list
                            if num_channels == 1
                                file_names = file_names_channel;
                            else
                                if size(file_names,1)~= size(file_names_channel,1)
                                    error('Different number of files defined as for previous channels')
                                end
                                file_names = [file_names,file_names_channel];
                            end
                        end
                        
                        if num_channels == 0
                            return
                        end
                        
                        % check if format is to be concatenated
                        if get(h.S1_concatenate,'Value')
                            % join
                            gui_info.tiff_structure{length(gui_info.tiff_structure)+1} = file_names;
                            % add name
                            name_idx = strfind(file_names{1},filesep);
                            gui_info.S1_list{length(gui_info.S1_list)+1} = file_names{1}(name_idx(end-1):end);
                            % add to GUI
                            set(h.S1_list,'String',gui_info.S1_list)
                        else
                            % concatenate
                            for fnum = 1:size(file_names,1)
                                gui_info.tiff_structure{length(gui_info.tiff_structure)+1} = file_names(fnum,:);
                                % add name
                                name_idx = strfind(file_names{fnum,1},filesep);
                                gui_info.S1_list{length(gui_info.S1_list)+1} = file_names{fnum,1}(name_idx(end-1):end);
                                % add to GUI
                                set(h.S1_list,'String',gui_info.S1_list)
                            end
                        end
                        
                        
                        
                end
                
        end
        
    end
end

%% FUNCTIONS TO READ AND WRITE SOURCE FILE 'gui_path.txt'

function data_path = read_path(sourcepath)
% open file
fh = fopen([sourcepath,'data_path.txt'],'r');
data_path = cell(0,1);
tline = fgetl(fh);
while ischar(tline)
    % jump comments and empty spaces
    if strcmp(tline,'')
        tline = fgetl(fh);
        continue
    end
    if strcmp(tline(1),'%')
        tline = fgetl(fh);
        continue
    end
    % add line to cell and read next
    data_path = cat(1,data_path,{tline});
    tline = fgetl(fh);
end
% close
fclose(fh);
% add current folder at the end
data_path{length(data_path)+1} = [pwd,filesep];
end

function save_path(sourcepath,line2add)
% open file, print new line and close
fh = fopen([sourcepath,'data_path.txt'],'a');
fprintf(fh,['\n',line2add,filesep]);
fclose(fh);
end

%% FUNCTION TO MAKE A LIST OF FILES TO PRE-PROCESS AND DO SO

function helper_preprocess(data_set,do_align,do_clip,max_shift,do_corr,do_chunkyMax)
% With this function the user can selet multiple tif data-sets to be
% pre-processed offline => Mean images and alignment are calculated and
% stored in the respective meta_data.mat. This is used to save time when
% opening the files to read create the ROIs with do_segmentation.m

%% READ FILES

% create tif object
clear obj
obj = tiff_reader(data_set);

% clip values?
if do_clip
    obj.do_clip(true)
    disp('Data clipped')
end

if do_align
    disp(['-- Calculating properties and alignment using channel ',num2str(size(data_set,2)),'.'])
    try
    obj.alignment_calculate(size(data_set,2),false,max_shift)
    % check if alignment is not shifting too much between frames
    catch
            obj.metadata_delete
            disp('-- Alignment failed: Too large x-y shifts found.')
            disp('   Calculating properties without alignment..')
            obj.properties_calculate(true,false)
    end
else
    disp('-- Calculating properties.')
    obj.properties_calculate(true,false)
end
if do_corr
    disp('-- Calculating correlation.')
   obj.pix_product_calculate(true,false); 
end
if do_chunkyMax
    disp('-- Calculating averaged maximum.')
    obj.chunkyMax_calculate(false);
end

disp('-- Saving results.')
obj.metadata_save


end
