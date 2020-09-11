function signal_export(file_list)
% Function to convert .mat files into Excel format
%%

% select files
if ischar(file_list)
    file_list = {file_list};
end
path_idx = strfind(file_list{1},filesep);
path_ini = file_list{1}(1:path_idx(end));
    
% check operating system
if isunix
    % path to toolbox for exporting in linux
    % https://www.mathworks.com/matlabcentral/fileexchange/38591-xlwrite--generate-xls-x--files-without-excel-on-mac-linux-win
    export_sheet_class = '/home/chepe/Science/Tools/Matlab_utilities/export_xls/';
    
    % => linux: add function xlwrite, that works with open office
    javaaddpath([export_sheet_class,'poi_library/poi-3.8-20120326.jar']);
    javaaddpath([export_sheet_class,'poi_library/poi-ooxml-3.8-20120326.jar']);
    javaaddpath([export_sheet_class,'poi_library/poi-ooxml-schemas-3.8-20120326.jar']);
    javaaddpath([export_sheet_class,'poi_library/xmlbeans-2.3.0.jar']);
    javaaddpath([export_sheet_class,'poi_library/dom4j-1.6.1.jar']);
    javaaddpath([export_sheet_class,'poi_library/stax-api-1.0.1.jar']);
    addpath(export_sheet_class)
    expfun = @xlwrite;
    file_sufix = '.xls';
    
elseif ispc || ismac
    % => windows or mac: use microsoft excel for export
    expfun = @xlswrite;
    file_sufix = '.xlsx';
    
else
    error('OS type not found.')
end

if length(file_list)>1
    do_master = true;
    
    [FileName,PathName] = uiputfile([path_ini,'master_file',file_sufix],'Save master file');
    if ~ischar(PathName)
        error('Error when determining master file name!');
    end
    % Sources - Readout
    master_cell = cell(length(file_list),2);
    master_file = [PathName,FileName];
else
    do_master = false;
end

%% Select files and loop over them

% loop
for file_ii = 1:length(file_list)
    source_file = file_list{file_ii};
    disp([' == Reading ',source_file])
    
    vars = whos('-file',source_file);
    if ~ismember('data', {vars.name}) || ~ismember('parameters', {vars.name})
        disp('WARNING: File does not contain data and parameters!')
        continue
    end
    
    if do_master
        master_cell{file_ii,1} = source_file;
    end
    
    %%  Make data structure
    
    data_info = [];
    
    % prepare data handle and read
    data_info.handle = matfile(source_file,'Writable',false);
    
    % get frame rate and num stacks
    data = data_info.handle.data(1,1);
    [~, data_info.ROIs] = size(data_info.handle,'data');
    data_info.frame_rate = data.frame_rate;
    data_info.num_stacks = length(data.signal);
    
    % get parameters
    tmp = load(source_file,'parameters');
    data_info.parameters = tmp.parameters;
    
    % prepare name of the file to save (same name as input, plus .xls) and
    % delete previous
    idx = strfind(source_file,'.mat');
    destination_file = [source_file(1:idx-1),file_sufix];
    if exist(destination_file,'file')
        delete(destination_file);
    end
    
    % save
    data_cell = [{'Source',data_info.handle.source_file;'Time stamp',datestr(data_info.handle.time_stamp);...
        'Frame rate [Hz]',data_info.frame_rate;'# channels',size(data.signal,3);'# frames',data_info.num_stacks;'# ROI',data_info.ROIs;[],[];'Parameters',[]};...
        fieldnames(data_info.parameters),struct2cell(data_info.parameters)];
    if ~do_master
        expfun(destination_file,data_cell,'Data');
    end
    
    %% Extract from each ROI
    
    data_cell = {'Label','ROI Type','Mean dF/F %','Std dF/F %','# Events','Events per minute','Mean event signal'};
    for roi_num = 1:data_info.ROIs
        data = data_info.handle.data(1,roi_num);
        
        data_cell{data.label+1,1} = data.label;
        data_cell{data.label+1,2} = data.type;
        if strcmp(data.type,'Neuron') || strcmp(data.type,'Glial cell')
            data_cell{data.label+1,3} = 100*mean(data.processed(1,:,1));
            data_cell{data.label+1,4} = std(100*data.processed(1,:,1));
            data_cell{data.label+1,5} = size(data.events,1);
            data_cell{data.label+1,6} = size(data.events,1)/(data_info.num_stacks/data_info.frame_rate/60);
            data_cell{data.label+1,7} = max([0 mean(data.processed(1,data.events(:,1),1))]);
        else
            data_cell{data.label+1,3} = mean(data.signal(1,:,1));
            data_cell{data.label+1,4} = std(data.signal(1,:,1));
        end
    end
    
    if ~do_master
        expfun(destination_file,data_cell,'Summary');
    else
        master_cell{file_ii,2} = data_cell;
    end
    
    %% Export individual results
    
    if ~do_master
        vars = whos('-file',source_file);
        if ismember('analysis_info', {vars.name})
            load(source_file,'analysis_info');
            load(source_file,'id_type_label_pos');
            % save processed results
            if ~isempty(analysis_info.correlation)
                data_cell = [{'ROI'},id_type_label_pos(:,3)';id_type_label_pos(:,3),mat2cell(analysis_info.correlation,ones(size(analysis_info.correlation,1),1),ones(1,size(analysis_info.correlation,1)))];
                expfun(destination_file,data_cell,'Correlation');
            end
            if ~isempty(analysis_info.ch_ratio)
                data_cell = [{'ROI Label','Initial','Final'};...
                    mat2cell(analysis_info.ch_ratio,ones(size(analysis_info.ch_ratio,1),1),ones(1,3))];
                expfun(destination_file,data_cell,'Channel ratio');
            end
            if ~isempty(analysis_info.part_rate)
                data_cell = [{'ROI Label','Cell centered %','Population centered %'};...
                    mat2cell(analysis_info.part_rate,ones(size(analysis_info.part_rate,1),1),ones(1,3))];
                expfun(destination_file,data_cell,'Participation rate');
            end
            
        end
        
        %% ready!
        disp(' == Finished.')
    end
    
end

%% Make master file
if do_master
    
    % delete previous master
    if exist(master_file,'file')
        delete(master_file);
    end
    
    % save list of source files
    data_cell = master_cell(:,1);
    expfun(master_file,data_cell,'Source files');
    
    % save list of ROIs
    roi_num = 0;
    for ii = 1:size(master_cell,1)
        roi_num = max([roi_num,size(master_cell{ii,2},1)])-1;
    end
    data_cell = cell(roi_num,1);
    for ii = 1:size(master_cell,1)
        for jj = 2:size(master_cell{ii,2},1)
            label = master_cell{ii,2}{jj,1};
            if ~isempty(label)
                data_cell{label,1} = master_cell{ii,2}{jj,2};
            end
        end
    end
    expfun(master_file,data_cell,'ROI type');
    
    % save mean dF/F
    data_cell = cell(roi_num,size(master_cell,1));
    for ii = 1:size(master_cell,1)
        for jj = 2:size(master_cell{ii,2},1)
            label = master_cell{ii,2}{jj,1};
            if ~isempty(label)
                data_cell{label,ii} = master_cell{ii,2}{jj,3};
            end
        end
    end
    expfun(master_file,data_cell,'Mean dFF');
    
    % save std dF/F
    data_cell = cell(roi_num,size(master_cell,1));
    for ii = 1:size(master_cell,1)
        for jj = 2:size(master_cell{ii,2},1)
            label = master_cell{ii,2}{jj,1};
            if ~isempty(label)
                data_cell{label,ii} = master_cell{ii,2}{jj,4};
            end
        end
    end
    expfun(master_file,data_cell,'Std dFF');
    
    % save event rate
    data_cell = cell(roi_num,size(master_cell,1));
    for ii = 1:size(master_cell,1)
        for jj = 2:size(master_cell{ii,2},1)
            label = master_cell{ii,2}{jj,1};
            if ~isempty(label)
                data_cell{label,ii} = master_cell{ii,2}{jj,6};
            end
        end
    end
    expfun(master_file,data_cell,'Events per minute');
    
    % save signal event
    data_cell = cell(roi_num,size(master_cell,1));
    for ii = 1:size(master_cell,1)
        for jj = 2:size(master_cell{ii,2},1)
            label = master_cell{ii,2}{jj,1};
            if ~isempty(label)
                data_cell{label,ii} = master_cell{ii,2}{jj,7};
            end
        end
    end
    expfun(master_file,data_cell,'Mean event signal');
    
end
end