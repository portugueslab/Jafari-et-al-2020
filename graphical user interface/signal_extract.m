function destination_file = signal_extract(input_file,parameters, suffix)
% use an exported ROI with the indications to load a set of tifs to read
% the data and export it
% A Version 7.3 MAT-file is generated, where the traces can be opened
% individually, saving memory space when doing the GUI analysis.

%% default parameters for processing

if nargin ==1
    parameters.NP_fraction = 0.6;
    parameters.Fo_range = 6;
    parameters.Fo_quantile = 0.25;
    parameters.Filter = 0.2;
    parameters.signal_threshold = 0.15;
    parameters.probability_threshold = 1;    
    suffix = 'extracted';
end

if nargin == 2
    suffix = 'extracted';
end

%%    
    opened = load(input_file);
    % check that everything that is needed is in the file
    if ~isfield(opened,'ROI') || ~isfield(opened,'data_class')
        disp('WARNING: File does not contain the necessary information to load data + ROI!')
        return
    end
    
    % prepare name of the file to save (same name as input, plus _extracted
    idx = strfind(input_file,'.mat');
    destination_file = [input_file(1:idx-1),'_',suffix,'.mat'];
    
    % check if destination file exist and load the parameters from it
    if exist(destination_file,'file')
        load(destination_file,'parameters')
    end
    
    %% load the data => data_info / parameters.frame_rate
    
    data_info = [];
    data_info.handle = tiff_reader(opened.data_class.files);
    
    % get data properties
    data_info.width = data_info.handle.width;
    data_info.height = data_info.handle.height;
    data_info.num_stacks = data_info.handle.num_stacks;
    data_info.num_channels = data_info.handle.num_channels;
    
    % change correction
    if isfield(opened.data_class,'x_correction')
        data_info.handle.set_x_correction(opened.data_class.x_correction)
    end
    % change clip value
    if isfield(opened.data_class,'clip_value')
        if any(opened.data_class.clip_value>0)
            data_info.handle.do_clip(true,opened.data_class.clip_value)
        end
    end
    % change frame rate
    if isfield(opened.data_class,'frame_rate')
        data_info.frame_rate = opened.data_class.frame_rate;
    else
        data_info.frame_rate = 15;
    end
    
    %% load the ROIs => ROI / ROI_count / ROI_label / ROI_type / ROI_filth
    
    % convert 2D to 3D array if necessary
    if ndims(opened.ROI)==2
        ROI_count = max(opened.ROI(:));
        ROI = cell(ROI_count,1);
        for roi_ii = 1:ROI_count
            [y,x] = find(opened.ROI == roi_ii);
            ROI{roi_ii} = [y,x];
        end
    else
        ROI_count = size(opened.ROI,3);
        ROI = cell(ROI_count,1);
        for roi_ii = 1:ROI_count
            [y,x] = find(opened.ROI(:,:,roi_ii));
            ROI{roi_ii} = [y,x];
        end
    end
    
    % load ROI labels or set to ROI number
    if isfield(opened,'ROI_label')
        ROI_label = opened.ROI_label;
    else
        ROI_label = 1:ROI_count;
    end
    
    % load ROI types or set to Neuron
    if isfield(opened,'ROI_type')
        ROI_type = opened.ROI_type;
    else
        ROI_type = repmat({'Neuron'},[ROI_count,1]);
    end
    
    % change filth ROI
    if isfield(opened,'ROI_filth')
        ROI_filth = opened.ROI_filth;
    else
        ROI_filth = false(data_info.height,data_info.width);
    end
    
    %% make background ROI (surrounding neuropile 2.5*radius)
    
    ROIisNP = strcmp(opened.ROI_type,'Neuropile');
    ROIisBV = strcmp(opened.ROI_type,'Blood vessel');
    
    % == make enlarged total ROI for background
    large_ROIs = false(data_info.height,data_info.width);
    for roi_ii = 1:ROI_count
        if ~strcmp(ROI_type{roi_ii},'Neuron') || strcmp(ROI_type{roi_ii},'Glial cell')
            continue
        end
        large_ROIs(sub2ind(size(large_ROIs),ROI{roi_ii}(:,1),ROI{roi_ii}(:,2))) = true;    
    end
    large_ROIs = imdilate(imfill(large_ROIs,'holes'),strel('disk',4));
    
    % calculate center of mass of the ROIs
    center_of_mass = zeros(ROI_count,3);
    for roi_ii=1:ROI_count
        ind_y = ROI{roi_ii}(:,1);
        ind_x = ROI{roi_ii}(:,2);
        center_of_mass(roi_ii,:) = [round(mean(ind_y)),round(mean(ind_x)),...
            mean(sqrt((ind_y-mean(ind_y)).^2+(ind_x-mean(ind_x)).^2))];
    end
    mean_radius = round(mean(center_of_mass(~ROIisNP & ~ROIisBV,3)));
    bg_region = [floor(-3*mean_radius):ceil(3*mean_radius)];
    
    % == get background AS SINGLE INDEX
    ROI_bg = cell(ROI_count,1);
    for roi_ii=1:ROI_count
        if ~ROIisNP(roi_ii)
            % make background ROI as squares 3x the mean radius
            bg_ROI = false(data_info.height,data_info.width);
            bg_ROI(...
                min([data_info.height*ones(size(bg_region));max([ones(size(bg_region));center_of_mass(roi_ii,1)+bg_region])]),...
                min([data_info.width*ones(size(bg_region));max([ones(size(bg_region));center_of_mass(roi_ii,2)+bg_region])]))=true;
            % remove cells and filth
            bg_ROI = bg_ROI.*~logical(large_ROIs + ROI_filth);
            % get indices
            ROI_bg{roi_ii} = find(bg_ROI);
        end
    end
    
    % == remove fith from ROI and CONVERT TO SINGLE INDEX
    for roi_ii = 1:ROI_count
        tmp_ROI = false(data_info.height,data_info.width);
        tmp_ROI(sub2ind(size(tmp_ROI),ROI{roi_ii}(:,1),ROI{roi_ii}(:,2))) = true; 
        tmp_ROI = tmp_ROI .* ~ROI_filth;
        ROI{roi_ii} = find(tmp_ROI);
    end
    
    %% Read the data (neuropile only on first channel!!)
    % Ignore the NaN in the calculation (from e.g. alignment)
    
    signal = zeros(ROI_count,data_info.num_stacks,data_info.num_channels);
    neuropile = zeros(ROI_count,data_info.num_stacks);
    
    % first channel
    for frame = 1:data_info.num_stacks
        img = double(data_info.handle.read_stack(frame,1));
        for roi_ii = 1:ROI_count
            tmp = img(ROI{roi_ii});
            signal(roi_ii,frame,1) = sum(tmp(~isnan(tmp)))/sum(~isnan(tmp));
            tmp = img(ROI_bg{roi_ii});
            neuropile(roi_ii,frame) = sum(tmp(~isnan(tmp)))/sum(~isnan(tmp));
        end
        
        % the rest of the channels
        for channel = 2:data_info.num_channels
            img = double(data_info.handle.read_stack(frame,channel));
            for roi_ii = 1:ROI_count
                tmp = img(ROI{roi_ii});
                signal(roi_ii,frame,channel) = sum(tmp(~isnan(tmp)))/sum(~isnan(tmp));
            end
        end
        
    end
    
    %% Convert to structure with independent info
    
    data = struct('label',[],'type',[],'frame_rate',[],'x_y_pos',[],'ROI_pixels',[],...
        'signal',[],'neuropile',[],'processed',[],'events',[]);
    
    for roi_ii=1:ROI_count
        data(roi_ii).label = ROI_label(roi_ii);
        data(roi_ii).type = ROI_type{roi_ii};
        data(roi_ii).frame_rate = data_info.frame_rate;
        data(roi_ii).x_y_pos = center_of_mass(roi_ii,[2 1]);
        data(roi_ii).ROI_pixels = size(ROI{roi_ii},1);
        
        data(roi_ii).signal = signal(1,:,:);
        signal(1,:,:) = [];
        
        if ~ROIisNP(roi_ii)
            data(roi_ii).neuropile = neuropile(1,:);
        end
        neuropile(1,:) = [];
        
    end
    
    %% process data and estimate events with default parameters
    
    data_range.raw.min = [];
    data_range.raw.max = [];
    data_range.processed.min = [];
    data_range.processed.max = [];
    
    for roi_ii=1:ROI_count
        % do only for Neurons and Glial cells
        if ROIisNP(roi_ii)||ROIisBV(roi_ii)
            continue
        end
        % process
        [signal_processed,signal_events] = signal_process(data(roi_ii),parameters,true,true);
        % add to results
        data(roi_ii).processed = signal_processed;
        data(roi_ii).events = signal_events;
        
        % get ranges
        data_range.raw.min = min([data_range.raw.min min(data(roi_ii).neuropile) min(data(roi_ii).signal(1,:,1))]);
        data_range.raw.max = max([data_range.raw.max max(data(roi_ii).neuropile) max(data(roi_ii).signal(1,:,1))]);
        data_range.processed.min = min([data_range.processed.min min(data(roi_ii).processed)]);
        data_range.processed.max = max([data_range.processed.max max(data(roi_ii).processed)]);
    end
    
    %% save data and parameters
    
    id_type_label_pos = [mat2cell((1:ROI_count)',ones(ROI_count,1)) ROI_type mat2cell(ROI_label,ones(ROI_count,1)) mat2cell(center_of_mass(:,[2 1]),ones(ROI_count,1))];
    time_stamp = now;
    source_file = input_file;
    save(destination_file,'data','parameters','id_type_label_pos','source_file','data_range','time_stamp','-v7.3')
    
end