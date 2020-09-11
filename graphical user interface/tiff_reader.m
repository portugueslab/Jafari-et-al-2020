classdef tiff_reader < handle
    %% class that handles experimental data in tiff format
    % Manage the recording of an imaging session that has multiple tif files,
    % e.g. for a recording split in different sets or multiple channels.
    % The class also handles the x-y alignment of the data using one of the
    % channel as reference (default is last one added) and
    % loads/calculates/loads the recording properties (mean, var, max, min)
    % such that it can be used for the analysis, e.g. ROI definition.
    %
    % HOW TO USE: Example
    % - obj = tiff_reader(path_to files)
    %       select the fileS for the first channel, then second, etc until
    %       the GUI is cancelled
    % - obj.alignment_calculate() 
    %       use last added channel to calculate alignment within the file
    % - obj.alignment_reference() 
    %       use last added channel to calculate alignment with a reference file
    % - obj.metadata_save
    %       save the properties and alignment such that it does not have to be calculated
    %       again (saved in meta_data.mat in the same folder as the data)
    % - img = obj.read_stack(frame_number,channel)
    %       read a specific frame and a specific channel. All corrections
    %       are applied
    % - img = obj.mean
    %       get the mean of the recording. Img has all channels
    %       included, i.e size(img) = [dimx, dimy, channels]
    %
    % METHOD LIST:
    % - tiff_reader (see function for constructor options)
    % - read_stack
    % - mean / variance / maximum / minimum / maxValue / minValue / xyShift
    % - alignment_set / alignment_calculate / alignment_reference
    % - properties_calculate
    % - metadata_delete / metadata_save / metadata_load
    %
    % PARAMETERS:
    % - alignment_channel: which channel is used for alignment (def. last one added)
    %
    % chepe@nld.ds.mpg.de
    %%
    
    properties (SetAccess = private)
        tiff_info
        width
        height
        num_channels
        num_stacks
        x_correction = 0
        clip_value = [];
        
    end
    properties (Hidden)
        stack_to_file
        
        alignment_channel
        img_T = [];
        
        % properties of all the files together
        % as calculated when doing alignment or when the individual sets
        % are joined
        img_mean = [];
        img_variance = [];
        img_maximum = [];
        img_minimum = [];
        img_minValue = [];
        img_maxValue = [];
        img_averagedMaximum = [];
        img_correlation = [];
        
        % properties of each set individually.
        % used for internal purposes, e.g. saving, adding channels without
        % calculating from scratch
        set_mean = {};
        set_variance = {};
        set_averagedMaximum = {};
        set_maximum = {};
        set_minimum = {};
        set_minValue = {};
        set_maxValue = {};
        set_pix_prod = {};
    end
    
    %%
    methods
        
        %% Constructor and destructor
        % NOTES:
        % The class is initiated either with a cell list of files (rows = files, columns = channel)
        % or with a graphical user interface
        
        function obj = tiff_reader(input_variable)
            % input_variable is either:
            % input is a single string = path for the GUI to open
            % input is a cell = structure of .tif file names for each channel
            
            % ----------- Parse input variable, which is either:
            % a single string or nothing = path for the GUI to open
            % a cell = structure of .tif file names for each channel
            if nargin==0
                input_variable = [];
            end
            file_names = obj.make_file_structure(input_variable);

            % ------------ Check data
            % Do all files have the same frame size?
            tmp=imfinfo(file_names{1});
            width = tmp(1).Width;
            height = tmp(1).Height;
            for ii=2:numel(file_names)
                tmp=imfinfo(file_names{ii});
                if tmp(1).Width~=width || tmp(1).Height~=height
                    error('Not all files have the same x-y dimensions')
                end
            end
            
            % Do all channels have the same number of stacks?
            for ii=1:size(file_names,1)
                tmp=imfinfo(file_names{ii,1});
                for jj=2:size(file_names,2)
                    tmp2=imfinfo(file_names{ii,jj});
                    if length(length(tmp))~=length(length(tmp2))
                        error('Not all files have the same number of stacks')
                    end
                end
            end
            
            % ------------  Open files and extract info
            
            % open files to read and create tiff_info structure.
            obj.tiff_info = struct();
            for ii=1:size(file_names,1)
                for jj=1:size(file_names,2)
                    tmp=imfinfo(file_names{ii,jj});
                    obj.tiff_info(ii,jj).tifflib=Tiff(tmp(1).Filename,'r');
                    obj.tiff_info(ii,jj).Filename=tmp(1).Filename;
                    if length(tmp)>1
                        obj.tiff_info(ii,jj).Stacks=length(tmp);
                    else
                        disp('Warning: Only 1 stack found in the file. If file siz >4GB try splitting it into multiple files.')
                        obj.tiff_info(ii,jj).Stacks=1;
                    end
                end
            end
            
            % Get the total number of stacks and a mapping "stack -> corresponding file"
            obj.num_stacks = 0;
            obj.stack_to_file = [];
            for ii=1:size(obj.tiff_info,1)
                % [file_number stack_ini stack_end]
                obj.stack_to_file = [obj.stack_to_file; ii,obj.num_stacks+[1,obj.tiff_info(ii,1).Stacks]];
                obj.num_stacks = obj.num_stacks + obj.tiff_info(ii,1).Stacks;
            end
            
            disp('-- Tiff files opened...')
            
            % ------------  Load properties
            tmp=imfinfo(file_names{1});
            obj.width = tmp(1).Width;
            obj.height = tmp(1).Height;
            obj.num_channels = size(file_names,2);
            
            % try to load file properties (for each file individually)
            metadata_delete(obj)            
            obj.alignment_channel = obj.num_channels;
            obj.metadata_load;
            
        end
        function file_names = make_file_structure(obj,input_variable)
            % Construct list of file names from input
            
            if isempty(input_variable)
                % use GUI to choose files and start from current directory
                file_names = [];
                ini_path = pwd;
            elseif iscell(input_variable)
                % file_names cell is given in input
                file_names = input_variable;
            elseif ischar(input_variable)
                % input is either path for GUI or a single file
                if strcmp(input_variable(end-3:end),'.tif')
                    file_names = {input_variable};
                else
                    file_names = [];
                    ini_path = input_variable;
                end
            else
                error('The input to tiff_reader is either a string of the path for the GUI or a cell of files {rows = file_num, columns = channel}.')
            end
            
            % select the file_names using a GUI, one channel at the time,
            % where for each channel multiple files can be combined
            if isempty(file_names)
                
                while true
                    
                    % get file names for next channel
                    [FileName,PathName] = uigetfile({[ini_path,'/*.tif']},'MultiSelect','on','Select .TIF files');
                    ini_path = PathName;
                    
                    % stop if canceled
                    if isnumeric(FileName)
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
                    if isempty(file_names)
                        file_names = file_names_channel;
                    else
                        % check that the number of files matches
                        if size(file_names,1)~= size(file_names_channel,1)
                            error('Different number of files defined as for previous channels')
                        end
                        file_names = [file_names,file_names_channel];
                    end
                end
            end
        end
        function delete(obj)
            % close all the tif files
            for fileNum = 1:numel(obj.tiff_info)
                obj.tiff_info(fileNum).tifflib.close;
            end
            
        end
        
        %% Correct zig-zag in x pixels
        function set_x_correction(obj,value)
            % determine how the x and y values are shifted
            obj.x_correction = round(value);
        end      
        function img = do_x_correction(obj,img,fill_values)
            % determine how the x and y values are shifted
            if obj.x_correction>0
                img(1:2:end,1+obj.x_correction:end,:) = img(1:2:end,1:end-obj.x_correction,:);
                img(1:2:end,1:obj.x_correction-1,:) = fill_values;
            end
            if obj.x_correction<0
                img(1:2:end,1:end+obj.x_correction,:) = img(1:2:end,1-obj.x_correction:end,:);
                img(1:2:end,end+obj.x_correction+1:end,:) = fill_values;
            end
        end
        
        %% Clip frames to remove shot noise
        function do_clip(obj,test,clip_value)
            % the third argument is to manually set the clip value.
            if ~test
                obj.clip_value = [];
            else
                if nargin==3 && numel(clip_value)==obj.num_channels
                    obj.clip_value = clip_value(:);
                else
                    tmp = obj.mean + 3*obj.variance;
                    obj.clip_value = min(reshape(tmp,[obj.width*obj.height obj.num_channels]),[],1)';
                end
            end
        end
        
        %% Read stack and file properties
        % NOTES:
        % - When reading the stacks the current alignment is used
        % - The properties of the files are only calculated once they are
        % required. The functions used for the calculation are below
        
        function output_image = read_stack(obj,stack_original,channel_to_extract)
            % read a stack from the files.
            % Input:
            % stack_original = stack # from the total stacks (corresponding stack in file is calculated)
            % chanel_to_extract = from which channel the stack is obtained (def.1)
            % fill_values = in case stack is aligned, how are empty spaces filled
            
            % remove warnings
            warning off
            
            % if not given use object defined channel to extract
            if nargin<3
                channel_to_extract = 1;
            end
            if nargin<2
                stack_original = 1;
            end
            
            % Check the corresponding file for the stack and the
            % stack number within
            file_number = find(stack_original>=obj.stack_to_file(:,2) & stack_original<=obj.stack_to_file(:,3));
            if isempty(file_number)
                disp('Required stack is not available. Returning empty stack.')
                return
            end
            stack = stack_original - (obj.stack_to_file(file_number,2)-1);
            
            % Go to stack and read
            obj.tiff_info(file_number,channel_to_extract).tifflib.setDirectory(stack);
            img = double(obj.tiff_info(file_number,channel_to_extract).tifflib.read);
            % if multi-dimensional, use only first
            img = img(:,:,1);
            img = obj.do_x_correction(img,NaN);
            
            % clip if value exist
            if ~isempty(obj.clip_value)
               img(img>obj.clip_value(channel_to_extract)) = obj.clip_value(channel_to_extract); 
            end
            
            % Shift image for alignment
            if sum(abs(obj.img_T(stack_original,:)))==0
                output_image = img;
            else
                output_image = nan(size(img));
                output_image(...
                    (1+max(0,obj.img_T(stack_original,1))):(size(img,1)+min(0,obj.img_T(stack_original,1))) ,...
                    (1+max(0,obj.img_T(stack_original,2))):(size(img,2)+min(0,obj.img_T(stack_original,2))) ) ...
                    = img( (1-min(0,obj.img_T(stack_original,1))):(size(img,1)-max(0,obj.img_T(stack_original,1))) ,...
                    (1-min(0,obj.img_T(stack_original,2))):(size(img,2)-max(0,obj.img_T(stack_original,2))) );
            end
        end
        
        function img_output = mean(obj)
            if isempty(obj.img_mean);obj.properties_calculate(false);end
            img_output = obj.do_x_correction(obj.img_mean,NaN);
        end
        function img_output = variance(obj)
            if isempty(obj.img_variance);obj.properties_calculate(false);end
            img_output =  obj.do_x_correction(obj.img_variance,NaN);
        end
        function img_output = maximum(obj)
            if isempty(obj.img_maximum);obj.properties_calculate(false);end
            img_output =  obj.do_x_correction(obj.img_maximum,NaN);
            if ~isempty(obj.clip_value)
                for ii=1:obj.num_channels
                    tmp = img_output(:,:,ii);
                    tmp(tmp>obj.clip_value(ii)) = obj.clip_value(ii);
                    img_output(:,:,ii) = tmp;
                end
            end
        end
        function img_output = minimum(obj)
            if isempty(obj.img_minimum);obj.properties_calculate(false);end
            img_output =  obj.do_x_correction(obj.img_minimum,NaN);
        end
        function img_output = correlation(obj)
            if isempty(obj.img_correlation);obj.correlation_calculate(false);end
            img_output =  obj.do_x_correction(obj.img_correlation,NaN);
        end
         function img_output = averagedMaximum(obj)
            if isempty(obj.img_averagedMaximum);obj.chunkyMax_calculate(false);end
            img_output =  obj.do_x_correction(obj.img_averagedMaximum,NaN);
        end
        function value = maxValue(obj)
            if isempty(obj.img_maxValue);obj.properties_calculate(false);end
            value =  obj.img_maxValue;
            if ~isempty(obj.clip_value)
                value(value(:)>obj.clip_value) = obj.clip_value(value(:)>obj.clip_value);
            end
        end
        function value = minValue(obj)
            if isempty(obj.img_minValue);obj.properties_calculate(false);end
            value =  obj.img_minValue;
        end
        function T = xyShift(obj)
            T = obj.img_T;
        end
        
        %% Alignment functions
        % While the alignment is calculated the image properties are also
        % obtained (mean, std, etc) of ALL the files together (not
        % individually). To store this properties check the "meta data"
        % methods below
        
        function alignment_set(obj,T)
            % set alignment array externally
            
            % check dimensions
            if ndims(T)>2 || size(T,2)~=2
                error('Input alignment array should have 2D dimensions [frames x 2] ')
            end
            if size(T,1)~= obj.num_stacks
               error('Length of alignment array doesn''t match number of frames') 
            end
            
            % delete previous properties
            obj.metadata_delete
            
            % add alignment array
            obj.img_T = T;
            
        end
        function alignment_calculate(obj,alignment_channel,show_msg,max_shift)
            if nargin==1
                alignment_channel = obj.alignment_channel;
                show_msg = true;
                max_shift = NaN;
            elseif nargin ==2
                show_msg = true;
                max_shift = NaN;
                if alignment_channel>obj.num_channels || alignment_channel<=0
                    error(['The total number of channels is ',num2str(obj.num_channels)])
                end
            elseif nargin==3
                max_shift = NaN;
            end
            
            % erase previous alignment
            obj.img_T = zeros(obj.num_stacks,2);
            
            % calculate alignment
            if show_msg
                h_box = msgbox(['Calculating alignment using channel ',num2str(alignment_channel),'...'],'tiff_reader');
                drawnow
            end
            r = sbxalign(obj,1:obj.num_stacks,alignment_channel,max_shift);
            if show_msg
                close(h_box)
            end

            % set file properties (from all files together!)
            obj.img_mean=[];obj.img_variance=[];obj.img_minimum=[];obj.img_maximum=[];
            obj.img_maxValue=[];obj.img_minValue=[];img_correlation=[];
            for ch=1:obj.num_channels
                obj.img_mean = cat(3,obj.img_mean,r.p{1,ch});
                obj.img_variance = cat(3,obj.img_variance,r.p{2,ch}/(r.n-1));
                obj.img_minimum = cat(3,obj.img_minimum,r.p{3,ch});
                obj.img_maximum = cat(3,obj.img_maximum,r.p{4,ch});
                obj.img_minValue = [obj.img_minValue;min(min(r.p{3,ch}))];
                obj.img_maxValue = [obj.img_maxValue;max(max(r.p{4,ch}))];
            end
            obj.img_T = r.S;
            
            % save alignment channel used
            obj.alignment_channel = alignment_channel;
            
            % delete set properties, as they might have changed with the
            % alignment
            obj.set_mean = cell(size(obj.tiff_info));
            obj.set_variance = cell(size(obj.tiff_info));
            obj.set_maximum = cell(size(obj.tiff_info));
            obj.set_minimum = cell(size(obj.tiff_info));
            obj.set_minValue = cell(size(obj.tiff_info));
            obj.set_maxValue = cell(size(obj.tiff_info));
            obj.set_pix_prod = cell(size(obj.tiff_info,1),1); % only first channel!
            obj.set_averagedMaximum = cell(size(obj.tiff_info,1),1); % only first channel!
            
            % if there is only one file for each channel, fill set
            % properties
            if size(obj.tiff_info,1)==1
                for chnum = 1:size(obj.tiff_info,2)
                    obj.set_mean{1,chnum} = obj.img_mean(:,:,chnum);
                    obj.set_variance{1,chnum} = obj.img_variance(:,:,chnum);
                    obj.set_maximum{1,chnum} = obj.img_maximum(:,:,chnum);
                    obj.set_minimum{1,chnum} = obj.img_minimum(:,:,chnum);
                    obj.set_minValue{1,chnum} = obj.img_minValue(chnum);
                    obj.set_maxValue{1,chnum} = obj.img_maxValue(chnum);
                end
            end
            
        end
        function alignment_reference(obj,alignment_channel)
           if nargin==1
               alignment_channel = obj.alignment_channel;        
           end
            
           % open file(s) class
           obj_h = tiff_reader;
           if obj_h.num_channels>1
              disp('The reference can have only one channel.')
              return
           end
           % get reference
           ref_img = obj_h.mean;
           clear obj_h;
           
           % check that size of the reference image is the same as the
           % opened file
           if ~(size(ref_img,1)==obj.height && size(ref_img,2)==obj.width)
               disp('WARNING: Reference dimensions do not match current files')
               disp('         Fitting to set dimensions...')
               
               dx = obj.width - size(ref_img,2);
               dy = obj.height - size(ref_img,1);
               
               if dx>0
                   ref_img = cat(2,zeros(size(ref_img,1),floor(dx/2)),...
                       ref_img,...
                       zeros(size(ref_img,1),ceil(dx/2)));
               elseif dx<0
                   ref_img =  ref_img(:,max(1,floor(abs(dx)/2)):size(ref_img,2)-ceil(abs(dx)/2));
               end
               
               if dy>0
                   ref_img = cat(1,zeros(floor(dy/2),size(ref_img,2)),...
                       ref_img,...
                       zeros(ceil(dy/2),size(ref_img,2)));
               elseif dy<0
                   ref_img =  ref_img(max(1,floor(abs(dy)/2)):size(ref_img,1)-ceil(abs(dy)/2),:);
               end
           end
           
           % use the image to align current mean
           current_mean = obj.mean;
           [shift_y,shift_x] = fftalign(current_mean(:,:,alignment_channel),ref_img);
           
           % display results and shift alignment
           disp(['X shift = ',num2str(shift_x),' Y shift = ',num2str(shift_y)])
           obj.img_T = bsxfun(@plus,obj.img_T,[shift_y,shift_x]);

            % save alignment channel used
           obj.alignment_channel = alignment_channel; 
        end                  
        
        %% Properties functions
        % The properties of the class are obtained either by combining
        % the info of the individual files (e.g. set_mean) or by
        % calculating them frame by frame.
        
        function properties_calculate(obj,do_all,show_msg)
            % do_all = logical, if all files are re-calculated or only the
            % missing ones
            if nargin == 1
                do_all = true;
                show_msg = true;
            elseif nargin ==2
                show_msg = true;
            end
            
            if show_msg
                h_box = msgbox('Calculating stack properties...','tiff_reader');
                drawnow
            end
            
            % calculate all or missing file properties
            for ii = 1:size(obj.tiff_info,1)
                for jj = 1:size(obj.tiff_info,2)
                    if isempty(obj.set_mean{ii,jj}) || do_all
                        
                        % make running averages for the file
                        ch_mean=zeros(obj.height,obj.width);
                        set_M2=zeros(obj.height,obj.width);
                        ch_minimum=zeros(obj.height,obj.width);
                        ch_maximum=zeros(obj.height,obj.width);
                        
                        stacks2use = obj.stack_to_file(ii,2):obj.stack_to_file(ii,3);
                        for stack_num = 1:length(stacks2use)
                            stack = obj.read_stack(stacks2use(stack_num),jj);
                            ch_minimum=min(cat(3,ch_minimum,stack),[],3);
                            ch_maximum=max(cat(3,ch_maximum,stack),[],3);
                            % update mean and variance online (see: B. P. Welford (1962))
                            delta = stack - ch_mean;
                            ch_mean = ch_mean + delta/stack_num;
                            set_M2 = set_M2 + delta.*(stack - ch_mean);
                        end
                        ch_variance = set_M2/(length(stacks2use)-1);
                        
                        % add to data
                        obj.set_mean{ii,jj} = ch_mean;
                        obj.set_variance{ii,jj} = ch_variance;
                        obj.set_maximum{ii,jj} = ch_maximum;
                        obj.set_minimum{ii,jj} = ch_minimum;
                        obj.set_maxValue{ii,jj} = max(ch_maximum(:));
                        obj.set_minValue{ii,jj} = min(ch_minimum(:));
                    end
                end
            end
            
            % combine single file properties to obtain total properties
            obj.img_mean = zeros(obj.height,obj.width,obj.num_channels);
            obj.img_variance = zeros(obj.height,obj.width,obj.num_channels);
            obj.img_maximum = zeros(obj.height,obj.width,obj.num_channels);
            obj.img_minimum = zeros(obj.height,obj.width,obj.num_channels);
            obj.img_minValue = zeros(obj.num_channels,1);
            obj.img_maxValue = zeros(obj.num_channels,1);
            
            % go through channels
            for channel = 1:obj.num_channels
                
                    % combine data sets
                    ch_minimum = obj.set_minimum{1,channel};
                    ch_maximum = obj.set_maximum{1,channel};
                    ch_variance = obj.set_variance{1,channel};
                    ch_mean = obj.set_mean{1,channel};
                    sum_stacks = (obj.stack_to_file(1,3)-obj.stack_to_file(1,2)+1);
                    
                    for file_num = 2:size(obj.tiff_info,1)
                        % min and max
                        ch_minimum=min(cat(3,obj.set_minimum{file_num,channel},ch_minimum),[],3);
                        ch_maximum=max(cat(3,obj.set_maximum{file_num,channel},ch_maximum),[],3);
                        % get mean and variance of new set
                        Vy = obj.set_variance{file_num,channel};
                        ny =  (obj.stack_to_file(file_num,3)-obj.stack_to_file(file_num,2)+1);
                        Y = obj.set_mean{file_num,channel};
                        % calculate total variance
                        weighted_mean = (ch_mean*sum_stacks + Y*ny)/(sum_stacks+ny);
                        ch_variance = ( sum_stacks*(ch_variance+(ch_mean-weighted_mean).^2) + ny*(Vy+(Y-weighted_mean).^2) )/(sum_stacks + ny);  
                        % update mean and total stacks
                        ch_mean = weighted_mean;
                        sum_stacks = sum_stacks + ny;
                    end
                    
                % add to channel parameters
                obj.img_mean(:,:,channel) = ch_mean;
                obj.img_variance(:,:,channel) = ch_variance;
                obj.img_maximum(:,:,channel) = ch_maximum;
                obj.img_minimum(:,:,channel) = ch_minimum;
                obj.img_minValue(channel) = min(ch_minimum(:));
                obj.img_maxValue(channel) = max(ch_maximum(:));
            end
            
            % close the message box
            if show_msg
                close(h_box)
            end
        end
        function correlation_calculate(obj,do_all,show_msg)
            % Formula used for X and Y pixels:
            % corr(X,Y) = ( E[XY]-E[X]E[Y] )/( std[X]std[Y] )
            if nargin == 1
                do_all = true;
                show_msg = true;
            elseif nargin ==2
                show_msg = true;
            end
            
            % === get sum of multiplications for part E[XY]
            if any(cellfun(@isempty,obj.set_pix_prod))
                obj.pix_product_calculate(do_all,show_msg)
            end      
            % Join products to get correlation between pixels
            E_XY = obj.set_pix_prod{1};
            for ii = 2:size(obj.tiff_info,1)
                E_XY = E_XY +  obj.set_pix_prod{ii};
            end           
            % calculate mean
            E_XY = E_XY./obj.num_stacks;
            
            % === Extract mean and std
            nn_indices = obj.get_nn_pix_ind();
            % E[Z]
            E_X = obj.mean; E_X = E_X(:,:,1); E_X = E_X(:);
            E_Y = zeros(size(E_XY));
            E_Y(~isnan(nn_indices)) = E_X(nn_indices(~isnan(nn_indices)));
            % std[Z]
            std_X = sqrt(obj.variance); std_X = std_X(:,:,1); std_X = std_X(:);
            std_Y = zeros(size(E_XY));
            std_Y(~isnan(nn_indices)) = std_X(nn_indices(~isnan(nn_indices)));
            
            % === Calculate correlation
            % corr(X,Y) = ( E[XY] - E[X]E[Y] )/( std[X]std[Y] )
            corr_values = (E_XY - bsxfun(@times,E_X,E_Y))./bsxfun(@times,std_X,std_Y);
            corr_values(isnan(corr_values)) = 0;
                        
            % == Get image with average NN correlation
            obj.img_correlation =  zeros(obj.height,obj.width);
            obj.img_correlation(:) = sum(corr_values(:,2:end),2)./(sum(~isnan(nn_indices(:,2:end)),2));            
        end
        function nn_indices = get_nn_pix_ind(obj)
            % Neighboring pixel indexing in a +/-1 window
            range = 1;
            
            [idx_y,idx_x]=ind2sub([obj.height obj.width],1:obj.height*obj.width);
            nn_indices = NaN(obj.height*obj.width,(2*range+1)^2);
            
            nn_idx = 0;
            for move_y = [0 -range:1:-1 1:range]
                for move_x = [0 -range:1:-1 1:range]
                    nn_idx = nn_idx +1;
                    
                    idx_y_shift = idx_y + move_y;
                    idx_x_shift = idx_x + move_x;
                    
                    in_frame = (idx_y_shift<=obj.height & idx_y_shift>=1 & idx_x_shift<=obj.width & idx_x_shift>=1);
                    nn_indices(in_frame,nn_idx) = sub2ind([obj.height obj.width],idx_y_shift(in_frame),idx_x_shift(in_frame));
                end
            end
        end
        function pix_product_calculate(obj,do_all,show_msg)
            % calculate the correlation between neighboring pixels. For
            % each set the product of the values are stored.
            % corr(X,Y) = ( E[XY]-E[X]E[Y] )/( std[X]std[Y] )
            % Here the sum E[XY] is calculated
            
            if show_msg
                h_box = msgbox('Calculating neighbor pixel correlation...','tiff_reader');
                drawnow
            end
            
            % === Go over frames and sum pixel products
            nn_indices = obj.get_nn_pix_ind();
            Y = NaN(size(nn_indices));
            for ii = 1:size(obj.tiff_info,1)
                if isempty(obj.set_pix_prod{ii}) || do_all
                    obj.set_pix_prod{ii} = zeros(size(nn_indices));
                    stacks2use = obj.stack_to_file(ii,2):obj.stack_to_file(ii,3);
                    for stack_num = 1:length(stacks2use)
                        % read first channel
                        X = obj.read_stack(stacks2use(stack_num),1);
                        Y(~isnan(nn_indices)) = X(nn_indices(~isnan(nn_indices)));
                        % sum XY
                        obj.set_pix_prod{ii} = obj.set_pix_prod{ii} + bsxfun(@times,X(:),Y);
                    end
                end
            end
            
            % ===  close the message box
            if show_msg
                close(h_box)
            end
        end
        function chunkyMax_calculate(obj,show_msg)
            if nargin == 1
                show_msg = true;
            end
            
            % Calculate averaged maximum
            mean_samples = 25;
            chunkyMax = zeros(obj.height,obj.width);
            
            % 3) Make starting average
            rolling_template = zeros(obj.height,obj.width);
            frame_count = zeros(obj.height,obj.width);
            
            for frame_ii=1:mean_samples
                S = obj.read_stack(frame_ii,1);
                rolling_template(~isnan(S)) = rolling_template(~isnan(S)) + S(~isnan(S));
                frame_count(~isnan(S)) = frame_count(~isnan(S)) + 1;
            end
            
            % 4) Go through the frames, 
            
            inform_stack = ceil(linspace(1,obj.num_stacks,11));
            for frame_ii=1:obj.num_stacks
                
                % inform how we are doing
                [~,ind] = intersect(inform_stack,frame_ii);
                if ~isempty(ind) && show_msg
                    disp(['     ',num2str((ind-1)*10),'%'])
                end
                
                % move mean
                template_range = frame_ii + mean_samples*[-1 1];
                
                % remove first image in range
                if template_range(1)>0
                    S = obj.read_stack(template_range(1),1);
                    rolling_template(~isnan(S)) = rolling_template(~isnan(S)) - S(~isnan(S));
                    frame_count(~isnan(S)) = frame_count(~isnan(S)) - 1;
                end
                % add new image to range
                if template_range(2)<=obj.num_stacks
                    S = obj.read_stack(template_range(2),1);
                    rolling_template(~isnan(S)) = rolling_template(~isnan(S)) + S(~isnan(S));
                    frame_count(~isnan(S)) = frame_count(~isnan(S)) + 1;
                end
                
                % update max
                chunkyMax = max(cat(3,chunkyMax,rolling_template./frame_count),[],3);
                
                % check if it is time to save and restart chunkyMax
                file_ind = find(ismember(obj.stack_to_file(:,3),frame_ii));
                if ~isempty(file_ind)
                    % save and restart
                    obj.set_averagedMaximum{file_ind} = chunkyMax;
                    chunkyMax = zeros(obj.height,obj.width);
                end
            end
            
            % join to get image
            obj.img_averagedMaximum = zeros(obj.height,obj.width);
            for file_ind=1:size(obj.stack_to_file,1)
                obj.img_averagedMaximum = max(cat(3,obj.img_averagedMaximum,...
                    obj.set_averagedMaximum{file_ind}),[],3);
            end
            
        end
        %% Meta data: read and save
        % In the file meta_data.mat the properties and alignment of the
        % images are stored in a cell.
        
        function metadata_delete(obj)
            % properties
            obj.img_mean = [];
            obj.img_variance = [];
            obj.img_maximum = [];
            obj.img_minimum = [];
            obj.img_minValue = [];
            obj.img_maxValue = [];
            obj.img_correlation = [];
            
            obj.set_mean = cell(size(obj.tiff_info));
            obj.set_variance = cell(size(obj.tiff_info));
            obj.set_maximum = cell(size(obj.tiff_info));
            obj.set_averagedMaximum = cell(size(obj.tiff_info,1),1);
            obj.set_minimum = cell(size(obj.tiff_info));
            obj.set_minValue = cell(size(obj.tiff_info));
            obj.set_maxValue = cell(size(obj.tiff_info));
            obj.set_pix_prod = cell(size(obj.tiff_info,1),1);
            
            % alignment
            obj.img_T = zeros(obj.num_stacks,2);
        end
        function metadata_load(obj)
            
            % delete previous properties
            obj.metadata_delete
            
            % open meta data file
            separators=strfind(obj.tiff_info(1).Filename,filesep);
            path=obj.tiff_info(1).Filename(1:separators(end));
            if exist([path,'meta_data.mat'],'file')
                load([path,'meta_data.mat'],'stack_properties');
                if ~exist('stack_properties','var')
                    disp('Stack properties not found in meta_data file.')
                    return
                end
            else
                disp('File meta_data.mat not found in folder!')
                return
            end
            
            % === PROPERTIES
            
            % fill array
            for file_num=1:size(obj.tiff_info,1)
                for ch_num=1:size(obj.tiff_info,2)
                    fname = obj.tiff_info(file_num,ch_num).Filename(separators(end)+1:end);
                    [~,ind_t]=intersect(stack_properties(:,1),fname);
                    if isempty(ind_t)
                        disp('Not all files are stored in meta_data.mat!')
                        return
                    end
                    obj.set_mean{file_num,ch_num} = stack_properties{ind_t,2}.mean;
                    obj.set_variance{file_num,ch_num} = stack_properties{ind_t,2}.variance;
                    obj.set_maximum{file_num,ch_num} = stack_properties{ind_t,2}.maximum;
                    obj.set_minimum{file_num,ch_num} = stack_properties{ind_t,2}.minimum;
                    obj.set_minValue{file_num,ch_num} = stack_properties{ind_t,2}.minValue;
                    obj.set_maxValue{file_num,ch_num} = stack_properties{ind_t,2}.maxValue;
                    % if any of the files has a correction, use it
                    if isfield(stack_properties{ind_t,2},'x_correction')
                        obj.set_x_correction(stack_properties{ind_t,2}.x_correction);
                    end
                    if isfield(stack_properties{ind_t,2},'clip_value')
                        obj.clip_value = stack_properties{ind_t,2}.clip_value;
                    end
                    % only for first channel
                    if ch_num==1
                        % get pixel product
                        if isfield(stack_properties{ind_t,2},'set_pix_prod')
                            obj.set_pix_prod{file_num} = stack_properties{ind_t,2}.set_pix_prod;
                        end
                        % get averaged maximum
                        if isfield(stack_properties{ind_t,2},'averagedMaximum')
                            obj.set_averagedMaximum{file_num} = stack_properties{ind_t,2}.averagedMaximum;
                        end
                    end
                end
            end
            % calculate correlation image if all pixel products calculated
            if ~any(cellfun(@isempty,obj.set_pix_prod))
                obj.correlation_calculate(false,false)
            end
            % calculate averagedMaximum image if all files are calculated
            if ~any(cellfun(@isempty,obj.set_averagedMaximum))
                % join to get image
                obj.img_averagedMaximum = zeros(obj.height,obj.width);
                for file_ind=1:size(obj.stack_to_file,1)
                    obj.img_averagedMaximum = max(cat(3,obj.img_averagedMaximum,...
                        obj.set_averagedMaximum{file_ind}),[],3);
                end
            end
            % === ALIGNMENT
            
            %------ !! Back compatibility issue!!
            % If 3rd column of stack properties is a cell, then it is the
            % old verion of meta_data and the alignment data is stored in
            % alignment_data. Convert to new version.
            if iscell(stack_properties{1,3})
                load([path,'meta_data.mat'],'alignment_data');
                if ~exist('alignment_data','var')
                    for stack_ii = 1:size(stack_properties,1)
                        stack_properties{stack_ii,3} = [];
                    end
                else
                    for stack_ii = 1:size(stack_properties,1)
                        [~,ind_t]=intersect(alignment_data(:,1),{stack_properties{stack_ii,1}(1:end-8)}); %_ch1.tif
                        if ~isempty(ind_t)
                            stack_properties{stack_ii,3} = alignment_data{ind_t,2};
                        end
                    end
                end
            end
            %------ !! END: Back compatibility issue!!
            
            % check that all files in alignment channel are stored in the file
            idx_file = zeros(size(obj.tiff_info,1),1);
            for file_num=1:size(obj.tiff_info,1)
                fname = obj.tiff_info(file_num,obj.alignment_channel).Filename(separators(end)+1:end);
                [~,ind_t]=intersect(stack_properties(:,1),fname);
                if isempty(ind_t)
                    disp(['File ',fname,' not found in alignment_data.mat!'])
                    return
                end
                idx_file(file_num) = ind_t;
            end
            
            % combine transform array
            T = [];
            for file_num=1:size(obj.tiff_info,1)
                Tprime = stack_properties{idx_file(file_num),3};
                if isempty(Tprime)
                    Tprime = zeros(obj.tiff_info(file_num,1).Stacks,2);
                end
                T = [T;Tprime];
            end
            obj.img_T = T;

            % === 
            % join
            obj.properties_calculate(false,false)
            disp('Properties and alignment loaded.')
            
        end
        function metadata_save(obj)
            % Generates or complements a Nx3 cell structure with the
            % properties and alignment of all files in the class
            
             % calculate only the missing properties
            obj.properties_calculate(false)
            
            % fill a cell with [file name, properties structure, all files in the list]
            stack_properties = cell(numel(obj.tiff_info),3);
            stack_ind = 0;
            for file_num=1:size(obj.tiff_info,1)
                for ch_num=1:size(obj.tiff_info,2)
                    
                    separators=strfind(obj.tiff_info(file_num,ch_num).Filename,filesep);
                    fname = obj.tiff_info(file_num,ch_num).Filename(separators(end)+1:end);
                    
                    stack_ind = stack_ind + 1;
                    stack_properties{stack_ind,1} = fname;
                    
                    stack_properties{stack_ind,2}.mean = obj.set_mean{file_num,ch_num};
                    stack_properties{stack_ind,2}.variance = obj.set_variance{file_num,ch_num};
                    stack_properties{stack_ind,2}.minimum = obj.set_minimum{file_num,ch_num};
                    stack_properties{stack_ind,2}.maximum = obj.set_maximum{file_num,ch_num};
                    stack_properties{stack_ind,2}.minValue = obj.set_minValue{file_num,ch_num};
                    stack_properties{stack_ind,2}.maxValue = obj.set_maxValue{file_num,ch_num};
                    stack_properties{stack_ind,2}.x_correction = obj.x_correction;
                    stack_properties{stack_ind,2}.clip_value = obj.clip_value;
                    
                    stack_properties{stack_ind,3} = obj.img_T(obj.stack_to_file(file_num,2):obj.stack_to_file(file_num,3),:);
                
                    if ch_num==1
                        stack_properties{stack_ind,2}.set_pix_prod = obj.set_pix_prod{file_num};
                        stack_properties{stack_ind,2}.averagedMaximum = obj.set_averagedMaximum{file_num};
                    end
                end
            end
            
            % open file to save
            tmp=strfind(obj.tiff_info(1).Filename,filesep);
            path=obj.tiff_info(1).Filename(1:tmp(end));
            
            % if exists, load and combine (overwrite info if files are doubled)
            if exist([path,'meta_data.mat'],'file')
                prior = load([path,'meta_data.mat'],'stack_properties');
                % delete repeated files in stack properties
                [~,stack_ind] = intersect(prior.stack_properties(:,1),stack_properties(:,1));
                prior.stack_properties(stack_ind,:) = [];
                % join
                stack_properties = [stack_properties;prior.stack_properties];               
                % save
                save([path,'meta_data.mat'],'stack_properties','-append')
            else
                % save directly
                save([path,'meta_data.mat'],'stack_properties')
            end
            
            disp('Properties saved')
            
            
        end
        
        %% Export function
        
        function export(obj)
                       
            % make general tag
            tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
            tagstruct.ImageLength = obj.height;
            tagstruct.ImageWidth = obj.width;
            tagstruct.RowsPerStrip = obj.tiff_info(1,1).tifflib.getTag('RowsPerStrip');
            tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
            tagstruct.SamplesPerPixel = 1;
            tagstruct.Compression = Tiff.Compression.None;
            
            data_class = class(obj.tiff_info(1,1).tifflib.read);
            switch data_class
                case {'uint8', 'uint16', 'uint32'}
                    tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
                case {'int8', 'int16', 'int32'}
                    tagstruct.SampleFormat = Tiff.SampleFormat.Int;
                case {'single', 'double', 'uint64', 'int64'}
                    tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
                otherwise
                    error('Error in determining SampleFormat');
            end 
            switch data_class
                case {'uint8', 'int8'}
                    tagstruct.BitsPerSample = 8;
                case {'uint16', 'int16'}
                    tagstruct.BitsPerSample = 16;
                case {'uint32', 'int32'}
                    tagstruct.BitsPerSample = 32;
                case {'single'}
                    tagstruct.BitsPerSample = 32;
                case {'double', 'uint64', 'int64'}
                    tagstruct.BitsPerSample = 64;
                otherwise
                    error('Error in determining BitsPerSample');
            end
            % transform to integers function
            eval(['trfun = @',data_class,';'])
            
            
            for file_num = 1:size(obj.tiff_info,1)
                for channel_num = 1:size(obj.tiff_info,2)
                    
                    disp(['Exporting set ',num2str(file_num),' channel ',num2str(channel_num),' ...'])
                    
                    % delete previous
                    [pathstr, fname, fext] = fileparts(obj.tiff_info(file_num,channel_num).Filename);
                    file_name = [pathstr,filesep,'export_',fname,fext];
                    if exist(file_name,'file')
                        delete(file_name)
                    end
                   
                    % make new file
                    tfile = Tiff(file_name, 'w');
                   
                    % add stack
                    stacks_to_use = obj.stack_to_file(file_num,2):obj.stack_to_file(file_num,3);
                    for stack = stacks_to_use
                        % get image
                        Im = trfun(obj.read_stack(stack,channel_num));
                        % save
                        tfile.setTag(tagstruct);
                        tfile.write(Im);
                        % expand
                        if stack ~= stacks_to_use(end)
                            tfile.writeDirectory();
                        end
                    end
                    
                    % close file
                    tfile.close();
                    
                    disp(['   ... done.'])
                end
            end
        end
    end
    
end

%% FUNCTIONS TO CALCULATE X-Y ALIGNMENT AND PROPERTIES
% based on:
% https://scanbox.wordpress.com/2014/03/20/recursive-image-alignment-and-statistics/

function r = sbxalign(obj,stacks_to_do,channel_align,max_shift)

if length(stacks_to_do)==1
    % if only one stack, assign initial values to parameters
    r.n = 1;
    r.S = [0 0];
    r.p=cell(4,obj.num_channels);
    for channel=1:obj.num_channels
        A = obj.read_stack(stacks_to_do,channel);
        r.p{1,channel} = A;
        r.p{2,channel} = zeros(size(A));
        r.p{3,channel} = A;
        r.p{4,channel} = A;
    end
else
    % split into two groups and run again (recursively)
    r_input = sbxalign(obj,stacks_to_do(1:floor(end/2)),channel_align,max_shift);
    r_goal = sbxalign(obj,stacks_to_do(floor(end/2)+1:end),channel_align,max_shift);
    
    % use the average of the selected channel to get alignment
    [shift_y,shift_x] = fftalign(r_input.p{1,channel_align},r_goal.p{1,channel_align},max_shift);
    
    % combine transforms and number of stacks
    r.S = [ones(size(r_input.S,1),1)*[shift_y shift_x] + r_input.S ; r_goal.S];
    r.n = r_input.n+r_goal.n;
    
    % shift and update properties for each channel
    for channel=1:obj.num_channels
        % shift
        if abs(shift_y)+abs(shift_x)>0
            for parameter=1:4
                r_input.p{parameter,channel} = shift_image(r_input.p{parameter,channel},[shift_y shift_x]);
            end
        end
        % online update (see: Chan, Tony F.; Golub, Gene H.; LeVeque, Randall J. (1979))
        r.p{1,channel} = (r_input.n*r_input.p{1,channel} + r_goal.n*r_goal.p{1,channel})/r.n;
        r.p{2,channel} = r_input.p{2,channel} + r_goal.p{2,channel} +...
            (r_goal.p{1,channel}-r_input.p{1,channel}).^2*r_input.n*r_goal.n/r.n;
        r.p{3,channel} = min(cat(3,r_input.p{3,channel},r_goal.p{3,channel}),[],3);
        r.p{4,channel} = max(cat(3,r_input.p{4,channel},r_goal.p{4,channel}),[],3);
    end
end
end

function [shift_y,shift_x] = fftalign(img_input,img_goal,max_shift)

% cut edges with NaN (from previous transforms) using the center of the
% images as reference
nan_y = isnan(img_input(:,ceil(size(img_input,2)/2)).*img_goal(:,ceil(size(img_input,2)/2)));
nan_x = isnan(img_input(ceil(size(img_input,1)/2),:).*img_goal(ceil(size(img_input,1)/2),:));
img_input = img_input(find(~nan_y),find(~nan_x));
img_goal = img_goal(find(~nan_y),find(~nan_x));

% cut both images to be squared using the minimum size of img_input as reference
N = floor(min(size(img_input))/2);
yidx = floor(size(img_input,1)/2) -N + 1 : floor(size(img_input,1)/2) + N;
xidx = floor(size(img_input,2)/2) -N + 1 : floor(size(img_input,2)/2) + N;
img_input = img_input(yidx,xidx);
img_goal = img_goal(yidx,xidx);

% Calculate cross correlation and find position of maximum
C = fftshift(real(ifft2(fft2(img_input).*fft2(rot90(img_goal,2)))));
if ~isnan(max_shift)
    % restrict search shifts
    [xi,yi] = meshgrid((1:2*N)-N+1);
    mask = (xi.^2+yi.^2) < (max_shift^2);
    C = C.*mask;
end
[ind_y,ind_x] = find(C == max(C(:)));

% Use the position of the maximum relative to center to determine shift
shift_y = N - ind_y;
shift_x = N - ind_x;

end

function img_output = shift_image(img_input,shift)
% shift = [y x] and fill empty spaces with NaN
img_output = NaN*zeros(size(img_input));
img_output( (1+max(0,shift(1))):(size(img_input,1)+min(0,shift(1))) , (1+max(0,shift(2))):(size(img_input,2)+min(0,shift(2))) ) = ...
    img_input( (1-min(0,shift(1))):(size(img_input,1)-max(0,shift(1))) , (1-min(0,shift(2))):(size(img_input,2)-max(0,shift(2))) );
end

