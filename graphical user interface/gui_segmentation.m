function gui_segmentation(start_type,start_info)
% Graphical user interface for ROI definition of 2photon microscopy data.
% - Handle tif files
% - Define and handle regions of interest (ROI)
%
% Required files:
%
% tiff_reader.m = class that opens a set of tif files as handle and reads it
%
% segmentation_gui_draw.m = makes the GUI
% segmentation_gui_draw.mat = used to make the GUI
% segmentation_gui_draw.fig = to re-compile segmentation_gui_draw.m (open in GUIDE and export)
%
% segmentation_ROI_detect.m = semi automatic ROI detection
%
% INPUT:
% start_type : 'data' or 'ROI' to load data or load saved ROIs+data
%
% Created by
% Juan Daniel Florez Weidinger - chepe@nld.ds.mpg.de

%% DEFINE DATA STRUCTURES AND DRAW GUI

if nargin<2
    ini_path =  [pwd,filesep];
    start_type = 'gui';
else
    
    switch start_type
        
        case 'gui'
            ini_path = start_info;
            
        case 'file_structure'
            file_structure = start_info;
            ini_path = file_structure{1};
            sep_idx = strfind(ini_path,filesep);
            ini_path = ini_path(1:sep_idx(end));
            
        case 'exported'
            ini_path =  start_info;
    end
    
end

% === GUI properties:
gui_info = [];

gui_info.control.path = ini_path;
gui_info.control.roi_export_ready = false;

gui_info.data.channel = [false;false;false];
gui_info.data.threshold = [0 0.5;0 1;0 1];
gui_info.data.frame = 1;
gui_info.data.do_clip = false;
gui_info.data.x_correction = 0;
gui_info.data.frame_rate = [];

gui_info.roi.label_count = 0;
gui_info.roi.threshold = 0.25;
gui_info.roi.percent_in = 0.25;
gui_info.roi.percent_out = 0.35;
gui_info.roi.type = 'Neuron';
gui_info.roi.list = {};
gui_info.roi.selected = 0;
gui_info.roi.new_active = false;
gui_info.roi.XYshif_active = false;
gui_info.roi.radius = 10;

gui_info.image.display = 'Average';
gui_info.image.handle = [];
gui_info.image.filth_th = 1;
gui_info.image.roi = true;
gui_info.image.label = false;

% === Data (alignment + images are handled in tiff_reader)
data_info = [];
data_info.handle = [];
data_info.channel_available = [false false false];
data_info.width = [];
data_info.height = [];
data_info.num_stacks = [];
data_info.filth = [];

% === Image to display
image_info = [];
image_info.extracted = [];
image_info.plotted = [];
image_info.xlim = [];
image_info.ylim = [];

% === ROI
roi_info = [];
roi_info.count = 0;
roi_info.ROI_xy = {};
roi_info.ROI_label = [];
roi_info.ROI_type = {};
roi_info.handle_contour = {};
roi_info.handle_label = {};

%% INITIATE

switch start_type
    
    case 'gui'
        % load a data set to start the GUI
        load_data('start')
        % draw the interface with callbacks and get handles
        [hfig,h]=gui_create();
        
        % ---> plot the first image
        prepare_image()
        plot_image()
        plot_ROI()
    
    case 'file_structure'
        % load a data set
        load_data('external',file_structure)
        % draw the interface with callbacks and get handles
        [hfig,h]=gui_create();
        
        % ---> plot the first image
        prepare_image()
        plot_image()
        plot_ROI()
        
    case 'exported'
        disp('Loading exported file...')
        load_exported()
end

gui_resize()

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               V I E W
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------- GUI
    function [hfig,h]=gui_create()
        % Function that draws the GUI. Separated in 2 panels: Data and ROI
        % Returns hfig: handle to figure and h: structure of
        % handles to all interfaces with interface's tag as name
        % ---
        
        % close previous window
        
        h_fig = findobj('Tag','hfig');
        close(h_fig)
        
        % Draw the initial GUI using the results of GUIDE
        % hfig: main figure hangle
        % h   : structure with controls handles
        sourcepath = mfilename('fullpath');
        idx = strfind(sourcepath,filesep);
        hfig = open([sourcepath(1:idx(end)),'gui_segmentation.fig']);
        h = guihandles(hfig);
        
        % General figure
        set(hfig,{'Name','WindowScrollWheelFcn','ResizeFcn'},{'Data segmentation GUI',@gui_scroll,@gui_resize})
        % Axes to display image of microscope
        set(h.D_axes,{'HandleVisibility','NextPlot','xTick','yTick','Box','YDir','ButtonDownFcn','DataAspectRatio',...
            'XLimMode','YLimMode','ZLimMode','XTickMode','YTickMode','ZTickMode'...
            ,'XTickLabelMode','YTickLabelMode','ZTickLabelMode','TickDirMode','ClimMode','Clim','AlimMode','Alim','DrawMode'},...
            {'on','replacechildren',[],[],'on','reverse',@gui_click,[1 1 1],...
            'manual','manual','manual','manual','manual','manual','manual','manual','manual','manual','manual',[0 1],'manual',[0 1],'fast'})
        % image slider
        set(h.D_sld_frames,{'Min','Max','Value','SliderStep','CallBack','Enable'},...
            {1,data_info.num_stacks,1,[1 0.1*data_info.num_stacks]/data_info.num_stacks,@callback_data,'on'})
        
        % Change GUI properties
        % ------------------------------------------------   Context menu
        % ::: Make context menu to align and save data properties
        mh0 = uimenu(hfig,'Label','Data');
        mh1 = uimenu(mh0,'Label','Properties');
        h.D_prop_clear = uimenu(mh1,'Label','Clear','Callback',@callback_uimenu,'Enable','on');
        h.D_prop_calculate = uimenu(mh1,'Label','Calculate','Callback',@callback_uimenu,'Enable','on');
        h.D_prop_internal = uimenu(mh1,'Label','Align frame-by-frame','Callback',@callback_uimenu,'Enable','on');
        h.D_prop_reference = uimenu(mh1,'Label','Align to reference file','Callback',@callback_uimenu,'Enable','on');
        h.D_prop_save = uimenu(mh1,'Label','Save','Callback',@callback_uimenu,'Enable','on');
        if gui_info.data.do_clip
            h.D_clip = uimenu(mh0,'Label','Do clipping','Checked','on','Callback',@callback_uimenu);
        else
            h.D_clip = uimenu(mh0,'Label','Do clipping','Checked','off','Callback',@callback_uimenu);
        end
        h.D_info = uimenu(mh0,'Label','Show info','Callback',@callback_uimenu);
        uimenu(mh0,'Label','Import Data + ROI','CallBack',@load_exported);
        
        % ::: Make context menu for the display
        mh0 = uimenu(hfig,'Label','Display');
        h.D_displayType = uimenu(mh0,'Label','Image');
        h.D_type_average = uimenu(h.D_displayType,'Label','Average','Checked','on','Callback',@callback_uimenu,'Enable','on');
        h.D_type_maximum = uimenu(h.D_displayType,'Label','Maximum','Callback',@callback_uimenu,'Enable','on');
        h.D_type_variance = uimenu(h.D_displayType,'Label','Variance','Callback',@callback_uimenu,'Enable','on');        
        h.D_type_frame = uimenu(h.D_displayType,'Label','Frame','Callback',@callback_uimenu,'Enable','on');
        h.D_type_difference = uimenu(h.D_displayType,'Label','Difference','Callback',@callback_uimenu,'Enable','off');
        h.D_type_correlation = uimenu(h.D_displayType,'Label','Correlation','Callback',@callback_uimenu,'Enable','off');
        h.D_type_chunkyMax = uimenu(h.D_displayType,'Label','Chunky Max','Callback',@callback_uimenu,'Enable','off');
        if sum(data_info.channel_available)>1
            set(h.D_type_difference,'Enable','on')
        end
        if ~isempty(data_info.handle.img_correlation)
            set(h.D_type_correlation,'Enable','on')
        end
        if ~isempty(data_info.handle.img_averagedMaximum)
            set(h.D_type_chunkyMax,'Enable','on')
        end
        mh1 = uimenu(mh0,'Label','ROI');
        h.R_show = uimenu(mh1,'Label','Show contour','Checked','on','Callback',@callback_uimenu);
        h.R_ID = uimenu(mh1,'Label','Show label','Checked','off','Callback',@callback_uimenu);
        h.R_do_zoom =  uimenu(mh0,'Label','Zoom to ROI','checked','off','Callback',@callback_uimenu);
        
        % :: Make context menu for the ROI generation
        mh0 = uimenu(hfig,'Label','ROI');
        uimenu(mh0,'Label','Set thresholds','Callback',{@callback_parameters,'detection_parameters'});
        h.R_invert = uimenu(mh0,'Label','Invert signal','Checked','off','Callback',@callback_uimenu);
        h.R_seedAll = uimenu(mh0,'Label','Re-seed all cells','Enable','off','Callback',@ROI_seed_all);
        
        % :: Make context menu to make figures with ROIs
        h.R_doPlot = uimenu(hfig,'Label','Plot');
        uimenu(h.R_doPlot,'Label','ROI number','Enable','on','Callback',{@ROI_figure,'ID'})
        uimenu(h.R_doPlot,'Label','ROI contour','Enable','off','Callback',{@ROI_figure,'ROI'})
        uimenu(h.R_doPlot,'Label','ROI background','Enable','off','Callback',{@ROI_figure,'Background'})
        
        %----------------------------------------------   Data panel  -> D
        % Channel display properties
        gui_info.data.channel = data_info.channel_available;
        textN = {'green','red','blue'};
        for chN=1:3
            if gui_info.data.channel(chN)
                enable = 'on';
            else
                enable = 'off';
            end
            set(h.(['D_',textN{chN},'_check']),{'CallBack','Value','Enable'},{@callback_data,gui_info.data.channel(chN),enable})
            set(h.(['D_',textN{chN},'_sld_low']),{'CallBack','Enable','Min','Max','Value','SliderStep'},{@callback_data,enable,0,100,gui_info.data.threshold(chN,1)*100,[0.01 0.2]})
            set(h.(['D_',textN{chN},'_sld_high']),{'CallBack','Enable','Min','Max','Value','SliderStep'},{@callback_data,enable,0,100,gui_info.data.threshold(chN,2)*100,[0.01 0.2]})
            set(h.(['D_',textN{chN},'_txt_low']),{'CallBack','Enable','String'},{@callback_data,enable,gui_info.data.threshold(chN,1)*100})
            set(h.(['D_',textN{chN},'_txt_high']),{'CallBack','Enable','String'},{@callback_data,enable,gui_info.data.threshold(chN,2)*100})
        end
        % handles to alignment and properties
        set(h.D_load,'CallBack',@callback_data)
        set(h.D_frame_rate,{'CallBack','String'},{@callback_data,num2str(gui_info.data.frame_rate)})
        % handles to x pixel correction
        set(h.D_xCor_sld,{'CallBack','Min','Max','Value','SliderStep'},{@callback_data,-49,50,gui_info.data.x_correction,[0.01 0.1]})
        set(h.D_xCor_txt,{'CallBack','String'},{@callback_data,num2str(gui_info.data.x_correction)})
        % get initial plot limits from size of image
        image_info.xlim = [1 data_info.width];
        image_info.ylim = [1 data_info.height];
        
        set(h.D_filth,{'CallBack','Value'},{@callback_data,0})
        
        %-------------------------------------------------- ROI panel  -> R
        % ROI list and buttons
        set(h.R_import,{'CallBack'},{@callback_ROI})
        set(h.R_list,{'CallBack','String'},{@callback_ROI,''})
        set(h.R_remove,{'CallBack'},{@callback_ROI})
        set(h.R_sort,{'CallBack'},{@callback_ROI})
        set(h.R_new,{'CallBack','Value'},{@callback_ROI,0})
        set(h.R_reseed,{'CallBack','Enable'},{{@ROI_make_new,'RESEED'},'off'})
        set(h.R_process_export,{'CallBack'},{@callback_ROI})
        set(h.R_plusLabel,{'CallBack'},{@callback_ROI})
        set(h.R_shift,{'CallBack','Value'},{@callback_ROI,0})
        set(h.R_align,'CallBack',@callback_ROI)
        set(h.R_do_ring,'Value',0)
        set(h.R_clear,{'CallBack'},{@callback_ROI})
        % DRAW ROI selection method
        set(h.R_make_type,{'SelectionChangeFcn','SelectedObject'},{@callback_ROI,h.R_draw})
        set(h.R_ROI_type,{'SelectionChangeFcn','SelectedObject'},{@callback_ROI,h.R_type_neuron})
        set(h.R_include,{'Enable'},{'off'})
        set(h.R_exclude,{'Enable'},{'off'})
        set(h.R_overwrite,{'CallBack','Value'},{@callback_ROI,0})       
    end

    function gui_scroll(hObject,eventdata)
        % Function than handles scrolling on the axes. The axes of the
        % mouse cursor position are found and depending on that a different
        % scroll function is set. Zooming happens in the direction of the
        % current mouse cursor position.
        % For zooming out of the figure with double click, check the
        % buttonFcnDown from the respective image!
        % ---
        % check in which figure the scrolling was done;
        
        % check if it is on the data axis
        curpos = get(h.D_axes,'CurrentPoint');
        if ~(curpos(1,1)>=image_info.xlim(1) && curpos(1,1)<=image_info.xlim(2) && curpos(1,2)>=image_info.ylim(1) && curpos(1,2)<=image_info.ylim(2) && curpos(1,3)==1)
            return
        end
        
        % depending if a new ROI is being made, the behavior of the scroll changes
        if gui_info.roi.new_active
            % Define the search radius for semi-automatic ROI detection.
            % This is done by scrolling
            % ---
            % ROI radius increase/decrease in pixels
            delta_radius = 2;
            % increase or decrease
            if eventdata.VerticalScrollCount < 0
                gui_info.roi.radius = gui_info.roi.radius + delta_radius;
            elseif eventdata.VerticalScrollCount > 0
                gui_info.roi.radius = max([0,gui_info.roi.radius - delta_radius]);
            end
            % display the radius as a circle around pointer
            % axes(h.D_axes)
            % hold on
            % get pointer position in D_axes
            %                     p=get(h.D_axes,'CurrentPoint');
            %                     p = p(1,1:2);
            % plot radius and cell thickness for 0.2 seconds
            tmp_radius = plot(curpos(1,1)+gui_info.roi.radius*cos(0:0.01:2*pi)...
                ,curpos(1,2)+gui_info.roi.radius*sin(0:0.01:2*pi),'--m','LineWidth',2);
            pause(0.3)
            delete(tmp_radius)
        else
            % zoom in or out of the image
            % -------------------------
            % determine zoom factor
            zoom_factor = [0.75,1.25];
            % scroll up: ZOOM IN TO POINTER
            if eventdata.VerticalScrollCount < 0
                % scale
                new_width = zoom_factor(1)*diff(image_info.xlim)/2;
                new_height = zoom_factor(1)*diff(image_info.ylim)/2;
                % scroll down: ZOOM OUT TO POINTER
            elseif eventdata.VerticalScrollCount > 0
                % scale
                new_width = zoom_factor(2)*diff(image_info.xlim)/2;
                new_height = zoom_factor(2)*diff(image_info.ylim)/2;
                % if larger than original, set original
                if new_width>(data_info.width-1)/2 || new_height>(data_info.height-1)/2
                    image_info.xlim=[1 data_info.width];
                    image_info.ylim=[1 data_info.height];
                    return
                end
                % otherwise: DO NOTHING
            else
                return
            end
            % center display around cursor
            new_xlim = [curpos(1,1)-new_width,curpos(1,1)+new_width];
            new_ylim = [curpos(1,2)-new_height,curpos(1,2)+new_height];
            % shift for not going outside image
            if new_xlim(1)<1
                new_xlim = new_xlim + (1-new_xlim(1));
            elseif new_xlim(2)>data_info.width
                new_xlim = new_xlim + (data_info.width-new_xlim(2));
            end
            % shift y for not going outside image
            if new_ylim(1)<1
                new_ylim = new_ylim + (1-new_ylim(1));
            elseif new_ylim(2)>data_info.height
                new_ylim = new_ylim + (data_info.height-new_ylim(2));
            end
            % set and update new axes limits
            image_info.xlim = new_xlim;
            image_info.ylim = new_ylim;
            set(h.D_axes,'xlim',image_info.xlim)
            set(h.D_axes,'ylim',image_info.ylim)
        end
    end

    function gui_click(hObject,eventdata)
        % Function than handles clicks on the axes. In both cases the axes
        % are zoomed out
        switch get(hObject,'Tag')
            case 'D_axes'
                % check if double clicked to zoom out and leave
                if strcmp(get(hfig,'selectiontype'),'open')
                    image_info.xlim = [1 data_info.width];
                    image_info.ylim = [1 data_info.height];
                    
                    set(h.D_axes,'xlim',image_info.xlim)
                    set(h.D_axes,'ylim',image_info.ylim)
                end
        end
    end

    function gui_resize(hObject,eventdata)
        % Change position of elements when the window size changes
        pos_fig = get(hfig,'Position');
        
        % Control panel position
        pos_controlPanel = get(h.P_controls,'Position');
        set(h.P_controls,'Position',...
            [pos_controlPanel(1), pos_fig(4)/2 - pos_controlPanel(4)/2  ,pos_controlPanel(3:4)]);
        
        blank_space = 1*abs(pos_controlPanel(1));
        ini_pos_x = pos_controlPanel(1)+ pos_controlPanel(3);
        
        % image axis position
        set(h.D_axes,'Position',...
            [ini_pos_x + blank_space, blank_space, pos_fig(3)-ini_pos_x-2*blank_space,pos_fig(4)-2*blank_space])
        
    end

    function gui_update(updateTag)
        % Function to update the GUI (i.e. enable/disable functions).
        % 'updateTag' detemines which part of the GUI has to be modified and
        % ----
        switch updateTag
            case {'ROI_exclude','ROI_include'}
                % ADD or REMOVE parts of the ROI
                if gui_info.roi.selected == 0
                    set(h.R_exclude,'Enable','off');
                    set(h.R_include,'Enable','off');
                else
                    set(h.R_exclude,'Enable','on');
                    set(h.R_include,'Enable','on');
                end
            case 'ROI_button'
                % change buttonDownFcn for the image axis
                h_tmp = findobj('Tag','D_axes');
                if get(h.R_new,'Value') == get(h.R_new,'Max')
                    set(h_tmp,'ButtonDownFcn',{@ROI_make_new,'NEW'})
                else
                    set(h_tmp,'ButtonDownFcn',@gui_click)
                end
                % change color
                if get(h.R_new,'Value') == get(h.R_new,'Max')
                    set(h.R_new,'BackgroundColor',[0.9333 0.8667 0.5098]);
                else
                    set(h.R_new,'BackgroundColor',[1 1 1]*0.941);
                end
                % change cursor and scroll function depending on selection type
                if gui_info.roi.new_active ...
                        && (strcmp(strtrim(get(get(h.R_make_type,'SelectedObject'),'String')),'CIRCLE') ...
                        || strcmp(strtrim(get(get(h.R_make_type,'SelectedObject'),'String')),'RING'))
                    setptr(hfig,'crosshair');
                else
                    setptr(hfig,'arrow');
                end
            case 'reseed_button'
                % Make the re-seed button only available if automatic
                % method is selected
                if gui_info.roi.selected==0
                    set(h.R_reseed,'Enable','off')
                    set(h.R_seedAll,'Enable','off')
                    return
                end
                switch get(get(h.R_make_type,'SelectedObject'),'Tag')
                    case {'R_drop','R_peak'}
                        set(h.R_reseed,'Enable','on')
                        set(h.R_seedAll,'Enable','on')
                    otherwise
                        set(h.R_reseed,'Enable','off')
                        set(h.R_seedAll,'Enable','off')
                end
            case 'IMG_slider'   % ========================================
                % change slider properties if FILTH_ROI is chosen
                switch strtrim(gui_info.image.display)
                    case 'Filth ROI'
                        set(h.D_sld_frames,'Min',0)
                        set(h.D_sld_frames,'Max',1)
                        set(h.D_sld_frames,'Value',gui_info.image.filth_th )
                        set(h.D_sld_frames,'SliderStep',[0.001 0.1])
                        set(h.D_sld_frames,'Enable','on')
                    otherwise
                        if data_info.num_stacks>1
                            set(h.D_sld_frames,'Enable','on')
                            set(h.D_sld_frames,'Max',data_info.num_stacks)
                            set(h.D_sld_frames,'Min',1)
                            set(h.D_sld_frames,'SliderStep',[1 0.1*data_info.num_stacks]/data_info.num_stacks)
                            set(h.D_sld_frames,'Value',gui_info.data.frame)
                        else
                            set(h.D_sld_frames,'Enable','off')
                            set(h.D_type_frame,'Enable','off')
                        end
                end
            case 'Export_ready'
                ch = get(h.R_doPlot,'Children');
                % ROI
                if gui_info.control.roi_export_ready
                    set(h.R_process_export,'String','EXPORT ROI')
                    set(h.R_process_export,'BackgroundColor',[0.6039 0.8039 0.1961])
                    set(ch(strcmp(get(ch,'Label'),'ROI contour')),'Enable','on')
                    set(ch(strcmp(get(ch,'Label'),'ROI background')),'Enable','on')
                else
                    set(h.R_process_export,'String','PROCESS ROI')
                    set(h.R_process_export,'BackgroundColor',[1 1 1]*0.941)
                    set(ch(strcmp(get(ch,'Label'),'ROI contour')),'Enable','on')
                    set(ch(strcmp(get(ch,'Label'),'ROI background')),'Enable','off')
                end
            case 'R_shift'
                if get(h.R_shift,'Value') == 1
                    set(h.R_shift,'BackgroundColor',[0.9333 0.8667 0.5098])
                    set(hfig,'KeyPressFcn', @ROI_shift);
                else
                    set(h.R_shift,'BackgroundColor',[1 1 1]*0.941)
                    set(hfig,'KeyPressFcn', '');
                end
        end
    end

% ------------------------------------------------------------------- PLOTS

    function plot_image()
        % Plot the current image (image_info) and the ROIs (roi_info).
        % ---
        % prepare image_info.plotted
        set(h.D_axes,'Visible','off')
        switch strtrim(gui_info.image.display)
            case 'Filth ROI'
                % ------ OPTION 1: PLOT FILTH ROI
                % find first channel that is active and clip to boundaries
                ch = find(gui_info.data.channel & data_info.channel_available,1);
                tmp = image_info.extracted(:,:,ch);
                img = (tmp - min(tmp(:)))/(max(tmp(:)) - min(tmp(:)));
                img(img>1)=1;
                img(img<0)=0;
                % make black and white picture of channel
                image_info.plotted=cat(3,img,img,img);
                % find which values in the original plotted image are above the threshold
                ind = find(image_info.plotted(:,:,ch)>gui_info.image.filth_th);
                % color found values in blue
                if ~isempty(ind)
                    image_info.plotted(ind)=0;
                    image_info.plotted(numel(img)+ind)=0;
                    image_info.plotted(2*numel(img)+ind)=1;
                end
                % add the pixels to the FILTH ROI
                data_info.filth(:) = false;
                data_info.filth(ind) = true;
            otherwise
                % ------ OPTION 2: PLOT IMAGE
                % list of channels to display
                % make empty 3-dim figure and fill channels after scaling
                image_info.plotted = zeros(data_info.height,data_info.width,3);
                ch_2_plot = find(gui_info.data.channel & data_info.channel_available);
                for ch = ch_2_plot
                    image_info.plotted(:,:,ch) = (image_info.extracted(:,:,ch) - gui_info.data.threshold(ch,1))...
                        /(gui_info.data.threshold(ch,2) - gui_info.data.threshold(ch,1));
                end
                % if only one channel to plot, repeat in other channels to obtain
                % black and white figure
                if numel(ch_2_plot)==1
                    % make 2d figure and scale
                    image_info.plotted(:,:,1) = image_info.plotted(:,:,ch_2_plot);
                    image_info.plotted(:,:,2) = image_info.plotted(:,:,ch_2_plot);
                    image_info.plotted(:,:,3) = image_info.plotted(:,:,ch_2_plot);
                elseif strcmp(gui_info.image.display,'Difference')
                    tmp = (image_info.plotted(:,:,ch_2_plot(2))-image_info.plotted(:,:,ch_2_plot(1)))./image_info.plotted(:,:,ch_2_plot(1));
                    % make 2d figure and scale
                    image_info.plotted(:,:,1) = tmp;
                    image_info.plotted(:,:,2) = tmp;
                    image_info.plotted(:,:,3) = tmp;
                end
                % clip image
                image_info.plotted(image_info.plotted<0) = 0;
                image_info.plotted(image_info.plotted>1) = 1;
                
        end
        % start or change the image
        if isempty(gui_info.image.handle)
            % make plot for the first time
            axes(h.D_axes)
            gui_info.image.handle = image([1 data_info.width],[1 data_info.height],image_info.plotted(:,:,[2 1 3]));
            set(gui_info.image.handle,'HitTest','off') % because mouse clicks go to original axes category
            % set axes limits
            set(h.D_axes,'xlim',image_info.xlim)
            set(h.D_axes,'ylim',image_info.ylim)
            % send image to bottom
            chH = get(h.D_axes,'Children');
            set(h.D_axes,'Children',[chH(2:end);chH(1)])
        else
            % change image
            set(gui_info.image.handle,'CData',image_info.plotted(:,:,[2 1 3]))
        end
        set(h.D_axes,'Visible','on')
    end

    function plot_ROI(type,old_roi,new_roi)
        % Plot the current image (image_info) and the ROIs (roi_info).
        if nargin==0
            type = 'draw all';
            old_roi = gui_info.roi.selected;
            new_roi = gui_info.roi.selected;
        elseif nargin==1
            old_roi = gui_info.roi.selected;
            new_roi = gui_info.roi.selected;
        end
        
        axes(h.D_axes)
        hold on
        
        % determine color to plot ROIs
        ch_2_plot = find(gui_info.data.channel & data_info.channel_available);
        if numel(ch_2_plot)==1
            roi_color = 'y';
        elseif strcmp(gui_info.image.display,'Difference')
            roi_color = 'y';
        else
            roi_color = 'w';
        end
        
        set(h.D_axes,'Visible','off')
        switch type
            
           case 'delete all'
                % remove all ROIs to make program run faster when adjusting
                % the image
                 for roi_ii=length(roi_info.handle_contour):-1:1
                    try
                        delete(roi_info.handle_contour{roi_ii})
                    end
                 end
                 % plot only current ROI
                 if ~gui_info.image.roi
                     plot_ROI('draw single',old_roi,new_roi)
                 end
            
            case 'draw single'
                if roi_info.count>0
                    ROI = ROI_get_2D(roi_info.ROI_xy{new_roi});
                    [~,roi_info.handle_contour{new_roi}] = contour(h.D_axes,ROI,[1 -1],'LineStyle','-','Color','r','LineWidth',2);
                end
            case 'draw all'
                % delete and remake all ROIs when e.g. a set is loaded
                plot_ROI('delete all',old_roi,new_roi)
                if roi_info.count>0 && ~strcmp(strtrim(gui_info.image.display),'Filth ROI') && gui_info.image.roi
                    for roi_ii = 1:roi_info.count
                        ROI = ROI_get_2D(roi_info.ROI_xy{roi_ii});
                        if roi_ii == gui_info.roi.selected
                            [~,roi_info.handle_contour{roi_ii}] = contour(h.D_axes,ROI,[1 -1],'LineStyle','-','Color','r','LineWidth',2);
                        else
                            switch roi_info.ROI_type{roi_ii}
                                case {'Neuron','Glial cell'}
                                    [~,roi_info.handle_contour{roi_ii}] = contour(h.D_axes,ROI,[1 -1],'LineStyle','-','Color',roi_color,'LineWidth',2);
                                case 'Neuropile'
                                    [~,roi_info.handle_contour{roi_ii}] = contour(h.D_axes,ROI,[1 -1],'LineStyle','-','Color',[0.5 1 0.5],'LineWidth',2);
                                case 'Blood vessel'
                                    [~,roi_info.handle_contour{roi_ii}] = contour(h.D_axes,ROI,[1 -1],'LineStyle','-','Color',[0 0.5 1],'LineWidth',2);
                            end
                        end
                        % see process of drawing
                        drawnow
                    end
                end
               
            case 'color all'
                % change the color of the ROIs when e.g. number of channels
                % changes
                if roi_info.count>0 && ~strcmp(strtrim(gui_info.image.display),'Filth ROI')
                    if ~gui_info.image.roi
                        plot_ROI('delete all',old_roi,new_roi)
                    else
                        for roi_ii = 1:roi_info.count
                            if (roi_ii ~= gui_info.roi.selected)
                                switch roi_info.ROI_type{roi_ii}
                                    case {'Neuron','Glial cell'}
                                        set(roi_info.handle_contour{roi_ii},'LineColor',roi_color);
                                    case 'Neuropile'
                                    case 'Blood vessel'
                                end
                            end
                        end
                    end
                end
                
            case 'color selection'
                % change the color when a different ROI is selected
                
                if ~gui_info.image.roi
                    plot_ROI('delete all',old_roi,new_roi)
                else
                    % Change color or previous
                    if old_roi ~= new_roi && old_roi > 0
                        switch roi_info.ROI_type{old_roi}
                            case {'Neuron','Glial cell'}
                                set(roi_info.handle_contour{old_roi},'LineColor',roi_color);
                            case 'Neuropile'
                                set(roi_info.handle_contour{old_roi},'LineColor',[0.5 1 0.5]);
                            case 'Blood vessel'
                                set(roi_info.handle_contour{old_roi},'LineColor',[0 0.5 1]);
                        end
                    end
                    % change color of new
                    set(roi_info.handle_contour{new_roi},'LineColor','r');
                end
                
            case 'delete single'
                % delete one ROI and change the color of previous one to
                % mark as selected
                % Delete
                
                if ~gui_info.image.roi
                    plot_ROI('delete all',old_roi,new_roi)
                else
                    delete(roi_info.handle_contour{old_roi})
                    roi_info.handle_contour(old_roi) = [];
                    
                    % mark the new ROI in red
                    if new_roi>0
                        % change color of new
                        set(roi_info.handle_contour{new_roi},'LineColor','r');
                    end
                end
                
            case 'change single'
                % modify a ROI and check colors of old and new one
                
                if ~gui_info.image.roi
                    plot_ROI('delete all',old_roi,new_roi)
                else
                    % Delete
                    if old_roi>0
                        delete(roi_info.handle_contour{old_roi})
                    end
                    % change color of previous
                    if old_roi ~= new_roi && old_roi > 0
                        ROI = ROI_get_2D(roi_info.ROI_xy{old_roi});
                        switch roi_info.ROI_type{old_roi}
                            case {'Neuron','Glial cell'}
                                [~,roi_info.handle_contour{old_roi}] = contour(h.D_axes,ROI,[1 -1],'LineStyle','-','Color',roi_color,'LineWidth',2);
                            case 'Neuropile'
                                [~,roi_info.handle_contour{old_roi}] = contour(h.D_axes,ROI,[1 -1],'LineStyle','-','Color',[0.5 1 0.5],'LineWidth',2);
                            case 'Blood vessel'
                                [~,roi_info.handle_contour{old_roi}] = contour(h.D_axes,ROI,[1 -1],'LineStyle','-','Color',[0 0.5 1],'LineWidth',2);
                        end
                    end
                    % mark the new ROI in red
                    ROI = ROI_get_2D(roi_info.ROI_xy{new_roi});
                    [~,roi_info.handle_contour{new_roi}] = contour(h.D_axes,ROI,[1 -1],'LineStyle','-','Color','r','LineWidth',2);
                    set(roi_info.handle_contour{new_roi},'HitTest','off')
                end
                
            otherwise
                
        end
        
        % display label of the ROI
        % delete previous
        for label_ii = 1:length(roi_info.handle_label)
            delete(roi_info.handle_label{label_ii})
        end
        if gui_info.image.label
            % make text of new
            for roi_ii=1:roi_info.count
                ind_y = roi_info.ROI_xy{roi_ii}(2,:);
                ind_x = roi_info.ROI_xy{roi_ii}(1,:);
                roi_info.handle_label{roi_ii} = text(round(mean(ind_x))-2,round(mean(ind_y))+2,num2str(roi_info.ROI_label(roi_ii)),'FontSize',10,'Color',[0 0.5 1],'FontWeight','bold');
            end
        end
        set(h.D_axes,'Visible','on')
    end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          C O N T R O L L E R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function callback_uimenu(hObject,eventdata)
        % ----DATA
        % clear used properties (and calculate again)
        if h.D_prop_clear==hObject
            data_info.handle.metadata_delete()
            gui_info.data.x_correction = 0;
            % UPDATE DATA AND GUI
            prepare_image()
            plot_image()
            return
        end
        % calculate properties using the same alignment
        if h.D_prop_calculate ==hObject
            data_info.handle.properties_calculate(true)
            % UPDATE DATA AND GUI
            prepare_image()
            plot_image()
            return
        end
        % save data properties in meta_data file
        if h.D_prop_save==hObject
            data_info.handle.metadata_save()
            return
        end
        % align frame by frame
        if h.D_prop_internal==hObject
            ch = find(gui_info.data.channel & data_info.channel_available,1);
            data_info.handle.alignment_calculate(ch)
            % UPDATE DATA AND GUI
            prepare_image()
            plot_image()
            return
        end
        % align to reference
        if h.D_prop_reference == hObject
            ch = find(gui_info.data.channel & data_info.channel_available,1);
            data_info.handle.alignment_reference(ch)
            % UPDATE DATA AND GUI
            prepare_image()
            plot_image()
            return
        end
        % do clipping
        if h.D_clip == hObject
            gui_info.data.do_clip = ~gui_info.data.do_clip;
            data_info.handle.do_clip(gui_info.data.do_clip)
            if gui_info.data.do_clip
                set(h.D_clip,'Checked','on');
            else
                set(h.D_clip,'Checked','off');
            end
            % UPDATE DATA AND GUI
            prepare_image()
            plot_image()
            return
        end
        %-- IMAGE DISPLAY
        % show ROIs
        if h.R_show==hObject
            if strcmp(get(h.R_show,'Checked'),'on')
                set(h.R_show,'Checked','off')
                gui_info.image.roi = false;
                plot_ROI('delete all')
            else
                set(h.R_show,'Checked','on')
                gui_info.image.roi = true;
                plot_ROI('draw all')
            end
            
            return
        end
        % show ROI label
        if h.R_ID==hObject
            if strcmp(get(h.R_ID,'Checked'),'on')
                set(h.R_ID,'Checked','off')
                gui_info.image.label = false;
            else
                set(h.R_ID,'Checked','on')
                gui_info.image.label = true;
            end
            plot_ROI('labels')
            return
        end
        % zoom to ROI when clicking
        if h.R_do_zoom==hObject 
            if strcmp(get(h.R_do_zoom,'Checked'),'on')
                set(h.R_do_zoom,'Checked','off')
            else
                set(h.R_do_zoom,'Checked','on')
            end
            return
        end
        % what image to show
        if get(hObject,'Parent')==h.D_displayType
            % uncheck all sisters
            h_list = get(h.D_displayType,'Children');
            for ii=1:length(h_list)
                set(h_list(ii),'Checked','off')
                % except for current one
                if h_list(ii)==hObject
                    set(h_list(ii),'Checked','on')
                end
            end
            gui_info.image.display = get(hObject,'Label');
            set(h.D_filth,'Value',0)
            % UPDATE DATA AND GUI
            gui_update('IMG_slider')
            prepare_image()
            plot_image()
            plot_ROI('color all')
            return
        end
        % -- ROI checker
        if hObject == h.R_invert
            if strcmp(get(h.R_invert,'Checked'),'off')
                set(h.R_invert,'Checked','on')
            else
                set(h.R_invert,'Checked','off')
            end
            return
        end
        % -- SHOW LOADED DATA INFO
        if h.D_info==hObject
            % display current set info
            fname = cell(1,length(data_info.handle.tiff_info));
            for fnum = 1:length(data_info.handle.tiff_info)
                fname{fnum} = ['Filename(',num2str(fnum),'): ',data_info.handle.tiff_info(fnum).Filename];
            end
            txt = [fname,{...
                ['Number of channels: ',num2str(sum(data_info.channel_available))],...
                ['Stacks: ',num2str(data_info.num_stacks)],...
                ['Pixels x: ',num2str(data_info.width)],...
                ['Pixels y: ',num2str(data_info.height)],...
                ['Minimum: ',num2str(data_info.handle.minValue')], ...
                ['Maximum: ',num2str(data_info.handle.maxValue')]}];
            msgbox(txt,'SET INFO')
            drawnow
        end
    end

    function callback_parameters(hObject,eventdata,parameter_type)
        
        switch parameter_type
            case 'detection_parameters'
                
                % close previous window
                h_params = findobj('Tag','h_params');
                close(h_params)
                h_params = figure('Tag','h_params','Name','Detection parameters','Position',[1165 184 500 100]);
                % make table (se below for callback R_act_params)
                t = uitable(h_params,'Tag','set_detect_params','Data',100*[gui_info.roi.threshold;gui_info.roi.percent_in;gui_info.roi.percent_out],'ColumnName',{''},'RowName',{'% signal drop','% ring inside','% ring outside'},...
                    'ColumnWidth',{100},'ColumnEditable',true,'CellEditCallback',{@callback_parameters,'set_detect_params'});
                ext = get(t,'Extent');pos =get(t,'Position');
                set(t,'Position',[pos(1) pos(2) ext(3) ext(4)]);
                
            case 'set_detect_params'
                switch eventdata.Indices(1)
                    case 1
                        gui_info.roi.threshold  = eventdata.NewData/100;
                    case 2
                        gui_info.roi.percent_in = eventdata.NewData/100;
                    case 3
                        gui_info.roi.percent_out = eventdata.NewData/100;
                end
        end
    end

    function callback_data(hObject,eventdata)
        % Callback for the Data panel. The used uicontrol tag (i.e. element
        % in structre 'h') determines the function of the callback. The
        % function deals mostly with tiff data importing and display.
        % ---
        gui_tag = get(hObject,'Tag');
        switch gui_tag
            case 'D_load'
                % load data
                load_data('change')
                plot_ROI('color all')
            case 'D_sld_frames'
                %gui_update('IMG_slider')
                switch strtrim(gui_info.image.display)
                    case 'Filth ROI'
                        gui_info.image.filth_th = get(hObject,'Value');
                    otherwise
                        gui_info.data.frame = round(get(hObject,'Value'));
                        gui_info.image.display = 'Frame';
                        % check the box with Frame display
                        callback_uimenu(h.D_type_frame,'change checkbox')
                end
                % UPDATE DATA AND GUI
                prepare_image()
                plot_image()
            case {'D_green_check','D_red_check','D_blue_check'}
                ch = find(strcmp({'D_green_check','D_red_check','D_blue_check'},gui_tag));
                gui_info.data.channel(ch)=get(h.(gui_tag),'Value')==1;
                % UPDATE DATA AND GUI
                prepare_image()
                plot_image()
                plot_ROI('color all')
            case {'D_green_txt_low','D_green_sld_low','D_green_txt_high','D_green_sld_high',...
                    'D_red_txt_low','D_red_sld_low','D_red_txt_high','D_red_sld_high',...
                    'D_blue_txt_low','D_blue_sld_low','D_blue_txt_high','D_blue_sld_high'}
                % get which channel and uicontrol is called
                tmp  = strfind(gui_tag,'_');
                channel = gui_tag(tmp(1)+1:tmp(2)-1);
                type = gui_tag(tmp(2)+1:tmp(3)-1);
                limit = gui_tag(tmp(3)+1:end);
                % update string - value pair in current object
                if strcmp(type,'txt')
                    set(hObject,'Value',str2double(get(hObject,'String')))
                else
                    set(hObject,'String',num2str(get(hObject,'Value')))
                end
                % get new value and bound
                if strcmp(limit,'low')
                    min_limit = 0;
                    max_limit = get(h.(['D_',channel,'_sld_high']),'Value');
                else
                    min_limit = get(h.(['D_',channel,'_sld_low']),'Value');
                    max_limit = 100;
                end
                newValue = min(max(min_limit,get(hObject,'Value')),max_limit);
                % udate gui
                set(hObject,'Value',newValue)
                set(hObject,'String',num2str(newValue))
                if strcmp(type,'sld')
                    co_tag = ['D_',channel,'_txt_',limit];
                else
                    co_tag = ['D_',channel,'_sld_',limit];
                end
                set(h.(co_tag),'Value',newValue)
                set(h.(co_tag),'String',num2str(newValue))
                % set plot characteristics
                gui_info.data.threshold(...
                    find(strcmp(channel,{'green','red','blue'})),...
                    find(strcmp(limit,{'low','high'}))) = newValue/100;
                % UPDATE DATA AND GUI
                plot_image()
            case {'D_xCor_sld','D_xCor_txt'}
                % slider to determine x_correction
                % get new value
                tmp  = strfind(gui_tag,'_');
                type = gui_tag(tmp(2)+1:end);
                if strcmp(type,'txt')
                    new_value = round(str2double(get(hObject,'String')));
                else
                    new_value = round(get(hObject,'Value'));
                end
                gui_info.data.x_correction = min([max([-49 new_value]) 50]);
                % set in tiff handle
                data_info.handle.set_x_correction(gui_info.data.x_correction)
                % UPDATE DATA AND GUI
                set(h.D_xCor_sld,{'Value'},{gui_info.data.x_correction})
                set(h.D_xCor_txt,{'String'},{num2str(gui_info.data.x_correction)})
                prepare_image()
                plot_image()
            case 'D_frame_rate'
                % change frame rate
                value = str2double(get(hObject,'String'));
                value = max(value,0);
                set(hObject,'String',num2str(value))
                gui_info.data.frame_rate = value;
            case 'D_filth'
                % select a threshold for filth
                if get(h.D_filth,'Value')==1
                    % check off all display types
                    h_list = get(h.D_displayType,'Children');
                    for ii=1:length(h_list)
                        set(h_list(ii),'Checked','off')
                    end
                    gui_info.image.display = 'Filth ROI';
                    % delete all ROI contours
                    plot_ROI('delete all')
                else
                    % check on average
                    set(h.D_type_average,'Checked','on')
                    gui_info.image.display = 'Average';
                    % plot all ROI contours
                    plot_ROI('draw all')
                end
                % UPDATE DATA AND GUI
                gui_update('IMG_slider')
                prepare_image()
                plot_image()
                
            otherwise
                warning('Data callback type not found!')
        end
    end

    function callback_ROI(hObject,eventdata)
        % Callback that imports/removes ROIs and determines how to display
        % them. The used uicontrol tag (i.e. element in structre 'h')
        % determines the function of the callback.
        % ---
        gui_tag = get(hObject,'Tag');
        switch gui_tag
            case 'R_import'
                ROI_import_export('import');
                plot_ROI('draw all')
                
            case 'R_process_export'
                % If ROI processed, save, otherwise, process
                if gui_info.control.roi_export_ready
                    ROI_import_export('export')
                else
                    % Split ROIs that are shared between cells and remove filth
                    ROI_process()
                    % UPDATE DATA AND GUI
                    plot_ROI('draw all')
                end
                
            case 'R_remove'
                % shorten list and choose the last roi for display
                if gui_info.roi.selected>0
                    % decrease number of ROIs if selected ROI is not last
                    % on the list (i.e just made)
                    if gui_info.roi.selected == roi_info.count
                        gui_info.roi.label_count = gui_info.roi.label_count - 1;
                    end
                    % update ROI display
                    roi_info.ROI_xy(gui_info.roi.selected)=[];
                    roi_info.ROI_label(gui_info.roi.selected)=[];
                    roi_info.ROI_type(gui_info.roi.selected)=[];
                    gui_info.roi.list(gui_info.roi.selected) = [];
                    roi_info.count = max(0,roi_info.count - 1);
                    plot_ROI('delete single',gui_info.roi.selected,roi_info.count)
                    gui_info.roi.selected = min([gui_info.roi.selected roi_info.count]);
                    set(h.('R_list'),'String',gui_info.roi.list)
                    set(h.('R_list'),'Value',gui_info.roi.selected)
                    % UPDATE DATA AND ROI
                    gui_info.control.roi_export_ready = false;
                    gui_update('Export_ready')
                    gui_update('ROI_exclude')
                    
                end
                
            case 'R_overwrite'
                % overwrite the current ROI with next new one
                if gui_info.roi.selected>0
                    % deactivate XY shift
                    if get(h.R_overwrite,'Value')
                        set(h.R_shift,'Value',0)
                    end
                    gui_update('R_shift')
                else
                    set(h.R_overwrite,'Value',0);
                end
                gui_info.control.roi_export_ready = false;
                gui_update('Export_ready')
                
            case 'R_sort'
                % Put new labels to the ROIs
                % close previous window
                h_sort = findobj('Tag','hsort');
                close(h_sort)
                h_sort = figure('Tag','hsort','Name','Label sorter','Position',[1165 184 442 751]);
                % make table (se below for callback R_sort_table)
                t = uitable(h_sort,'Tag','R_sort_table','Data',roi_info.ROI_label,'ColumnName',{'ROI Label'},...
                    'ColumnWidth',{100},'ColumnEditable',true,'CellEditCallback',@callback_ROI);
                ext = get(t,'Extent');pos =get(t,'Position');
                set(t,'Position',[pos(1) pos(2) ext(3) ext(4)]);
                
            case 'R_sort_table'
                % Check that labels are not repeated when changing them
                if ismember(eventdata.NewData,roi_info.ROI_label)
                    disp('Label exists')
                    set(hObject,'Data',roi_info.ROI_label);
                else
                    roi_info.ROI_label = get(hObject,'Data');
                    gui_info.roi.label_count = max(roi_info.ROI_label);
                    % change name in list
                    ind = eventdata.Indices(1);
                    switch roi_info.ROI_type{ind}
                        case {'Neuron','Glial cell'}
                            gui_info.roi.list{ind} = ['<HTML><FONT color="black">ROI ',num2str(roi_info.ROI_label(ind))];
                        case 'Neuropile'
                            gui_info.roi.list{ind} = ['<HTML><FONT color="green">ROI ',num2str(roi_info.ROI_label(ind))];
                        case 'Blood vessel'
                            gui_info.roi.list{ind} = ['<HTML><FONT color="blue">ROI ',num2str(roi_info.ROI_label(ind))];
                    end
                    set(h.('R_list'),'String',gui_info.roi.list)
                    % mark ROIs as NOT ready to export
                    gui_info.control.roi_export_ready = false;
                    gui_update('Export_ready')
                    % change on plot
                    plot_ROI('labels')
                    % go back to selector
                    figure(findobj('Tag','hsort'))
                end
            case 'R_list'
                % Get selected ROI
                previous = gui_info.roi.selected;
                gui_info.roi.selected = get(hObject,'Value');
                % UPDATE DATA AND GUI
                switch roi_info.ROI_type{gui_info.roi.selected}
                    case 'Neuron'
                        set(h.R_ROI_type,'SelectedObject',h.R_type_neuron)
                    case 'Glial cell'
                        set(h.R_ROI_type,'SelectedObject',h.R_type_glia)
                    case 'Neuropile'
                        set(h.R_ROI_type,'SelectedObject',h.R_type_np)
                    case 'Blood vessel'
                        set(h.R_ROI_type,'SelectedObject',h.R_type_bv)
                end
                plot_ROI('color selection',previous,gui_info.roi.selected)
                % if zoom to axis check in menu, zoom-in
                % zoom in to the current ROI
                if strcmp(get(h.R_do_zoom,'Checked'),'on')
                    % get center of ROI and radius
                    ind_y = roi_info.ROI_xy{gui_info.roi.selected}(2,:);
                    ind_x = roi_info.ROI_xy{gui_info.roi.selected}(1,:);
                    center_x = round(mean(ind_x));
                    center_y = round(mean(ind_y));
                    radius = mean(sqrt((ind_y-mean(ind_y)).^2+(ind_x-mean(ind_x)).^2));
                    % get new width and height
                    new_width = 4*radius;
                    new_height = 4*radius;
                    new_width = new_width * data_info.width/data_info.height;
                    % center display around cursor
                    new_xlim = [center_x-new_width,center_x+new_width];
                    new_ylim = [center_y-new_height,center_y+new_height];
                    % shift for not going outside image
                    if new_xlim(1)<1
                        new_xlim = new_xlim + (1-new_xlim(1));
                    elseif new_xlim(2)>data_info.width
                        new_xlim = new_xlim + (data_info.width-new_xlim(2));
                    end
                    % shift y for not going outside image
                    if new_ylim(1)<1
                        new_ylim = new_ylim + (1-new_ylim(1));
                    elseif new_ylim(2)>data_info.height
                        new_ylim = new_ylim + (data_info.height-new_ylim(2));
                    end
                    % set and update new axes limits
                    image_info.xlim = new_xlim;
                    image_info.ylim = new_ylim;
                    set(h.D_axes,'xlim',image_info.xlim)
                    set(h.D_axes,'ylim',image_info.ylim)
                end
            case 'R_new'
                % deactivate XY shift
                if get(h.R_new,'Value')
                    set(h.R_shift,'Value',0)
                    gui_info.roi.new_active = true;
                else
                    gui_info.roi.new_active = false;
                end
                gui_update('R_shift')
                % change GUI propperties and change type to cell
                gui_update('ROI_button')
                if get(h.R_new,'Value') || gui_info.roi.selected == 0
                    set(h.R_ROI_type,'SelectedObject',h.R_type_neuron)
                else
                    switch roi_info.ROI_type{gui_info.roi.selected}
                        case 'Neuron'
                            set(h.R_ROI_type,'SelectedObject',h.R_type_neuron)
                        case 'Glial cell'
                            set(h.R_ROI_type,'SelectedObject',h.R_type_glia)
                        case 'Neuropile'
                            set(h.R_ROI_type,'SelectedObject',h.R_type_np)
                        case 'Blood vessel'
                            set(h.R_ROI_type,'SelectedObject',h.R_type_bv)
                    end
                end
            case 'R_make_type'
                % change GUI if R_new is on
                gui_update('reseed_button')
            case 'R_plusLabel'
                % Increase label of selected ROI by the largest +1
                if gui_info.roi.selected>0
                    % add 1 to GUI roi list
                    gui_info.roi.label_count = gui_info.roi.label_count +1;
                    roi_info.ROI_label(gui_info.roi.selected) = gui_info.roi.label_count;
                    switch roi_info.ROI_type{gui_info.roi.selected}
                        case {'Neuron','Glial cell'}
                            gui_info.roi.list{gui_info.roi.selected} = ['<HTML><FONT color="black">ROI ',num2str(gui_info.roi.label_count)];
                        case 'Neuropile'
                            gui_info.roi.list{gui_info.roi.selected} = ['<HTML><FONT color="green">ROI ',num2str(gui_info.roi.label_count)];
                        case 'Blood vessel'
                            gui_info.roi.list{gui_info.roi.selected} = ['<HTML><FONT color="blue">ROI ',num2str(gui_info.roi.label_count)];
                    end
                    % update GUI
                    set(h.R_list,'String',gui_info.roi.list)
                end
            case 'R_shift'
                % shift the position of ALL ROIs with the keyboard arrows
                if gui_info.roi.selected>0
                    if get(h.R_shift,'Value')
                        set(h.R_new,'Value',0)
                        set(h.R_overwrite,'Value',0)
                    end
                else
                    set(h.R_shift,'Value',0);
                end
                gui_info.control.roi_export_ready = false;
                gui_update('Export_ready')
                gui_update('ROI_button')
                gui_update('R_shift')
            case 'R_align'
                % find parameters to align the ROIs
                if gui_info.roi.selected>0
                    % do the alignment
                    ROI_align()
                    % plot new
                    gui_info.control.roi_export_ready = false;
                    gui_update('Export_ready')
                    plot_ROI('draw all')
                end
            case 'R_ROI_type'
                % change the type of the ROI (cell/neuropile/blood vessel)
                if gui_info.roi.selected>0
                    roi_info.ROI_type{gui_info.roi.selected} = get(get(h.R_ROI_type,'SelectedObject'),'String');
                    switch roi_info.ROI_type{gui_info.roi.selected}
                        case {'Neuron','Glial cell'}
                            gui_info.roi.list{gui_info.roi.selected} = ['<HTML><FONT color="black">ROI ',num2str(roi_info.ROI_label(gui_info.roi.selected))];
                        case 'Neuropile'
                            gui_info.roi.list{gui_info.roi.selected} = ['<HTML><FONT color="green">ROI ',num2str(roi_info.ROI_label(gui_info.roi.selected))];
                        case 'Blood vessel'
                            gui_info.roi.list{gui_info.roi.selected} = ['<HTML><FONT color="blue">ROI ',num2str(roi_info.ROI_label(gui_info.roi.selected))];
                    end
                    % update list
                    set(h.R_list,'String',gui_info.roi.list);
                    % selected color is red, does not need to be changed.
                    % as the type of the neurons changed, the ROIs have to be
                    % processed again and the the signal has to be re-calculated
                    gui_info.control.roi_export_ready = false;
                    gui_update('Export_ready')
                end
            case 'R_clear'    
                % clear all ROI data
                
                % remove from plot
                plot_ROI('delete all')
                
                % empty data
                gui_info.roi.label_count = 0;
                gui_info.roi.list = {};
                gui_info.roi.selected = 0;
                
                roi_info = [];
                roi_info.count = 0;
                roi_info.ROI_xy = {};
                roi_info.ROI_label = [];
                roi_info.ROI_type = {};
                roi_info.handle_contour = {};
                roi_info.handle_label = {};                
                % UPDATE DATA AND ROI
                gui_info.control.roi_export_ready = false;
                gui_update('Export_ready')
                gui_update('ROI_exclude')
                
                set(h.('R_list'),'String',gui_info.roi.list)
                    set(h.('R_list'),'Value',roi_info.count)
            otherwise
                warning('ROI callback type not found!')
        end
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              M O D E L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------------------------------------------------- DATA FUNCTIONS

    function load_data(type,file_structure)
        % open a given tif file, or open with GUI from initial path
        % type tells how it should be opened:
        % - 'start' : when the program starts
        % - 'change' : change the set opening a GUI to check
        % - 'external' : when a FILE_STRUCTURE is given to bypass the GUI
        
        try
            switch type
                case {'start','change'}
                    tmp_handle = tiff_reader(gui_info.control.path);
                case 'external'
                    tmp_handle = tiff_reader(file_structure);
            end
        catch ierr
            warning('Error while loading tiff files...')
            return
        end
        
        data_info = [];
        data_info.handle = tmp_handle;
        % set new path
        tmp=strfind(data_info.handle.tiff_info(1).Filename,filesep);
        gui_info.control.path = data_info.handle.tiff_info(1).Filename(1:tmp(end));
        % get data properties
        data_info.channel_available = [false,false,false];
        data_info.channel_available(1:data_info.handle.num_channels) = true;
        data_info.width = data_info.handle.width;
        data_info.height = data_info.handle.height;
        data_info.num_stacks = data_info.handle.num_stacks;
        gui_info.data.x_correction = sum(data_info.handle.x_correction);
        data_info.handle.do_clip(gui_info.data.do_clip)
        % make a new FILTH_ROI
        data_info.filth = false(data_info.height,data_info.width);
        
        % if only starting, interrupt the rest
        if ~strcmp(type,'change')
            return
        end
        
        % === check that ROI x-y size matches set
        max_pix = max(cat(2,roi_info.ROI_xy{:}),[],2);
        if roi_info.count>0 && (max_pix(1)>data_info.width || max_pix(2)>data_info.height)
            disp('WARNING: Frame dimensions do not match current ROIs')
            disp('         Pixels outside the image range will be removed when processing...')
        end
        
        % === UPDATE GUI
        % check number of channels and activate difference image
        if sum(data_info.channel_available)>1
            set(h.D_type_difference,'Enable','on')
        else
            set(h.D_type_difference,'Enable','off')
        end
        if ~isempty(data_info.handle.img_correlation)
            set(h.D_type_correlation,'Enable','on')
        else
            set(h.D_type_correlation,'Enable','off')
        end
        if ~isempty(data_info.handle.img_averagedMaximum)
            set(h.D_type_chunkyMax,'Enable','on')
        else
            set(h.D_type_chunkyMax,'Enable','off')
        end
        % change x_corection
        set(h.D_xCor_sld,{'Value'},{gui_info.data.x_correction})
        set(h.D_xCor_txt,{'String'},{num2str(gui_info.data.x_correction)})
        % make new channels available
        gui_info.data.channel = data_info.channel_available;
        textN = {'green','red','blue'};
        for chN=1:3
            if gui_info.data.channel(chN)
                enable = 'on';
            else
                enable = 'off';
            end
            set(h.(['D_',textN{chN},'_check']),{'Value','Enable'},{gui_info.data.channel(chN),enable})
            set(h.(['D_',textN{chN},'_sld_low']),{'Enable'},{enable})
            set(h.(['D_',textN{chN},'_sld_high']),{'Enable'},{enable})
            set(h.(['D_',textN{chN},'_txt_low']),{'Enable'},{enable})
            set(h.(['D_',textN{chN},'_txt_high']),{'Enable'},{enable})
        end
        %             set(h.D_sld_frames,{'Min','Max','Value','SliderStep'},...
        %                 {1,100,1,[1 0.1*data_info.num_stacks]/data_info.num_stacks})
        
        % change it to mark as not exportable
        gui_update('Export_ready')
        % update image
        prepare_image()
        plot_image()
        
    end

    function load_exported(hObject,eventdata)
        % Use an exported ROIs mat file to load both data and ROIs
        % get file
        [FileName,PathName] = uigetfile([gui_info.control.path,'*.mat'],'MultiSelect','off','Select EXPORTED file to open');
        file_to_open = [PathName,FileName];
        opened = load(file_to_open,'ROI','data_class');
        % check that everything that is needed is in the file
        if ~isfield(opened,'ROI') || ~isfield(opened,'data_class')
            % if _extracted file, check for source_file
            opened = load(file_to_open,'source_file');
            if ~isfield(opened,'source_file')
                error('File does not contain the necessary information to load data + ROI')
            end
            file_to_open = opened.source_file;
            opened = load(file_to_open,'ROI','data_class');
            if ~isfield(opened,'ROI') || ~isfield(opened,'data_class')
                error('File does not contain the necessary information to load data + ROI')
            end
        end
        % load full file
        opened = load(file_to_open,'ROI','data_class');
        % --------------------------------------------------
        % load the data
        load_data('external',opened.data_class.files)
        % change x correction
        if isfield(opened.data_class,'x_correction')
            gui_info.data.x_correction = opened.data_class.x_correction;
            data_info.handle.set_x_correction(gui_info.data.x_correction)
        end
        % change clip value
        if isfield(opened.data_class,'clip_value')
            if any(opened.data_class.clip_value>0)
                gui_info.data.do_clip = true;
                data_info.handle.do_clip(true,opened.data_class.clip_value)
            else
                gui_info.data.do_clip = false;
            end
        end
        % change frame rate
        if isfield(opened.data_class,'frame_rate')
            gui_info.data.frame_rate = opened.data_class.frame_rate;
        end
        % change filth ROI
        if isfield(opened,'ROI_filth')
            data_info.filth = opened.ROI_filth;
        end
        % --------------------------------------------------
        % Clear previous figure
        if exist('h_fig','var')
            delete hfig
        end
        clear h
        % draw the interface with callbacks and get handles
        [hfig,h]=gui_create();
        % --------------------------------------------------
        % load the ROIs
        ROI_import_export('import',file_to_open)
        % --------------------------------------------------
        % ---> plot the first image
        gui_info.image.handle = [];        
        prepare_image()
        plot_image()
        roi_info.handle_contour = {};
        roi_info.handle_label = {};
        plot_ROI()
    end

    function prepare_image()
        % Prepare the image to be displayed in the Data axis. The function
        % determines the type of image to be read (average,maximum,variance,frame, fith)
        % extracts it from data_info and thresholds it. The result is stored
        % in image_info.
        % ---
        % make empty 3D figure and fill by channels after fitting range to [0 1]
        image_info.extracted = zeros(data_info.height,data_info.width,3);
        
        if strcmp(strtrim(gui_info.image.display),'Frame')
            % get the general maximum and minimum pixel value for all sets
            total_min = data_info.handle.minValue;
            total_max = data_info.handle.maxValue;
            % get the frame from the selected file
            for ch = find(data_info.channel_available)
                tmp = data_info.handle.read_stack(gui_info.data.frame,ch);
                image_info.extracted(:,:,ch) = (tmp - total_min(ch))/(total_max(ch) - total_min(ch));
            end
            
        else
            % determine which image to show
            switch strtrim(gui_info.image.display)
                case 'Maximum'
                    tmp = data_info.handle.maximum;
                case {'Average','Difference'}
                    tmp = data_info.handle.mean;
                case 'Variance'
                    tmp = sqrt(data_info.handle.variance)./data_info.handle.mean;
                case 'Filth ROI'
                    tmp = data_info.handle.mean./sqrt(data_info.handle.variance);
                case 'Correlation'
                    tmp = data_info.handle.correlation;
                    % add to all channels
                    for ch=2:sum(data_info.channel_available)
                        tmp = cat(3,tmp,zeros(size(tmp)));
                    end
                 case 'Chunky Max'
                    tmp = data_info.handle.averagedMaximum;
                    % add to all channels
                    for ch=2:sum(data_info.channel_available)
                        tmp = cat(3,tmp,zeros(size(tmp)));
                    end
            end
            
            % normalize each channel and assign to final image
            for ch=1:sum(data_info.channel_available)
                tmp2 = tmp(:,:,ch);
                tmp2 = (tmp2-min(tmp2(:)))/(max(tmp2(:))-min(tmp2(:)));
                image_info.extracted(:,:,ch) = tmp2;
            end
        end
    end

% ----------------------------------------------------------- ROI FUNCTIONS

    function ROI_import_export(type,file_to_open)
        % Function to convert an imported ROI (2D image with numbers
        % marking ROI number) to the format used by the GUI to handle it
        % (3D boolean array) and back.
        switch type
            case 'import'
                % get file name and try opening it to extract ROI
                if nargin==1
                    [FileName,PathName] = uigetfile([gui_info.control.path,'*.mat'],'MultiSelect','off','Select file with saved ROIs');
                    file_to_open = [PathName,FileName];
                end
                
                try
                    opened = load(file_to_open);
                    
                    % use old exported ROIs
                    if isfield(opened,'roi_info')
                        opened = opened.roi_info;
                    end
                    
                    % get original size of ROI
                    [size_y,size_x,~] = size(opened.ROI);
                    
                    % convert 2D to 3D array if necessary
                    if ndims(opened.ROI)==2
                        tmp_ROI.count = max(opened.ROI(:));
                        tmp_ROI.ROI_xy = cell(tmp_ROI.count,1);
                        for roi_ii = 1:tmp_ROI.count
                            [pix_y,pix_x] = find(opened.ROI == roi_ii);
                            tmp_ROI.ROI_xy{roi_ii} = [pix_x(:)';pix_y(:)'];
                        end
                    else
                        tmp_ROI.count = size(opened.ROI,3);
                        tmp_ROI.ROI_xy = cell(tmp_ROI.count,1);
                        for roi_ii = 1:tmp_ROI.count
                            [pix_y,pix_x] = find(opened.ROI(:,:,roi_ii));
                            tmp_ROI.ROI_xy{roi_ii} = [pix_x(:)';pix_y(:)'];
                        end
                    end
                    
                    % load ROI labels or set to ROI number
                    if isfield(opened,'ROI_label')
                        tmp_ROI.ROI_label = opened.ROI_label;
                    else
                        tmp_ROI.ROI_label = 1:tmp_ROI.count;
                    end
                    
                    % load ROI types or set to Neuron
                    if isfield(opened,'ROI_type')
                        tmp_ROI.ROI_type = opened.ROI_type;
                    else
                        tmp_ROI.ROI_type = repmat({'Neuron'},[tmp_ROI.count,1]);
                    end
                    
                    % add handles to make figure
                    tmp_ROI.handle_contour = cell(tmp_ROI.count,1);
                    tmp_ROI.handle_label = cell(tmp_ROI.count,1);
                catch ierr
                    disp(ierr)
                    error('Could not load file with ROIs!')
                end
                
                
                % check that imported ROI x-y size matches current data
                if abs(data_info.width - size_x)>0 || abs(data_info.height - size_y)>0
                    disp('WARNING: ROI dimensions do not match current set!')
                    disp('         Shifting ROIs, pixels outside the frame will be removed when processing ...')
                    
                    dx = round((data_info.width - size_x)/2);
                    dy = round((data_info.height - size_y)/2);
                    
                    % shift ROIs
                    for roi_ii = 1:tmp_ROI.count
                        tmp_ROI.ROI_xy{roi_ii}(1,:) = tmp_ROI.ROI_xy{roi_ii}(1,:) + dx;
                        tmp_ROI.ROI_xy{roi_ii}(2,:) = tmp_ROI.ROI_xy{roi_ii}(2,:) + dy;
                    end
                end
                
                % delete previous and add new
                for roi_ii = 1:roi_info.count
                    delete(roi_info.handle_contour{roi_ii})
                end
                roi_info = tmp_ROI;
                
                % ======> modify GUI
                gui_info.roi.label_count = max(roi_info.ROI_label);
                gui_info.roi.selected = roi_info.count;
                gui_info.roi.list = cell(1,roi_info.count);
                for ind = 1:roi_info.count
                    switch roi_info.ROI_type{ind}
                        case {'Neuron','Glial cell'}
                            gui_info.roi.list{ind} = ['<HTML><FONT color="black">ROI ',num2str(roi_info.ROI_label(ind))];
                        case 'Neuropile'
                            gui_info.roi.list{ind} = ['<HTML><FONT color="green">ROI ',num2str(roi_info.ROI_label(ind))];
                        case 'Blood vessel'
                            gui_info.roi.list{ind} = ['<HTML><FONT color="blue">ROI ',num2str(roi_info.ROI_label(ind))];
                    end
                end
                % update GUI
                set(h.('R_list'),'String',gui_info.roi.list)
                set(h.('R_list'),'Value',roi_info.count)
                gui_info.control.roi_export_ready = false;
                gui_update('ROI_exclude')
                gui_update('Export_ready')
                
            case 'export'
                % check that the Frame Rate was set
                if isempty(gui_info.data.frame_rate)
                   msgbox('Before exporting define Frame Rate [Hz]','ERROR')
                   drawnow
                   return
                end
                % get file destination
                [FileName,PathName] = uiputfile([gui_info.control.path,'saved_ROIs.mat'],'Save defined ROIs');
                if ~isnumeric(FileName)
                    % Make 2D ROI and extract labels and types
                    ROI = zeros(data_info.height,data_info.width);
                    for roi_ii = 1:roi_info.count
                        ROI = ROI + roi_ii*ROI_get_2D(roi_info.ROI_xy{roi_ii});
                    end
                    ROI_label = roi_info.ROI_label;
                    ROI_type = roi_info.ROI_type;
                    ROI_filth = data_info.filth;
                    % extract variables with which the data can be reloaded
                    data_class = [];
                    files = cell(size(data_info.handle.tiff_info));
                    for tii=1:size(data_info.handle.tiff_info,1)
                       for tjj=1:size(data_info.handle.tiff_info,2)
                           files{tii,tjj} = data_info.handle.tiff_info(tii,tjj).Filename;
                       end
                    end
                    data_class.files = files;
                    data_class.x_correction = data_info.handle.x_correction;
                    data_class.clip_value = data_info.handle.clip_value;
                    data_class.frame_rate = gui_info.data.frame_rate;
                    % add a snapshot of the image
                    data_img = image_info.plotted;
                    % save
                    time_stamp = now;
                    save([PathName,FileName],'ROI','ROI_label','ROI_type','ROI_filth','data_class','data_img','time_stamp');
                end
        end
    end

    function ROI_make_new(hObject,eventdata,caller_function)
        % Function used to define new ROI. For semi-automatic detection the
        % code in 'ROI_detect.m is used'. In this case the scroll function
        % changes to display the search radius of the cell boundary (see
        % ROI_scroll_radius function). The signal drop/increase that
        % determines the boundary is set in function
        % 'ROI_define_parameters'.
        % Accessed via NEW ROI button or left click on axis
        % ---
        
        roi_selection_type = get(get(h.R_make_type,'SelectedObject'),'Tag');
        % set ROI type to 'Neuron' for new region of interest
        if ~(strcmp(caller_function,'RESEED') || get(h.R_overwrite,'Value') || strcmp(roi_selection_type,'R_exclude') || strcmp(roi_selection_type,'R_include'))         
            set(h.R_ROI_type,'SelectedObject',h.R_type_neuron)
        end
        % make ROI
        switch roi_selection_type
            case 'R_draw'
                % draw the ROI and double click inside to accept
                ROI_handle = imfreehand(h.D_axes,'Closed',true);
                ROI_coordinates=wait(ROI_handle);
                delete(ROI_handle);
                new_ROI = poly2mask(ROI_coordinates(:,1), ROI_coordinates(:,2),data_info.height,data_info.width);
                new_ROI = imfill(new_ROI,'holes');
                if get(h.R_do_ring,'Value')==1
                    % m
                    radius = mean(sqrt(sum(bsxfun(@minus,ROI_coordinates,mean(ROI_coordinates,1)).^2,2)));
                    % dilate
                    out_ROI = imdilate(new_ROI,strel('disk',ceil(gui_info.roi.percent_out*radius)));
                    % erode
                    in_ROI = imerode(new_ROI,strel('disk',ceil(gui_info.roi.percent_in*radius)));
                    % join
                    new_ROI = out_ROI;
                    new_ROI(in_ROI) = false;
                end
            case 'R_include'
                % start with selected ROI
                new_ROI = ROI_get_2D(roi_info.ROI_xy{gui_info.roi.selected});
                % draw region to exclude from ROI
                ROI_handle = imfreehand(h.D_axes,'Closed',true);
                ROI_coordinates=wait(ROI_handle);
                delete(ROI_handle);
                include_ROI = poly2mask(ROI_coordinates(:,1), ROI_coordinates(:,2),data_info.height,data_info.width);
                % combine and check that size matches to add
                new_ROI(include_ROI) = true;
            case 'R_exclude'
                % start with selected ROI
                new_ROI = ROI_get_2D(roi_info.ROI_xy{gui_info.roi.selected});
                % draw region to exclude from ROI
                ROI_handle = imfreehand(h.D_axes,'Closed',true);
                ROI_coordinates=wait(ROI_handle);
                delete(ROI_handle);
                exclude_ROI = poly2mask(ROI_coordinates(:,1), ROI_coordinates(:,2),data_info.height,data_info.width);
                % combine and check that size matches to add
                new_ROI(exclude_ROI) = false;
            case {'R_drop','R_peak'}
                % semi automatic ROI detection happens on the first ACTIVE
                % channel.
                % get current point
                switch caller_function
                    case 'NEW'
                        % get cursor location to detect ROI
                        tmp=get(h.D_axes,'CurrentPoint');
                        center_x=tmp(1,1);
                        center_y=tmp(1,2);
                    case 'RESEED'
                        % use center of mass of current ROI to detect again
                        center_x = round(mean(roi_info.ROI_xy{gui_info.roi.selected}(1,:)));
                        center_y = round(mean(roi_info.ROI_xy{gui_info.roi.selected}(2,:)));
                end
                % find first channel that is active, get image and invert if necessary
                img = image_info.plotted(:,:,find(gui_info.data.channel & data_info.channel_available,1));
                if strcmp(get(h.R_invert,'Checked'),'on')==1
                    img = 1 - img;
                end
                % get ROI boundary
                ROI_coordinates = segmentation_ROI_detect(roi_selection_type,img,[center_x,center_y],gui_info.roi.radius,gui_info.roi.threshold);
                new_ROI = poly2mask(ROI_coordinates(:,1), ROI_coordinates(:,2),data_info.height,data_info.width);
                new_ROI = imfill(new_ROI,'holes');
                if get(h.R_do_ring,'Value')==1 || strcmp(roi_selection_type,'R_peak')
                    % m
                    radius = mean(sqrt(sum(bsxfun(@minus,ROI_coordinates,mean(ROI_coordinates,1)).^2,2)));
                    % dilate
                    out_ROI = imdilate(new_ROI,strel('disk',ceil(gui_info.roi.percent_out*radius)));
                    % erode
                    in_ROI = imerode(new_ROI,strel('disk',ceil(gui_info.roi.percent_in*radius)));
                    % join
                    new_ROI = out_ROI;
                    new_ROI(in_ROI) = false;
                end
            otherwise
                return
        end
        % check that the ROI contains at least 15 pixels
        if sum(new_ROI(:))<15
            return
        end
        % if the ROI has to be overwritten, added or removed, change current
        if gui_info.roi.selected>0 && (strcmp(caller_function,'RESEED') || get(h.R_overwrite,'Value') || strcmp(roi_selection_type,'R_exclude') || strcmp(roi_selection_type,'R_include'))
            [pix_y,pix_x]=find(new_ROI);
            roi_info.ROI_xy{gui_info.roi.selected} = [pix_x(:)';pix_y(:)'];
            plot_ROI('change single',gui_info.roi.selected,gui_info.roi.selected)
        else
            % store new ROI in data strctures and add to ROI list
            roi_info.count = roi_info.count +1;
            gui_info.roi.label_count = gui_info.roi.label_count +1;
            [pix_y,pix_x]=find(new_ROI);
            roi_info.ROI_xy{roi_info.count} = [pix_x(:)';pix_y(:)'];
            roi_info.ROI_label = cat(1,roi_info.ROI_label(:),gui_info.roi.label_count);
            roi_info.ROI_type = [roi_info.ROI_type ; {get(get(h.R_ROI_type,'SelectedObject'),'String')}];
            switch roi_info.ROI_type{roi_info.count}
                case 'Neuron'
                    gui_info.roi.list{roi_info.count} = ['<HTML><FONT color="black">ROI ',num2str(gui_info.roi.label_count)];
                case 'Neuropile'
                    gui_info.roi.list{roi_info.count} = ['<HTML><FONT color="green">ROI ',num2str(gui_info.roi.label_count)];
                case 'Blood vessel'
                    gui_info.roi.list{roi_info.count} = ['<HTML><FONT color="blue">ROI ',num2str(gui_info.roi.label_count)];
            end
            % add only the last ROI
            plot_ROI('change single',gui_info.roi.selected,roi_info.count)
            % select latest ROI and update GUI if not overwritten
            gui_info.roi.selected = roi_info.count;
            set(h.R_list,'Value',roi_info.count)
            set(h.R_list,'String',gui_info.roi.list)
            
        end
        % UPDATE
        gui_info.control.roi_export_ready = false;
        gui_update('ROI_exclude')
        gui_update('Export_ready')
        
    end

    function ROI_seed_all(hObject,eventdata)
       % go through all ROIs that are cells and re-seed with center of mass
       
       selected_bk = gui_info.roi.selected;
       
       for roi_ii = 1:roi_info.count
           switch roi_info.ROI_type{roi_ii}
               case {'Neuron','Glial cell'}
                   gui_info.roi.selected = roi_ii;
                   ROI_make_new([],[],'RESEED')
           end
       end
       
       gui_info.roi.selected = selected_bk;
       plot_ROI()
        
    end

    function ROI_shift(hObject,eventdata)
        % shift all the ROIs with the arrows of the keyboard
        % detect pressed arrow

        switch eventdata.Key
            case 'rightarrow'
                for roi_ii=1:roi_info.count
                    roi_info.ROI_xy{roi_ii}(1,:) = roi_info.ROI_xy{roi_ii}(1,:) + 1;
                    set(roi_info.handle_contour{roi_ii},'ZData',circshift(get(roi_info.handle_contour{roi_ii},'ZData'),[0 1]));
                end
            case 'leftarrow'
                for roi_ii=1:roi_info.count
                    roi_info.ROI_xy{roi_ii}(1,:) = roi_info.ROI_xy{roi_ii}(1,:) - 1;
                    set(roi_info.handle_contour{roi_ii},'ZData',circshift(get(roi_info.handle_contour{roi_ii},'ZData'),[0 -1]));
                end
            case 'uparrow'
                for roi_ii=1:roi_info.count
                    roi_info.ROI_xy{roi_ii}(2,:) = roi_info.ROI_xy{roi_ii}(2,:) - 1;
                    set(roi_info.handle_contour{roi_ii},'ZData',circshift(get(roi_info.handle_contour{roi_ii},'ZData'),[-1 0]));
                end
            case 'downarrow'
                for roi_ii=1:roi_info.count
                    roi_info.ROI_xy{roi_ii}(2,:) = roi_info.ROI_xy{roi_ii}(2,:) + 1;
                    set(roi_info.handle_contour{roi_ii},'ZData',circshift(get(roi_info.handle_contour{roi_ii},'ZData'),[1 0]));
                end
            case 'a'
                roi_info.ROI_xy{gui_info.roi.selected}(1,:) = roi_info.ROI_xy{gui_info.roi.selected}(1,:) - 1;
                plot_ROI('change single',gui_info.roi.selected,gui_info.roi.selected)
            case 'd'
                roi_info.ROI_xy{gui_info.roi.selected}(1,:) = roi_info.ROI_xy{gui_info.roi.selected}(1,:) + 1;
                plot_ROI('change single',gui_info.roi.selected,gui_info.roi.selected)
            case 's'
                roi_info.ROI_xy{gui_info.roi.selected}(2,:) = roi_info.ROI_xy{gui_info.roi.selected}(2,:) + 1;
                plot_ROI('change single',gui_info.roi.selected,gui_info.roi.selected)
            case 'w'
                roi_info.ROI_xy{gui_info.roi.selected}(2,:) = roi_info.ROI_xy{gui_info.roi.selected}(2,:) - 1;
                plot_ROI('change single',gui_info.roi.selected,gui_info.roi.selected)
            otherwise
                return
        end
        % mark stuff as not ready to export
        gui_info.control.roi_export_ready = false;
        % update GUI
        gui_update('Export_ready')
        
    end

    function ROI_align()
        % Opens a GIU to find the parameters that align the ROI
        % The parameters to vary are rotation angle (+ center point),
        % scaling and x-y translation.
        roi_info.ROI_xy = segmentation_ROI_alignment(image_info.plotted(:,:,[2 1 3]),roi_info.ROI_xy);
    end
        

    function ROI_process()
        % This function processes the ROI before extracting the signal. It
        % splits shared ROIs, removes the FILTH_ROI, and defines a
        % neuropile ROI to each to normalize signal. A figure with the
        % final ROIs is opened for the user to store it
        % ---
        if gui_info.control.roi_export_ready
            return
        end
        % Sort ROIs to have Neuropile and Blood vessels at the end
        ind = [find(strcmp(roi_info.ROI_type,'Neuron'));find(strcmp(roi_info.ROI_type,'Glial cell'));find(strcmp(roi_info.ROI_type,'Neuropile'));find(strcmp(roi_info.ROI_type,'Blood vessel'))];
        roi_info.ROI_xy = roi_info.ROI_xy(ind);
        roi_info.ROI_label = roi_info.ROI_label(ind);
        roi_info.ROI_type = roi_info.ROI_type(ind);
        gui_info.roi.list = gui_info.roi.list(ind);
        set(h.('R_list'),'String',gui_info.roi.list)
        
        % 0) chop ROIs that move across borders because of shifting
        for roi_ii = 1:roi_info.count
            roi_info.ROI_xy{roi_ii}(1,roi_info.ROI_xy{roi_ii}(1,:)<1) = 1;
            roi_info.ROI_xy{roi_ii}(1,roi_info.ROI_xy{roi_ii}(1,:)>data_info.width) = data_info.width;
            
            roi_info.ROI_xy{roi_ii}(2,roi_info.ROI_xy{roi_ii}(2,:)<1) = 1;
            roi_info.ROI_xy{roi_ii}(2,roi_info.ROI_xy{roi_ii}(2,:)>data_info.height) = data_info.height;
        end
        
        % 1) Split the shared ROIs by assigning the pixels to the cell
        % with the nearest center of mass.
        
        % calculate center of mass [y x radius] of the ROIs for the plots
        center_of_mass = zeros(roi_info.count,3);
        for roi_ii=1:roi_info.count
            ind_y = roi_info.ROI_xy{roi_ii}(2,:);
            ind_x = roi_info.ROI_xy{roi_ii}(1,:);
            center_of_mass(roi_ii,:) = [round(mean(ind_y)),round(mean(ind_x)),...
                mean(sqrt((ind_y-mean(ind_y)).^2+(ind_x-mean(ind_x)).^2))];
        end
        
        % sum all ROIs that are not Neuropile
        sum_ROIs = false(data_info.height,data_info.width);
        for roi_ii = 1:roi_info.count
            sum_ROIs = sum_ROIs + ROI_get_2D(roi_info.ROI_xy{roi_ii});
        end
        
        % find pixels that are in more than one ROI
        [rep_y,rep_x] = find(sum_ROIs>1);
        % assign to closes center of mass
        if ~isempty(rep_y)
            % find center of mass closer to overlapping pixels
            [~,roi_ind]=min(bsxfun(@minus,center_of_mass(:,1)',rep_y).^2+bsxfun(@minus,center_of_mass(:,2)',rep_x).^2,[],2);
            % remove combined pixels from all ROI's
            for roi_ii=1:roi_info.count
                list = ismember(roi_info.ROI_xy{roi_ii}',[rep_x(roi_ind~=roi_ii),rep_y(roi_ind~=roi_ii)],'rows');
                roi_info.ROI_xy{roi_ii}(:,list) = [];
            end
        end
        
        % mark ROIs as ready to export
        gui_info.control.roi_export_ready = true;
        gui_update('Export_ready')
    end

    function ROI_figure(hObject,eventdata,type)
        % make figures of the ROIs
        
        % calculate center of mass of the ROIs for the plots
        center_of_mass = zeros(roi_info.count,3);
        for roi_ii=1:roi_info.count
            ind_y = roi_info.ROI_xy{roi_ii}(2,:);
            ind_x = roi_info.ROI_xy{roi_ii}(1,:);
            center_of_mass(roi_ii,:) = [round(mean(ind_y)),round(mean(ind_x)),...
                mean(sqrt((ind_y-mean(ind_y)).^2+(ind_x-mean(ind_x)).^2))];
        end
        
        switch type
            case 'ID'  % PLOT NUMBERS ON 2D PLOT OF ROIS
                
                % close previous if exists
                close(findobj('type','figure','name','Cell number'))
                scale = 1000/max(data_info.width);
                
                % Current view with ROI labels
                figure('Name','Cell number','Color','white','numbertitle','off',...
                    'position',[50 50 data_info.width*scale data_info.height*scale]);
                % draw current view
                subplot('Position',[0 0 1 1]);
                image(image_info.plotted(:,:,[2 1 3]));
                axis image off% xy
                hold on
                for roi_ii = 1:roi_info.count
                    switch roi_info.ROI_type{roi_ii}
                        case {'Neuron','Glial cell'}
                            text(center_of_mass(roi_ii,2),center_of_mass(roi_ii,1),num2str(roi_info.ROI_label(roi_ii)),'FontSize',12,'Color','y','FontWeight','bold')
                        case 'Neuropile'
                            text(center_of_mass(roi_ii,2),center_of_mass(roi_ii,1),num2str(roi_info.ROI_label(roi_ii)),'FontSize',12,'Color',[0.5 1 0.5],'FontWeight','bold')
                        case 'Blood vessel'
                            text(center_of_mass(roi_ii,2),center_of_mass(roi_ii,1),num2str(roi_info.ROI_label(roi_ii)),'FontSize',12,'Color',[0 0.5 1],'FontWeight','bold')
                    end
                end
                
            case 'ROI'     % PLOT CONTOURS
                
                % close previous if exists
                close(findobj('type','figure','name','ROIs'))
                % Processed ROIs and neuropile
                scale = 1000/max(data_info.width);
                figure('Name','ROIs','Color','white','numbertitle','off',...
                    'position',[50 50 data_info.width*scale data_info.height*scale]);
                subplot('Position',[0 0 1 1]);
                image(image_info.plotted(:,:,[2 1 3]));
                axis image off
                hold on
                % plot background only of cells and otherwise contours
                for roi_ii=1:roi_info.count
                    ROI = ROI_get_2D(roi_info.ROI_xy{roi_ii});
                    switch roi_info.ROI_type{roi_ii}
                        case {'Neuron','Glial cell'}
                            contour(ROI,[-1 1],'LineStyle','-','Color','w','LineWidth',2)
                        case 'Neuropile'
                            contour(ROI,[-1 1],'LineStyle','-','Color',[0.5 1 0.5]*0.5,'LineWidth',2)
                        case 'Blood vessel'
                            contour(ROI,[-1 1],'LineStyle','-','Color',[0 0.5 1]*0.5,'LineWidth',2)
                    end
                end
                
            case 'Background' % PLOT ROIS AND BG
                
                % close previous if exists
                close(findobj('type','figure','name','Processed ROIs'))
                % Processed ROIs and neuropile
                scale = 1000/max(data_info.width);
                figure('Name','Processed ROIs','Color','white','numbertitle','off',...
                    'position',[50 50 data_info.width*scale data_info.height*scale]);
                subplot('Position',[0 0 1 1]);
                image(image_info.plotted(:,:,[2 1 3]));
                axis image off
                hold on
                
                % Make background ROIs
                sum_ROIs = false(data_info.height,data_info.width);
                for roi_ii = 1:roi_info.count
                    if ~strcmp(roi_info.ROI_type{roi_ii},'Neuron') || strcmp(roi_info.ROI_type{roi_ii},'Glial cell')
                        continue
                    end
                    sum_ROIs = sum_ROIs + ROI_get_2D(roi_info.ROI_xy{roi_ii});
                end
                sum_ROIs = imdilate(imfill(logical(sum_ROIs),'holes'),strel('disk',4));
                % makes bg_ROI as squares 2.5x the mean radius and remove filt of
                % ROI in processed_ROI
                mean_radius = round(mean(center_of_mass(:,3)));
                region = [floor(-3*mean_radius):ceil(3*mean_radius)];
                
                % check which ROIs come from cells and plot them
                h_tmp=imagesc(logical(sum_ROIs).*~data_info.filth+2*data_info.filth);
                axis image off
                colormap([0 0 0;1 0 1;0 1 0])
                set(h_tmp, 'AlphaData', 0.2);
                % plot background only of cells and otherwise contours
                sum_ROIs = imdilate(imfill(logical(sum_ROIs),'holes'),strel('disk',4));
                for roi_ii=1:roi_info.count
                    switch roi_info.ROI_type{roi_ii}
                        case {'Neuron','Glial cell'}
                            bg_ROI = false(data_info.height,data_info.width);
                            bg_ROI(...
                                min([data_info.height*ones(size(region));max([ones(size(region));center_of_mass(roi_ii,1)+region])]),...
                                min([data_info.width*ones(size(region));max([ones(size(region));center_of_mass(roi_ii,2)+region])]))=true;
                            bg_ROI = bg_ROI.*~logical(sum_ROIs+data_info.filth);
                            
                            contour(bg_ROI,[-1 1],'LineStyle','-','Color','w','LineWidth',2)
                            text(center_of_mass(roi_ii,2)-2,center_of_mass(roi_ii,1)+2,num2str(roi_info.ROI_label(roi_ii)),'FontSize',10,'Color','y','FontWeight','bold')
                        case 'Neuropile'
                            ROI = ROI_get_2D(roi_info.ROI_xy{roi_ii});
                            contour(ROI,[-1 1],'LineStyle','-','Color',[0.5 1 0.5]*0.5,'LineWidth',2)
                            text(center_of_mass(roi_ii,2)-2,center_of_mass(roi_ii,1)+2,num2str(roi_info.ROI_label(roi_ii)),'FontSize',10,'Color',[0.5 1 0.5],'FontWeight','bold')
                        case 'Blood vessel'
                            ROI = ROI_get_2D(roi_info.ROI_xy{roi_ii});
                            contour(ROI,[-1 1],'LineStyle','-','Color',[0 0.5 1]*0.5,'LineWidth',2)
                            text(center_of_mass(roi_ii,2)-2,center_of_mass(roi_ii,1)+2,num2str(roi_info.ROI_label(roi_ii)),'FontSize',10,'Color',[0 0.5 1],'FontWeight','bold')
                    end
                end
        end
    end

    function ROI = ROI_get_2D(x_y_pix)
        % remove pixels that are outside of the image (shifting)
        x_y_pix(1,x_y_pix(1,:)<1) = 1;
        x_y_pix(1,x_y_pix(1,:)>data_info.width) = data_info.width;
        
        x_y_pix(2,x_y_pix(2,:)<1) = 1;
        x_y_pix(2,x_y_pix(2,:)>data_info.height) = data_info.height;
        
        % construct 2D ROI
        ROI = false(data_info.height,data_info.width);
        ROI(sub2ind([data_info.height,data_info.width],x_y_pix(2,:),x_y_pix(1,:))) = true;
    end
end