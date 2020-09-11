function gui_analyze(ini_path)

if nargin==0
    ini_path =  [pwd,filesep];
end

% === GUI properties:

gui_info = [];
gui_info.control.path = ini_path;

% default processing parameters
gui_info.parameters = [];
gui_info.parameters.NP_fraction = 0.6;
gui_info.parameters.Fo_range = 6;
gui_info.parameters.Fo_quantile = 0.25;
gui_info.parameters.Filter = 0.2;
gui_info.parameters.signal_threshold = 0.15;
gui_info.parameters.probability_threshold = 1;

% data panel
gui_info.data.linked_file = [];
gui_info.data.count = [];
gui_info.data.list = {};
gui_info.data.selected = [];

% graphs
gui_info.graph.stimulus = [];
gui_info.graph.xlim = [0 1];
gui_info.graph.ylim_raw = [0 1];
gui_info.graph.ylim_dFF = [0 1];
gui_info.graph.ylim_prob = [0 1];
gui_info.graph.manipulation = 'X-Zoom';
gui_info.graph.limit = 'current';

% plot options
gui_info.analysis.reprocess = true;
gui_info.analysis.add_dFF = true;
gui_info.analysis.activity_maximum = 6;
gui_info.analysis.activity_bins = 12;
gui_info.analysis.correlation_minimum = 0;
gui_info.analysis.correlation_maximum = 1;
gui_info.analysis.participation_bin_sec = 2;

% === data + ROI

data_info = [];

% common to all the set
data_info.handle = [];
data_info.source_file = [];
data_info.num_stacks = [];
data_info.frame_rate = [];
data_info.ROI_label = [];
data_info.ROI_type = {};
data_info.ROI_x_y = [];
data_info.range_raw = [];
data_info.range_processed = [];

% for the individual ROI
data = [];

% === analysis results

analysis_info = [];

analysis_info.calculate = true;
analysis_info.all_events = [];
analysis_info.correlation = [];
analysis_info.ch_ratio = [];
analysis_info.participation_rate = [];

%% INITIATE

% ---> draw the interface with callbacks and get handles
[hfig,h]=gui_create();

% load first set
callback_data(h.D_load,[])
gui_resize()

%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               V I E W
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------- GUI
    function [hfig,h]=gui_create()
        
        % Draw the initial GUI using the results of GUIDE
        % hfig: main figure hangle
        % h   : structure with controls handles
        %hfig = process_gui_draw_build;
        sourcepath = mfilename('fullpath');
        idx = strfind(sourcepath,filesep);
        hfig = open([sourcepath(1:idx(end)),'gui_analyze.fig']);
        h = guihandles(hfig);
        
        % General figure
        set(hfig,{'Name','WindowScrollWheelFcn','ResizeFcn'},{'Signal Analysis GUI',@gui_scroll,@gui_resize})
        % axes to graph extracted signal
        set(h.A_raw,{'HandleVisibility','ButtonDownFcn','NextPlot','TickDir','xTickLabel','yTick','xlim','ylim','Box'},...
            {'on',@gui_click,'replacechildren','out',{},0:0.1:100,[0 1],[0 1],'on'})   
        set(h.A_dFF,{'HandleVisibility','ButtonDownFcn','NextPlot','TickDir','xTickLabel','yTick','xlim','ylim','Box'},...
            {'on',@gui_click,'replacechildren','out',{},0:0.1:100,[0 1],[0 1],'on'})
        set(h.A_prob,{'HandleVisibility','ButtonDownFcn','NextPlot','TickDir','xTick','xTickLabel','yTick','xlim','ylim','Box'},...
            {'on',@gui_click,'replacechildren','out',[0:100],{0:100},0:0.1:100,[0 1],[0 1],'on'})
        
        % Change GUI properties
        % ------------------------------------------------   Context menu
        % :: Make context to add stimulus info to signal
        h.M_data = uimenu(hfig,'Label','Data');
        mh1 = uimenu(h.M_data,'Label','Linked file');
        uimenu(mh1,'Tag','D_link_add','Label','Add','Callback',@callback_uimenu,'Enable','on');
        uimenu(mh1,'Tag','D_link_delete','Label','Delete','Callback',@callback_uimenu,'Enable','on');
        uimenu(h.M_data,'Tag','D_info','Label','Show info','Callback',@callback_uimenu);
        uimenu(h.M_data,'Label','Export','Tag','export_excel','Enable','on','Callback',@callback_uimenu);
        
        % :: Make context menu for extra plots
        hm0 = uimenu(hfig,'Label','Plot');
        uimenu(hm0,'Tag','show_ROI','Label','ROIs','Callback',@callback_uimenu);
        uimenu(hm0,'Tag','current_dFF','Label','Current dF/F','Callback',@callback_uimenu);
        uimenu(hm0,'Tag','all_dFF','Label','All dF/F','Callback',@callback_uimenu);
        
        % :: Make context menu for analysis
        hm0 = uimenu(hfig,'Label','Analysis');
        
        mh1 = uimenu(hm0,'Label','Signal correlation','Enable','on');
        uimenu(mh1,'Label','Show','Enable','on','Callback',{@analysis_show,'correlation'})
        uimenu(mh1,'Label','Parameters','Callback',{@analysis_parameters,'correlation_parameters'})
        uimenu(mh1,'Label','Colorbar','Enable','on','Callback',{@analysis_parameters,'correlation_colorbar'})
        
        mh1 = uimenu(hm0,'Label','Activity map','Enable','on');
        uimenu(mh1,'Label','Show','Enable','on','Callback',{@analysis_show,'activity'});
        uimenu(mh1,'Label','Parameters','Callback',{@analysis_parameters,'activity_parameters'})
        uimenu(mh1,'Label','Colorbar','Enable','on','Callback',{@analysis_parameters,'activity_colorbar'})
        
        mh1 = uimenu(hm0,'Label','Activity raster','Enable','on');
        uimenu(mh1,'Label','Show','Enable','on','Callback',{@analysis_show,'raster'})
        uimenu(mh1,'Label','Mark selected dF/F','Checked','on','Callback',{@analysis_parameters,'raster'});
        
        mh1 = uimenu(hm0,'Label','Participation rate','Enable','on');
        uimenu(mh1,'Label','Show','Enable','on','Callback',{@analysis_show,'Participation rate'})
        uimenu(mh1,'Label','Parameters','Enable','on','Callback',{@analysis_parameters,'participation'})
        
        h.R_ratio = uimenu(hm0,'Label','Channel ratio','Enable','off','Callback',{@analysis_show,'Channel ratio'});
        
        %----------------------------------------------   Data panel  -> D
        
        set(h.D_load,'CallBack',@callback_data)
        set(h.D_list,{'CallBack','String'},{@callback_data,''})
        set(h.D_update,{'CallBack'},{@callback_data})
        set(h.D_keep,'Value',1)
        
        %----------------------------------------------- Parameters panel  ->S
        % Re-process button
        set(h.S_reprocess,{'Value','Callback'},{gui_info.analysis.reprocess,@callback_signal})
        
        % Neuropile panel
        set(h.S_NP_ROI,{'Enable','CallBack'},{'off',@callback_signal})
        set(h.S_NP_fraction,{'String','CallBack'},{num2str(gui_info.parameters.NP_fraction),@callback_signal})
        % Fo and filter panel
        set(h.S_Filter,{'String','CallBack'},{num2str(gui_info.parameters.Filter),@callback_signal})
        set(h.S_Fo_range,{'String','CallBack'},{num2str(gui_info.parameters.Fo_range),@callback_signal})
        set(h.S_Fo_quantile,{'String','CallBack'},{num2str(gui_info.parameters.Fo_quantile),@callback_signal})
        % event detection panel
        set(h.S_signal_threshold,{'String','CallBack'},{num2str(gui_info.parameters.signal_threshold),@callback_signal})
        set(h.S_probability_threshold,{'String','CallBack'},{num2str(gui_info.parameters.probability_threshold),@callback_signal})
        
    end

    function gui_update(updateTag,value)
        if nargin==1
            value = true;
        end
        
        switch updateTag
            case 'parameters'
                % when new values are set and the GUI has to be updated
                list = {'NP_fraction','Fo_range','Fo_quantile','Filter','signal_threshold','probability_threshold'};
                for l_ii = 1:length(list)
                    set(h.(['S_',list{l_ii}]),'String',num2str(gui_info.parameters.(list{l_ii})));
                end
                
            case 'ROI_list'
                % display a list of ROIs and select last one
                gui_info.data.count = length(data_info.ROI_label);
                gui_info.data.list = cell(1,gui_info.data.count);
                for ind = 1:gui_info.data.count
                    switch data_info.ROI_type{ind}
                        case {'Neuron','Glial cell'}
                            gui_info.data.list{ind} = ['<HTML><FONT color="black">ROI ',num2str(data_info.ROI_label(ind))];
                        case 'Neuropile'
                            gui_info.data.list{ind} = ['<HTML><FONT color="green">ROI ',num2str(data_info.ROI_label(ind))];
                        case 'Blood vessel'
                            gui_info.data.list{ind} = ['<HTML><FONT color="blue">ROI ',num2str(data_info.ROI_label(ind))];
                    end
                end
                % update GUI
                gui_info.data.selected = 1;
                set(h.D_list,'String',gui_info.data.list)
                set(h.D_list,'Value',gui_info.data.selected)
                % activate BV fraction button
                if sum(strcmp(data_info.ROI_type,'Blood vessel'))>0
                    set(h.S_NP_ROI,'Enable','on')
                else
                    set(h.S_NP_ROI,'Enable','off')
                end
                
            case 'axes'
                axes(h.A_raw)
               % ylabel('Raw signal')
                xlim(gui_info.graph.xlim)
                set(gca,'XTick',0:round(gui_info.graph.xlim(end)))
                ylim(gui_info.graph.ylim_raw)
                set(gca,'YTick',0:10:round(gui_info.graph.ylim_raw(end)))
                
                axes(h.A_dFF)
               % ylabel('dF/F')
                xlim(gui_info.graph.xlim)
                set(gca,'XTick',0:round(gui_info.graph.xlim(end)))
                ylim(gui_info.graph.ylim_dFF)
                set(gca,'YTick',-1:0.1:round(gui_info.graph.ylim_dFF(end)))
                
                axes(h.A_prob)
               % ylabel('Event likelihood')
               % xlabel('time [min]')
                xlim(gui_info.graph.xlim)
                set(gca,'XTick',0:round(gui_info.graph.xlim(end)))
                ylim(gui_info.graph.ylim_prob)
                set(gca,'YTick',0:round(gui_info.graph.ylim_prob(end)))
                
            case 'update_color'
                if value
                    set(h.D_update,'BackgroundColor',[0.9333 0.8667 0.5098])
                else
                    set(h.D_update,'BackgroundColor',[1 1 1]*0.941)
                end
                
            case 'signal_panel'
                 if gui_info.analysis.reprocess
                    set(h.P_signal,'visible','on')
                 else
                    set(h.P_signal,'visible','off')
                end
        end
    end

    function gui_scroll(hObject,eventdata)
        
        axis_to_check = {'raw','dFF','prob'};
        zoom_factor = 0.1;
        
        % start with empty cursor location
        cursor_location = [];
        % check all graphs to get current point
        for ax=1:3
            if ~isempty(cursor_location)
                break
            end
            ylm = gui_info.graph.(['ylim_',axis_to_check{ax}]);
            tmp = get(h.(['A_',axis_to_check{ax}]),'CurrentPoint');
            if tmp(1,1)>=gui_info.graph.xlim(1) && tmp(1,1)<=gui_info.graph.xlim(2) && tmp(1,2)>=ylm(1) && tmp(1,2)<=ylm(2) && tmp(1,3)==1
                cursor_location = [tmp([1 3]),ax];
            end
        end
        % if scroll outside of axes, return
        if isempty(cursor_location)
            return
        end
        
        % check type of zoom
        switch gui_info.graph.manipulation
            case 'X-Zoom'
                % get current limit
                new_width = diff(gui_info.graph.xlim);
                % measure the scroll
                if eventdata.VerticalScrollCount < 0
                    % scroll up: ZOOM IN TO POINTER
                    new_width = (1-zoom_factor)*diff(gui_info.graph.xlim);
                elseif eventdata.VerticalScrollCount > 0
                    % scroll down: ZOOM OUT TO POINTER
                    new_width = (1+zoom_factor)*diff(gui_info.graph.xlim);
                    % if larger than original, set original
                    if new_width>data_info.num_stacks/data_info.frame_rate/60
                        new_width = data_info.num_stacks/data_info.frame_rate/60;
                    end
                end
                % center around cursor
                fraction = (cursor_location(1)-gui_info.graph.xlim(1))/(gui_info.graph.xlim(2)-gui_info.graph.xlim(1));
                new_xlim = [cursor_location(1)-fraction*new_width,cursor_location(1)+(1-fraction)*new_width];
                % shift x for not going outside graph
                if new_xlim(1)<=0
                    new_xlim = new_xlim -new_xlim(1);
                elseif new_xlim(2)>data_info.num_stacks/data_info.frame_rate/60
                    new_xlim = new_xlim + (data_info.num_stacks/data_info.frame_rate/60-new_xlim(2));
                end
                % set limits
                gui_info.graph.xlim = new_xlim;
                
            case 'Y-Zoom'
                % get current limit
                ylim = gui_info.graph.(['ylim_',axis_to_check{cursor_location(3)}]);
                new_height = diff(ylim);
                % measure the scroll
                if eventdata.VerticalScrollCount < 0
                    % scroll up: ZOOM IN TO POINTER
                    new_height = (1-zoom_factor)*new_height;
                elseif eventdata.VerticalScrollCount > 0
                    % scroll down: ZOOM OUT TO POINTER
                    new_height = (1+zoom_factor)*new_height;
                end
                % center around cursor
                fraction = (cursor_location(2)-ylim(1))/(ylim(2)-ylim(1));
                new_ylim = [cursor_location(2)-fraction*new_height,cursor_location(2)+(1-fraction)*new_height];
                % for probability axes, bound to zero
                if cursor_location(3)==3
                    new_ylim = new_ylim-new_ylim(1);
                end
                % set limits and update figure
                gui_info.graph.(['ylim_',axis_to_check{cursor_location(3)}]) = new_ylim;
        end
        
        % update the axes of the plot
        gui_update('axes')
    end

    function gui_click(hObject,eventdata)
        % control graphs uisng mouse clicks (scrolling in another function!)
        switch get(hfig,'selectiontype')
            case {'normal','open'}
                % == switch between limits for current data or for global
                if strcmp(gui_info.graph.limit,'current');
                    plot_range('fit all')
                else
                    plot_range('fit current')
                end
                
            case {'alt','extend'}
                % === switch between X and Y zoom
                if strcmp(gui_info.graph.manipulation,'X-Zoom')
                    gui_info.graph.manipulation = 'Y-Zoom';
                else
                    gui_info.graph.manipulation = 'X-Zoom';
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
        
        % Figure 1 panel position
        handles = [h.A_prob,h.A_dFF,h.A_raw];
        y_start = 3*blank_space;
        y_size = (pos_fig(4)-7*blank_space)/3;
        for f_ii = 1:3
            set(handles(f_ii),'Position',...
                [ini_pos_x + 6*blank_space, y_start, pos_fig(3)-ini_pos_x - 8*blank_space, y_size]);
            y_start = y_start + y_size + blank_space;
        end
    end

% ------------------------------------------------------------------- PLOTS
    function plot_graph()
        
        time_axis = (1:data_info.num_stacks)/data_info.frame_rate/60;
        
        %====  RAW signal
        axes(h.A_raw)
        cla
        plot_stim([0 300])
        hold on
        plot(time_axis,data.signal(1,:,1),'Color',[1 1 1]*0.3)
        if ~isempty(data.neuropile)
            plot(time_axis,data.neuropile(1,:,1),'Color',[0 0 1]*0.8)
        end
        ylabel('Raw signal')
        
        %====  Show thresholds
        axes(h.A_dFF)
        cla
        plot_stim([-1 20])
        hold on
        plot([0 time_axis(end)],[1 1]*gui_info.parameters.signal_threshold,'--k')
        ylabel('dFF')
        
        axes(h.A_prob)
        cla
        plot([0 time_axis(end)],[1 1]*gui_info.parameters.probability_threshold,'--k')
        ylabel('Likelihood')
        xlabel('time [min]')
        
        if strcmp(data.type,'Neuropile') || strcmp(data.type,'Blood vessel')
            return
        end
        
        %====  Processed signal + events
        axes(h.A_dFF)
        hold on
        plot(time_axis,data.processed(1,:,1),'Color',[1 1 1]*0.3)
        if ~isempty(data.events)
            plot(time_axis(data.events(:,1)),data.processed(data.events(:,1)+1),'ko','MarkerFaceColor','g')
        end

        %==== Probability of event
        axes(h.A_prob)
        hold on
        if isfield(data,'signal_bumps')
            plot([time_axis(data.signal_bumps.frame);time_axis(data.signal_bumps.frame)],...
                [zeros(size(data.signal_bumps.prob)) data.signal_bumps.prob]','Color',[1 1 1]*0.3)
        end
        
    end

    function plot_stim(plot_limits)
        if isempty(gui_info.graph.stimulus)
            return
        end
        % convert stimulus into minutes
        stim = gui_info.graph.stimulus(:,1:2)/60;
        for stim_ii = 1:size(gui_info.graph.stimulus,1)
            fill(...
                [stim(stim_ii,1) stim(stim_ii,2) stim(stim_ii,2) stim(stim_ii,1)],...
                [plot_limits(1) plot_limits(1) plot_limits(2) plot_limits(2)],...
                'r','FaceColor',[1 0 0],'FaceAlpha',0.2,'EdgeColor',0.6*[1 1 1]);
        end
        
    end

    function plot_range(type)
        switch type
            case 'fit current'
                gui_info.graph.limit = 'current';
                % change limits of graph given the data
                gui_info.graph.xlim = [0 data_info.num_stacks/data_info.frame_rate/60];
                gui_info.graph.ylim_raw = [min([data.signal(1,:,1) data.neuropile]) max([data.signal(1,:,1) data.neuropile])];
                gui_info.graph.ylim_dFF = [min(data.processed) max([data.processed gui_info.parameters.signal_threshold+0.1])];
                if isfield(data,'signal_bumps')
                    gui_info.graph.ylim_prob = [0 max(data.signal_bumps.prob)];
                else
                    gui_info.graph.ylim_prob = [0 2*gui_info.parameters.probability_threshold];
                end
                gui_update('axes')
                
            case 'fit all'
                gui_info.graph.limit = 'all';
                % update the axes of the plot
                gui_info.graph.xlim = [0 data_info.num_stacks/data_info.frame_rate/60];
                gui_info.graph.ylim_raw = data_info.range_raw;
                gui_info.graph.ylim_dFF = [data_info.range_processed(1) max([data_info.range_processed(2) gui_info.parameters.signal_threshold])];
                if isfield(data,'signal_bumps')
                    gui_info.graph.ylim_prob = [0 max(data.signal_bumps.prob)];
                else
                    gui_info.graph.ylim_prob = [0 2*gui_info.parameters.probability_threshold];
                end
                gui_update('axes')
        end
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          C O N T R O L L E R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function callback_uimenu(hObject,eventdata)
        
        gui_tag = get(hObject,'Tag');
        
        switch gui_tag
            
            % === Data
            case 'D_link_add'
                % Add linked file
                [FileName,PathName] = uigetfile([gui_info.control.path,'*.mat'],'MultiSelect','off','Select stimulus file to open');
                gui_info.graph.stimulus = process_stim_file([PathName,FileName],data_info.num_stacks,data_info.frame_rate);
                plot_graph()
                
            case 'D_link_delete'
                % Delete linked file
                gui_info.graph.stimulus = [];
                plot_graph()
                
            case 'D_info'
                % show info
                d_info = load(data_info.source_file);
                % display current set info
                fname = d_info.data_class.files(:);
                fname = [cell(2,1);fname];
                fname{1} = 'Files: ';
                params = data_info.handle.parameters;
                txt = [fname;{...
                    ['Number of channels: ',num2str(size(data.signal,3))];...
                    ['Stacks: ',num2str(data_info.num_stacks)];...
                    ['Frame rate: ',num2str(data_info.frame_rate)]};...
                    ['  '];...
                    ['Analysis paramaters'];...
                    ['  '];...
                    ['NP correction: ',num2str(params.NP_fraction)];...
                    ['Fo range: ',num2str(params.Fo_range)];...
                    ['Fo quantile: ',num2str(params.Fo_quantile)];...
                    ['Exp.Filter constant: ',num2str(params.Filter)];...
                    ['Event threshold signal: ',num2str(params.signal_threshold)];...
                    ['Event threshold likelihood: ',num2str(params.probability_threshold)];...
                    ];
                msgbox(txt,'SET INFO')
                drawnow
                
            case 'show_ROI'
                % display the position of the ROIs
                
                % load image
                load(data_info.source_file,'data_img');
                % close previous if exists
                close(findobj('type','figure','name','ROIs'))
                scale = 1000/size(data_img,2);
                
                % Current view with ROI labels
                figure('Name','ROIs','Color','white','numbertitle','off',...
                    'position',[50 50 size(data_img,2)*scale size(data_img,1)*scale]);
                subplot('Position',[0 0 1 1]);
                image(0.75*data_img(:,:,[2 1 3]));
                axis image off% xy
                hold on
                
                % Label the ROIs
                for roi_num = 1:gui_info.data.count
                    plot(data_info.ROI_x_y(roi_num,1)+2,data_info.ROI_x_y(roi_num,2),'o','LineWidth',1,'Color','w')
                    switch data_info.ROI_type{roi_num}
                        case {'Neuron','Glial cell'}
                            text(data_info.ROI_x_y(roi_num,1)+5,data_info.ROI_x_y(roi_num,2),num2str(data_info.ROI_label(roi_num)),'FontSize',12,'Color','y','FontWeight','bold')
                        case 'Neuropile'
                            text(data_info.ROI_x_y(roi_num,1)+5,data_info.ROI_x_y(roi_num,2),num2str(data_info.ROI_label(roi_num)),'FontSize',12,'Color',[0.5 1 0.5],'FontWeight','bold')
                        case 'Blood vessel'
                            text(data_info.ROI_x_y(roi_num,1)+5,data_info.ROI_x_y(roi_num,2),num2str(data_info.ROI_label(roi_num)),'FontSize',12,'Color',[0 0.5 1],'FontWeight','bold')
                    end
                end
                
            case 'current_dFF'
                % plot current view
                figure('Name',['ROI ',num2str(data_info.ROI_label(gui_info.data.selected))],'Color','white','numbertitle','off',...
                    'position',[50 50 1000 300]);
                subplot('Position',[0.075 0.15 0.9 0.75])
                time_axis = (1:data_info.num_stacks)/data_info.frame_rate/60;
                hold on
                xlim(gui_info.graph.xlim)
                plot_stim([-1 20])                
                plot(time_axis,data.processed(1,:,1),'Color',[1 1 1]*0.3)
                if ~isempty(data.events)
                    plot(time_axis(data.events(:,1)),data.processed(data.events(:,1)+1),'ko','MarkerFaceColor','g')
                end
                ylabel('dF/F')
                xlabel('time [min]')
                title(['ROI ',num2str(data_info.ROI_label(gui_info.data.selected))])
                set(gca,'XTick',1:round(time_axis(end)))
                ylim(gui_info.graph.ylim_dFF)
                set(gca,'YTick',-1:0.1:20)
                set(gca,'TickDir','out','TickLength',[0.005 0.005])
                box on

            case 'all_dFF'

                % make a figure with all ROIs
                figure('Name','All ROIs','Color','white','numbertitle','off',...
                    'position',[50 50 1000 20*max(data_info.ROI_label)]);
                subplot('Position',[0.075 0.1 0.9 0.875])
                time_axis = (1:data_info.num_stacks)/data_info.frame_rate/60;
                 xlim([0 max(time_axis)])
                 hold on
                plot_stim([-10 max(data_info.ROI_label)+10])
                % grid lines
                plot(repmat(0:time_axis(end),[2 1]),repmat([-1;max(data_info.ROI_label)+5],[1 length(0:time_axis(end))]),'-','LineWidth',0.5,'Color',[1 1 1]*0.9)
                % go through ROIs
                maxVal = 0;
                for roi_num = 1:gui_info.data.count
                    if any(strcmp(data_info.ROI_type{roi_num},{'Neuron','Glial cell'}))
                        % load data, process and plot
                        data = data_info.handle.data(1,roi_num);
                        process_data(true,false)
                        trace = data_info.ROI_label(roi_num)+data.processed(1,:,1)/(3*gui_info.parameters.signal_threshold);
                        maxVal = max([maxVal,max(trace)]);
                        plot(time_axis,trace ,'Color',[1 1 1]*0.3)
                        drawnow
                    end
                end
                % return to marked loaded
                data = data_info.handle.data(1,gui_info.data.selected);
                % beautify
                ylabel('ROI')
                xlabel('time [min]')
                set(gca,'XTick',1:round(time_axis(end)))
                ylim([0 ceil(maxVal)+1])
                set(gca,'YTick',1:max(data_info.ROI_label))
                %set(gca,'YTickLabel',{''})
                set(gca,'TickDir','out','TickLength',[0.005 0.005])
                box on

            case 'export_excel'
                % export current file to excel
                signal_export({data_info.handle.Properties.Source})
        end
    end

    function callback_data(hObject,eventdata)
        gui_tag = get(hObject,'Tag');
        switch gui_tag
            case 'D_load'
                % load data set
                [FileName,PathName] = uigetfile([gui_info.control.path,'*_extracted.mat'],'MultiSelect','off','Select file to open');
                gui_info.control.path = PathName;
                % load common info
                tmp = load([PathName,FileName],'id_type_label_pos','parameters','source_file','data_range');
                data_info.ROI_type = tmp.id_type_label_pos(:,2);
                data_info.ROI_label = cell2mat(tmp.id_type_label_pos(:,3));
                data_info.ROI_x_y = cell2mat(tmp.id_type_label_pos(:,4));
                gui_update('ROI_list')
                % load parameters and show parameters
                if get(h.D_keep,'Value')==1
                    gui_info.parameters = tmp.parameters;
                end
                gui_update('parameters')
                data_info.source_file = tmp.source_file;
                % get range for plots
                data_info.range_raw = [tmp.data_range.raw.min tmp.data_range.raw.max];
                data_info.range_processed = [tmp.data_range.processed.min tmp.data_range.processed.max];
                % prepare data handle and read
                data_info.handle = matfile([PathName,FileName],'Writable',true);
                data = data_info.handle.data(1,gui_info.data.selected);
                % get frame rate and num stacks
                data_info.frame_rate = data.frame_rate;
                data_info.num_stacks = length(data.signal);
                % if two channels, activate figure option
                if size(data.signal,3)>1
                    set(h.R_ratio,'Enable','on')
                else
                    set(h.R_ratio,'Enable','off')
                end
                % do plots
                gui_info.graph.stimulus = [];
                process_data(true,true)
                plot_range('fit current')
                plot_graph()
                % erase analysis results
                analysis_info.calculate = true;
                
            case 'D_list'
                % Get selected ROI and load data
                gui_info.data.selected = get(hObject,'Value');
                data = data_info.handle.data(1,gui_info.data.selected);
                % do plots
                process_data(true,true)
                plot_graph()
                
            case 'D_update'
                % calculate all events and save
                gui_update('update_color',true)
                drawnow
                disp('====== UPDATING ALL ROIs WITH GIVEN PARAMETERS')
                for roi_ii=1:gui_info.data.count
                    data = data_info.handle.data(1,roi_ii);
                    process_data(true,true)
                    if isfield(data,'signal_bumps')
                        data = rmfield(data,'signal_bumps');
                    end
                    data_info.handle.data(1,roi_ii)=data;
                end
                % save the parameters
                data_info.handle.parameters = gui_info.parameters;
                % load the previously selected ROI
                data = data_info.handle.data(1,gui_info.data.selected);
                disp('====== DONE!')
                disp('====== CALCULATING PROPERTIES')
                analysis_info.calculate = true;
                tmp = gui_info.analysis.reprocess;
                gui_info.analysis.reprocess = false;
                analysis_calculate(true) % waitbar is off
                gui_info.analysis.reprocess = tmp;
                time_stamp = now;
                save(data_info.handle.Properties.Source,'analysis_info','time_stamp','-append')
                disp('====== DONE!')
                gui_update('update_color',false)
        end
    end

    function callback_signal(hObject,eventdata)
        gui_tag = get(hObject,'Tag');
        switch gui_tag
            
            case 'S_reprocess'
                % re_calculate properties for figures
                if get(h.S_reprocess,'Value')==0
                    gui_info.analysis.reprocess = false;
                else
                    gui_info.analysis.reprocess = true;
                end
                % delete previous results
                analysis_info.calculate = true;
                % change panel access
                gui_update('signal_panel')
            
            case 'S_NP_ROI'
                % determine fraction for neuropile correction using ROIs
                % around blood vessels. This option is only activated when
                % there are such defined ROIs and the signal was already
                % extracted
                ROI_BV = find(strcmp(data_info.ROI_type,'Blood vessel'));
                % sum the activity divided by the neuropile
                values = zeros(length(ROI_BV),1);
                for roi_ii = 1:length(ROI_BV)
                    data = data_info.handle.data(1,ROI_BV(roi_ii));
                    values(roi_ii) = mean(squeeze(data.signal(1,:,1)./data.neuropile(1,:,1)));
                end
                % return to old data
                data = data_info.handle.data(1,gui_info.data.selected);
                % set mean as value and display
                gui_info.parameters.NP_fraction = round(100*mean(values))/100;
                set(h.S_NP_fraction,'String',num2str(gui_info.parameters.NP_fraction));
                % do plots
                process_data(true,true)
                plot_range('fit current')
                plot_graph()
                
            case {'S_NP_fraction','S_Fo_range','S_Fo_quantile','S_Filter'}
                gui_info.parameters.(gui_tag(3:end)) = max([0 str2double(get(h.(gui_tag),'String'))]);
                set(h.(gui_tag),'String',num2str(gui_info.parameters.(gui_tag(3:end))));
                % do plots
                process_data(true,true)
                plot_range('fit current')
                plot_graph()
                
            case {'S_signal_threshold','S_probability_threshold'}
                gui_info.parameters.(gui_tag(3:end)) = str2double(get(h.(gui_tag),'String'));
                set(h.(gui_tag),'String',num2str(gui_info.parameters.(gui_tag(3:end))));
                % do plots
                process_data(false,true)
                plot_range('fit current')
                plot_graph()
                
            otherwise
                warning('Signal callback type not found!')
                return
        end
        
        % if not read from stored data, empty results
        if gui_info.analysis.reprocess
            analysis_info.calculate = true;
        end
        
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              M O D E L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function process_data(do_process,do_events)
        % process data
        if ~(strcmp(data.type,'Neuropile') || strcmp(data.type,'Blood vessel'))
            
            % re process data
            if gui_info.analysis.reprocess
            
            % calculate with new parameters and store
            [signal_processed,signal_events,signal_bumps] = signal_process(data,gui_info.parameters,do_process,do_events);
            if do_process
                data.processed = signal_processed;
            end
            if do_events
                data.events = signal_events;
                data.signal_bumps = signal_bumps;
            end
            
            else
                % remove signal_bumps to avoid plotting
                if isfield(data,'signal_bumps')
                    data = rmfield(data,'signal_bumps');
                end
            end
        end
    end

    function analysis_calculate(waitbar_off)
        if nargin==0
           waitbar_off = false; 
        end
        
        % Read the data again and process it if required
        if analysis_info.calculate
            if ~waitbar_off
                h_waitbar = waitbar(0,'Reading data');
                drawnow
            end
            
            % prepare variables
            isCell = cellfun(@(x) sum(strcmp(x,{'Neuron','Glial cell'}))>0,data_info.ROI_type);
            analysis_info.all_events = cell(gui_info.data.count,1);
            all_signal_processed = zeros(gui_info.data.count,length(data.signal(1,:,1)));
            if size(data.signal,3)>1
                all_signal_raw = zeros(sum(isCell),size(data.signal,2),size(data.signal,3));
            end
            roi_ii=0;
            
            % go through ROIs
            
            for roi_num = 1:gui_info.data.count
                if ~waitbar_off
                    h_waitbar = waitbar(roi_num/gui_info.data.count,h_waitbar);
                end
                % load data and process
                data = data_info.handle.data(1,roi_num);
                if gui_info.analysis.reprocess
                    process_data(true,true)
                end
                % save
                switch data_info.ROI_type{roi_num}
                    case {'Neuron','Glial cell'}
                        all_signal_processed(roi_num,:) = data.processed;
                        analysis_info.all_events{roi_num} = data.events;
                        roi_ii=roi_ii+1;
                        if size(data.signal,3)>1
                            all_signal_raw(roi_ii,:,:) = data.signal;
                        end
                    case {'Neuropile','Blood vessel'}
                        all_signal_processed(roi_num,:) = data.signal(1,:,1);
                end
            end
            if ~waitbar_off
                close(h_waitbar)
            end
            % return to marked loaded
            data = data_info.handle.data(1,gui_info.data.selected);
            process_data(true,true)
            
            % = CORRELATION
            analysis_info.correlation = corrcoef(all_signal_processed');
            analysis_info.correlation(eye(size(analysis_info.correlation))==1)=0;
            clear all_signal_processed
            
            % == CHANNEL RATIO
            % calculate a Fo of the raw data
            analysis_info.ch_ratio = [];
            if size(data.signal,3)>1
                min_size_half = round((gui_info.parameters.Fo_range*data_info.frame_rate)/2);
                Fo = zeros(sum(isCell),size(data.signal,2),size(data.signal,3));
                for frame=1:data_info.num_stacks
                    tmp = sort(all_signal_raw(:,max(1,frame-min_size_half) : min(data_info.num_stacks,frame+min_size_half),:),2);
                    Fo(:,frame,:) = tmp(:,ceil(gui_info.parameters.Fo_quantile*size(tmp,2)),:);
                end
                clear all_signal_raw
                
                % get average of first 10 seconds
                average_length = ceil(10*data_info.frame_rate);
                if data_info.num_stacks>=average_length
                    average_length = data_info.num_stacks;
                end
                Fo = [mean(Fo(:,1:average_length,:),2) , mean(Fo(:,end:-1:(average_length-1),:),2)];
                % get ratio between first two channels
                Fo = Fo(:,:,1)./Fo(:,:,2);
                % add label and sort
                analysis_info.ch_ratio = sortrows([data_info.ROI_label(isCell),Fo],1);
            end
        end
        
        % The values in Participation rate depend on the parameters. Update
        % it if changed
        if analysis_info.calculate || isempty(analysis_info.participation_rate)
            % When neuron N spikes, how many of the other neurons also
            % spike? And, when any neuron M spikes, does N also spike?
     
            isCell = cellfun(@(x) sum(strcmp(x,{'Neuron','Glial cell'}))>0,data_info.ROI_type);
            % bin the events in seconds window
            events_bin = [];
            range = 1:gui_info.analysis.participation_bin_sec*data_info.frame_rate:data_info.num_stacks;
            for ii=find(isCell)'
                %length((histc(signal_info.events{tmp}(:,1),range)>0))
                if isempty(analysis_info.all_events{ii})
                    tmp = zeros(size(range));
                else
                    tmp = histc(analysis_info.all_events{ii}(:,1),range)>0;
                end
                events_bin = [events_bin;tmp(:)'];
            end
            % calculate participation ratio
            analysis_info.part_rate = zeros(sum(isCell),3);
            analysis_info.part_rate(:,1) = data_info.ROI_label(isCell);
            for ii = 1:sum(isCell)
                % get two vectors for the neuron of interest and the
                % rest
                events_bin_center = events_bin(ii,:);
                events_bin_rest = events_bin;
                events_bin_rest(ii,:)=[];
                events_bin_rest = mean(events_bin_rest,1);
                % calculate ratios
                analysis_info.part_rate(ii,2) = 100*mean(events_bin_rest(events_bin_center>0));
                analysis_info.part_rate(ii,3) = 100*mean(events_bin_center(events_bin_rest>0));
            end
            % sort the rate
            analysis_info.part_rate = sortrows(analysis_info.part_rate,1);
        end
        
        % change status such that it is not calculated again
        analysis_info.calculate = false;
    end

    function analysis_show(hObject,eventdata,figure_type)
        % make figures with the results of the data
        
        % calculate results
        analysis_calculate()
        
        % plot
        switch figure_type
            
            case 'activity' % ===== display events per minute in ROIs
                
                % make color map
                type = 'bluewhitered';
                num_bins = gui_info.analysis.activity_bins;
                cm = getColorMap(type,num_bins);
                
                % divide range in parts of colormap
                activity_range = linspace(0,gui_info.analysis.activity_maximum,gui_info.analysis.activity_bins);
                
                % go through ROIs with cells, calculate activity
                spikes_min = zeros(gui_info.data.count,1);
                for roi_num = 1:gui_info.data.count
                    if ~strcmp(data_info.ROI_type{roi_num},{'Neuron','Glial cell'})
                        continue
                    end
                    % load data and extract events
                    spikes_min(roi_num) = size(analysis_info.all_events{roi_num},1)/(data_info.num_stacks/data_info.frame_rate/60);
                end
                
                % make empty image for ROIs
                load(data_info.source_file,'data_img','ROI');
                im2plot = zeros(size(data_img));
                for roi_num = 1:gui_info.data.count
                    if ~any(strcmp(data_info.ROI_type{roi_num},{'Neuron','Glial cell'}))
                        continue
                    end
                    % find corresponding range
                    ind_range = find(activity_range<=spikes_min(roi_num),1,'last');
                    % change color of image
                    %filled_ROI = imfill(ROI==roi_num,'holes');
                    filled_ROI = ROI==roi_num;
                    im2plot = im2plot + cat(3,...
                        filled_ROI*cm(ind_range,1),...
                        filled_ROI*cm(ind_range,2),...
                        filled_ROI*cm(ind_range,3));
                end
                
                % make figure
                scale = 1000/size(data_img,2);
                close(findobj('type','figure','name','Events per minute'))
                figure('Name','Events per minute','Color','white','numbertitle','off',...
                    'position',[50 50 size(data_img,2)*scale size(data_img,1)*scale]);
                % plot current view
                subplot('Position',[0 0 1 1]);
                image(0.9*data_img(:,:,[2 1 3]));
                axis image off
                hold on
                % plot ROI activity
                h_tmp=image(im2plot);
                set(h_tmp, 'AlphaData', 0.5);
                % Get center of mass for labels
                for roi_num=1:gui_info.data.count
                    if ~strcmp(data_info.ROI_type{roi_num},{'Neuron','Glial cell'})
                        continue
                    end
                    % text(data_info.ROI_x_y(roi_num,1)-3,data_info.ROI_x_y(roi_num,2),num2str(round(100*spikes_min(roi_num))/100),'FontSize',10,'Color','w','Fontweight','bold')
                end
                
            case 'correlation'   % ===== display correlation between signals
                
                % make figure of size given by number of ROIs
                close(findobj('type','figure','name','Signal correlation'))
                figure('Name','Signal correlation','Color','white','numbertitle','off',...
                    'position',[50 50 (gui_info.data.count)*20 gui_info.data.count*20]);
                
                % make color map
                type = 'default';
                num_bins = 12;
                cm = getColorMap(type,num_bins);
                
                % divide range in parts of colormap
                correlation_range = linspace(gui_info.analysis.correlation_minimum,gui_info.analysis.correlation_maximum,size(cm,1));
                corr_img = 0.7*ones(max(data_info.ROI_label),max(data_info.ROI_label),3);
                for roi_num=1:gui_info.data.count
                    for jj = 1:gui_info.data.count
                        if roi_num~=jj
                            ind = min([find(analysis_info.correlation(roi_num,jj)<=correlation_range,1,'first') length(correlation_range)]);
                            corr_img(data_info.ROI_label(roi_num),data_info.ROI_label(jj),:) = cm(ind,:);
                        end
                    end
                end
                
                subplot('Position',[0.05 0.05 0.9 0.9]);
                image(corr_img);
                axis equal
                
                set(gca,'Xtick',[1:max(data_info.ROI_label)]-0.5)
                set(gca,'Xticklabel',[])
                set(gca,'Ytick',[1:max(data_info.ROI_label)]-0.5)
                set(gca,'Yticklabel',[])
                grid on
                set(gca,'GridLineStyle','-')
                
            case 'raster'    % ===== display raster plot + current dFF

                % make figure of size given by number of ROIs
                figure('Name','Raster Plot','Color','white','numbertitle','off',...
                    'position',[50 50 1000 20*max(data_info.ROI_label)]);
                if gui_info.analysis.add_dFF
                    mark_red = gui_info.data.selected;
                else
                    mark_red = NaN;
                end
                subplot('Position',[0.075 0.1 0.9 0.875])
                
                time_axis = (1:data_info.num_stacks)/data_info.frame_rate/60;
                xlim([0 max(time_axis)])
                hold on
                plot_stim([-3 max(data_info.ROI_label)+3])
                
                % guide lines every ROI and every 50 seconds
                plot(repmat([0;time_axis(end)],[1 max(data_info.ROI_label)-1]),repmat([2:max(data_info.ROI_label)]-0.5,[2 1]),'--','LineWidth',0.5,'Color',[1 1 1]*0.9)
                plot(repmat(0:time_axis(end),[2 1]),repmat([0.2;max(data_info.ROI_label)+0.8],[1 length(0:time_axis(end))]),'-','LineWidth',0.5,'Color',[1 1 1]*0.9)
                % get event times and plot
                for roi_num=1:gui_info.data.count
                    if ~isempty(analysis_info.all_events{roi_num})
                        t = time_axis(analysis_info.all_events{roi_num}(:,1))';
                            y = data_info.ROI_label(roi_num)*ones(size(t));
                        if roi_num == mark_red
                            plot([t,t]',[y-0.4 y+0.4]','-r','LineWidth',1.5)
                        else
                            plot([t,t]',[y-0.4 y+0.4]','-k','LineWidth',1.5)
                        end
                    end
                end
                % beautify
                ylabel('ROI')
                xlabel('time [min]')
                xlim(gui_info.graph.xlim)
                set(gca,'XTick',1:round(time_axis(end)))
                ylim([0 max(data_info.ROI_label)+1])
                set(gca,'YTick',1:max(data_info.ROI_label))
                %set(gca,'YTickLabel',{''})
                set(gca,'TickDir','out','TickLength',[0.005 0.005])
                box on
                
            case 'Participation rate'
                
                % Create the uitable
                close(findobj('type','figure','name','% Event Participation'))
                h_result = figure('Name','% Event Participation','Position',[1165 184 442 751]);
                t = uitable(h_result,'Data',analysis_info.part_rate,...
                    'ColumnName',{'Cell Label','Cell centered','Population centered'},...
                    'ColumnWidth',{120});
                ext = get(t,'Extent');pos =get(t,'Position');
                set(t,'Position',[pos(1) pos(2) ext(3) ext(4)]);
                
            case 'Channel ratio'   %  Table with initial and final ratio between first two channels
                
                % add the average of neurons
                ch_ratio = [analysis_info.ch_ratio;mean(analysis_info.ch_ratio,1)];
                ch_ratio = mat2cell(ch_ratio,ones(size(ch_ratio,1),1),[1 1 1]);
                ch_ratio{end,1} = 'Mean :';
                % Create the uitable
                % make figure with table
                close(findobj('type','figure','name','Channel ratio'))
                h_result = figure('Name','Channel ratio','Position',[1165 184 442 751]);
                t = uitable(h_result,'Data',ch_ratio,...
                    'ColumnName',{'Cell Label','Initial','Final'},...
                    'ColumnWidth',{100});
                ext = get(t,'Extent');pos =get(t,'Position');
                set(t,'Position',[pos(1) pos(2) ext(3) ext(4)]);
        end
    end

    function analysis_parameters(hObject,eventdata,figure_type)
        %-- set the parameters for the plots
        
        switch figure_type
            
            case 'activity_parameters'
                
                % close previous window
                h_params = findobj('Tag','h_params');
                close(h_params)
                h_params = figure('Tag','hactpar','Name','Activity parameters','Position',[1165 184 500 100]);
                % make table (se below for callback R_act_params)
                t = uitable(h_params,'Tag','set_act_params','Data',[gui_info.analysis.activity_maximum;gui_info.analysis.activity_bins],'ColumnName',{''},'RowName',{'Maximum','#Bins'},...
                    'ColumnWidth',{100},'ColumnEditable',true,'CellEditCallback',{@analysis_parameters,'set_act_params'});
                ext = get(t,'Extent');pos =get(t,'Position');
                set(t,'Position',[pos(1) pos(2) ext(3) ext(4)]);
                
            case 'raster'
                % re_calculate properties for figures
                if gui_info.analysis.add_dFF
                    gui_info.analysis.add_dFF = false;
                    set(hObject,'Checked','off')
                else
                    gui_info.analysis.add_dFF = true;
                    set(hObject,'Checked','on')
                end
                
            case 'set_act_params'
                switch eventdata.Indices(1)
                    case 1
                        gui_info.analysis.activity_maximum = eventdata.NewData;
                    case 2
                        gui_info.analysis.activity_bins = eventdata.NewData;
                end
                
            case 'activity_colorbar'  % ===== display colorbar of activity
                
                load(data_info.source_file,'data_img');
                scale = 1000/size(data_img,2);
                close(findobj('type','figure','name','activity colorbar'))
                figure('Name','activity colorbar','Color','white','numbertitle','off',...
                    'position',[50 50 size(data_img,2)*scale size(data_img,1)*scale*0.1]);
                
                % make color map from http://geog.uoregon.edu/datagraphics/color_scales.htm
                original_cm = [41 10 216; 38 77 255; 63 160 255; 114 217 255; 170 247 255; 224 255 255; 255 255 191; 255 224 153; 255 173 114; 247 109 94; 216 38 50; 165 0  33]/255;
                cm = interp1([1:12],original_cm,linspace(1,12,gui_info.analysis.activity_bins));
                % divide range in parts of colormap
                activity_range = linspace(0,gui_info.analysis.activity_maximum,gui_info.analysis.activity_bins);
                
                % Add a colorbar
                subplot('Position',[0 0 1 1]);
                color_bar = ones(round(size(data_img,1)*0.1),size(data_img,2),3);
                unit = floor(size(data_img,2)/length(activity_range));
                for cut = 1:length(activity_range)
                    color_bar(:,(cut-1)*unit+[1:unit],1)=cm(cut,1);
                    color_bar(:,(cut-1)*unit+[1:unit],2)=cm(cut,2);
                    color_bar(:,(cut-1)*unit+[1:unit],3)=cm(cut,3);
                end
                h_tmp=image(color_bar);
                set(h_tmp, 'AlphaData', 0.7);
                for roi_num=1:length(activity_range)
                    text((roi_num-1)*unit+unit/2,size(color_bar,1)*2/3,num2str(round(10*activity_range(roi_num))/10),'FontSize',12,'Color','k')
                end
                axis image off
                
            case 'correlation_parameters'
                
                h_params = findobj('Tag','hactpar');
                close(h_params)
                h_params = figure('Tag','hactpar','Name','Correlation parameters','Position',[1165 184 500 100]);
                % make table (se below for callback set_corr_params)
                t = uitable(h_params,'Tag','set_corr_params','Data',[gui_info.analysis.correlation_minimum;gui_info.analysis.correlation_maximum],'ColumnName',{''},'RowName',{'Minimum','Maximum'},...
                    'ColumnWidth',{100},'ColumnEditable',true,'CellEditCallback',{@analysis_parameters,'set_corr_params'});
                ext = get(t,'Extent');pos =get(t,'Position');
                set(t,'Position',[pos(1) pos(2) ext(3) ext(4)]);
                
            case 'set_corr_params'
                switch eventdata.Indices(1)
                    case 1
                        gui_info.analysis.correlation_minimum = eventdata.NewData;
                    case 2
                        gui_info.analysis.correlation_maximum = eventdata.NewData;
                end
                
            case 'correlation_colorbar'  % ===== display colorbar of correlation
                
                close(findobj('type','figure','name','Correlation Colorbar'))
                figure('Name','Correlation Colorbar','Color','white','numbertitle','off',...
                    'position',[50 50 6*20 (gui_info.data.count)*20]);
                
                % make color map from http://geog.uoregon.edu/datagraphics/color_scales.htm
                original_cm = [41 10 216; 38 77 255; 63 160 255; 114 217 255; 170 247 255; 224 255 255; 255 255 191; 255 224 153; 255 173 114; 247 109 94; 216 38 50; 165 0  33]/255;
                cm = interp1([1:12],original_cm,linspace(1,12,25));
                
                subplot('Position',[0.05 0.05 0.5 0.9]);
                color_bar = ones(size(cm,1),3,3);
                for roi_num = 1:size(cm,1)
                    color_bar(roi_num,:,1)=cm(roi_num,1);
                    color_bar(roi_num,:,2)=cm(roi_num,2);
                    color_bar(roi_num,:,3)=cm(roi_num,3);
                end
                image(color_bar);
                axis xy
                set(gca,'Xtick',[])
                set(gca,'Ytick',[])
                text(4,0.75,num2str(round(10*gui_info.analysis.correlation_minimum)/10),'FontSize',15,'Color','k')
                text(4,size(cm,1)-0.25,num2str(round(10*gui_info.analysis.correlation_maximum)/10),'FontSize',15,'Color','k')
                
            case 'participation'
                
                h_params = findobj('Tag','hactpar');
                close(h_params)
                h_params = figure('Tag','hactpar','Name','Participation rate parameters','Position',[1165 184 500 100]);
                % make table (se below for callback set_part_params)
                t = uitable(h_params,'Tag','set_corr_params','Data',gui_info.analysis.participation_bin_sec,'ColumnName',{''},'RowName',{'Bining in seconds'},...
                    'ColumnWidth',{100},'ColumnEditable',true,'CellEditCallback',{@analysis_parameters,'set_part_params'});
                ext = get(t,'Extent');pos =get(t,'Position');
                set(t,'Position',[pos(1) pos(2) ext(3) ext(4)]);
                
            case 'set_part_params'
                gui_info.analysis.participation_bin_sec = eventdata.NewData;
                % delete stored value to calculate again
                analysis_info.participation_rate = [];
        end
        
        
    end
end