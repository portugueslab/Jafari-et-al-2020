function ROI_xy = segmentation_ROI_alignment(img,ROI_xy)
% Short GUI to find the parameters that align the ROI centers to the input
% image. The parameters to vary are rotation angle (+ center point),
% scaling and x-y translation. Keyboard instructions are detailed below.
% chepe@nld.ds.mpg.de

%% get center of ROIs to visually get rotation
roi_center_xy = zeros(2,length(ROI_xy));
for roi_ii=1:length(ROI_xy)
    roi_center_xy(:,roi_ii) = mean(ROI_xy{roi_ii},2);
end

% Parameters
theta = 0;
shift_xy = [0 0];
scale = 1;
center_xy = ceil([size(img,2)/2,size(img,1)/2]);
        
%% Make gui
fig_scale = 1000/size(img,2);
hfig = figure('Name','ROI alignment GIU',...%'WindowScrollWheelFcn',@chooseAngle,...
    'Position', [50 50 size(img,2)*fig_scale size(img,1)*fig_scale],'numbertitle','off',...
    'Color','w','menubar','none','numbertitle','off',...
    'WindowButtonDownFcn',@chooseCenter,'KeyPressFcn',@chooseShifts);

uimenu(hfig,'Label','Instructions','Callback',@showInstructions)
uimenu(hfig,'Label','Reset','Callback',@resetTransform)
uimenu(hfig,'Label','Apply','Callback',@applyTransform)
set(hfig,'Pointer','fullcrosshair')

% freeze main program while this one is still running
doPlot()
uiwait(hfig)

%% Instructions

    function showInstructions(src,evnt)
        txt = {'INSTRUCTIONS FOR CHOOSING ALIGNMENT PARAMETERS';...
            '  ';...
            'The goal is to find the translation, rotation and scaling to align';...
            'the ROI centers (magenta stars) to the background image';...
            '  ';...
            'Commands:';...
            '- Left mouse click to deterime rotation center point';...
            '- Keyboard "a" and "s" to vary rotation angle';...
            '- Keyboard "+" and "-" to determine scaling'
            '- Keyboard "arrows" to determine x-y shift';...
            '  (applied AFTER rotation+scaling)';...
            '- Click on "Apply" to align ROIs and return to main program';...
            '- Clicl on "Reset" to restart parameters'};
        msgbox(txt,'Instructions')
    end

%% Plot
    function doPlot
        
        % --- apply current transform:
        
        % Shift to center of rotation
        roi_center_xy_T = bsxfun(@minus,roi_center_xy,center_xy');
        % Rotate
        T=scale*[cos(theta/180*pi) sin(theta/180*pi) ;...
                -sin(theta/180*pi) cos(theta/180*pi)];
        roi_center_xy_T = T * roi_center_xy_T;
        % Shift back to original position
        roi_center_xy_T = bsxfun(@plus,roi_center_xy_T,center_xy');
        % Shift to new position
        roi_center_xy_T = bsxfun(@plus,roi_center_xy_T,shift_xy');
        
        % --- do plot
        
        cla
        subplot('Position',[0 0 1 1]);
        image(img)
        axis image off
        hold on
        plot(roi_center_xy_T(1,:),roi_center_xy_T(2,:),'*m')
        
        plot(center_xy(1)+scale*10*[0 cos(-theta/180*pi)],center_xy(2)+scale*10*[0 sin(-theta/180*pi)],'g-')
        plot(center_xy(1)+[0 shift_xy(1)],center_xy(2)+[0 shift_xy(2)],'g-o')
        plot(center_xy(1),center_xy(2),'g+')
    end

%%
    function resetTransform(src,evnt)
        theta = 0;
        scale = 1;
        shift_xy = [0 0];
        center_xy = ceil([size(img,2)/2,size(img,1)/2]);
        doPlot()
    end
    function chooseCenter(src,evnt)
        tmp = get(gca,'CurrentPoint');
        theta = 0;
        scale = 1;
        shift_xy = [0 0];
        center_xy = tmp(1,1:2);
        doPlot()
    end
%     function chooseAngle(src,evnt)
%         if evnt.VerticalScrollCount > 0
%             theta = theta - 1;
%         elseif evnt.VerticalScrollCount < 0
%             theta = theta + 1;
%         end
%         doPlot()
%     end
    function chooseShifts(src,evnt)
        switch evnt.Key
            case 'rightarrow'
                shift_xy(1) = shift_xy(1) + 1;
            case 'leftarrow'
                shift_xy(1) = shift_xy(1) - 1;
            case 'uparrow'
                shift_xy(2) = shift_xy(2) - 1;
            case 'downarrow'
                shift_xy(2) = shift_xy(2) + 1;
            case 'add'
                scale = scale + 0.01;
            case 'subtract'
                scale = scale - 0.01;
            case 'a'
                theta = theta - 1;
            case 's'
                theta = theta + 1;
        end
        doPlot()
    end
%%


%% When finished, apply transform
    function applyTransform(src,evnt)
        % Inverse mapping method:
        % S. Shah, J.K. Aggrawal, Pattern Recognition, 29 (11) (1996), pp. 1775–1778
        disp('Transform parameters set:')
        
        % start from rotated coordinages (3x larger frame than original)
        width = size(img,2);
        height = size(img,1);
        
        [ix,iy] = meshgrid(round(-width/2):round(1.5*width),round(-height/2):round(1.5*height));
        xy_rot = [ix(:)';iy(:)'];
        xy_rot = bsxfun(@minus,xy_rot,center_xy');
        
        % Rotate backwards
        T= 1/scale *[cos(theta/180*pi) -sin(theta/180*pi) ;...
                     sin(theta/180*pi) cos(theta/180*pi)];
        
        xy_orig = T * xy_rot;
        
        % shift indices again
        xy_rot = bsxfun(@plus,xy_rot,center_xy');
        xy_orig = bsxfun(@plus,xy_orig,center_xy');
        
        xy_rot = round(xy_rot);
        xy_orig = round(xy_orig);
        
        % apply to ROIs
        for roi_jj=1:length(ROI_xy)
            disp(['- Transforming ROI ',num2str(roi_jj),' of ',num2str(length(ROI_xy))])
            % check original indices that are inside a ROI
            ind = ismember(xy_orig',ROI_xy{roi_jj}','rows');
            ROI_xy{roi_jj} = xy_rot(:,ind);
            % do x-y shifts
            ROI_xy{roi_jj} = bsxfun(@plus,ROI_xy{roi_jj},shift_xy');
        end
        
        % close figure to continue with main program
        close(hfig)
    end
end