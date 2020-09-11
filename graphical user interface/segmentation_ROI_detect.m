function ROI_coordinates=segmentation_ROI_detect(type,img,cell_center,cell_radius,signal_drop)
% ORIGINAL VERSION BY TSAI-WEN CHEN
%   Created 01/20/2011, by Tsai-Wen Chen
%   Modified 01/21/2011, by Tsai-Wen Chen
%   Adapted 01/12/2014, by Chepe

%% reshape image to have same resolution in X and Y

expand_factor = [1 1];
size_original = size(img);
if abs(expand_factor(1)-expand_factor(2))>0
    % reshape expand_factor to have the minimum re-size equal to 1
    expand_factor = expand_factor/min(expand_factor);
    % resize image and cell center to have equal x-y proportions
    img=double(imresize(img,size_original.*expand_factor));
    cell_center = cell_center.*expand_factor([2 1]);
end

%% calculate radial profile of image around cell center

search_radius=1.3*cell_radius;
angle_samples=90;
profile_samples=50;

img_profile = zeros(profile_samples,angle_samples,3); %[signal,pix_x,pix_y]
for theta=1:angle_samples
    
    % get signal of corresponding pixels values and positions
    [cx,cy,max_profile] = improfile(img,...
        cell_center(1)+[0,search_radius*cos(theta*2*pi/angle_samples)],...
        cell_center(2)+[0,search_radius*sin(theta*2*pi/angle_samples)],...
        profile_samples,'bilinear');
    img_profile(:,theta,1) = max_profile;
    img_profile(:,theta,2:3)=[cx,cy];
end

%% get ROI
switch type
    
    case 'R_drop'
        
        % detect radial threshold where signal falls by 'signal_drop'
        signal_center = max(max(img_profile(1:5,:,1)));
        signal_outside = min(img_profile(:,:,1));
        threshold = signal_drop*(signal_center-signal_outside) + signal_outside;
        
        % convert profile to pixel count from threshold
        [~,threshold_cross]=min(abs(bsxfun(@minus,img_profile(:,:,1),threshold)),[],1);
        drop_profile=-1*abs(bsxfun(@minus,meshgrid(1:profile_samples,1:angle_samples)',threshold_cross))+profile_samples;
        
        % optimize path
        cell_contour = round(find_path(drop_profile));
        
        
    case 'R_peak'
        
         % limit profile around first maximum
        search_range=cell_radius*0.3;
        max_profile = img_profile(:,:,1);
        for ii=1:angle_samples
            % find 50% increase threshold
            [globalmax,globalind]=max(img_profile(:,ii,1));
            threshold=(globalmax-mean(img_profile(1:3,ii,1)))/2+mean(img_profile(1:3,ii,1));
            % find first signal peak above the threshold
            rpeak=find((img_profile(2:end-1,ii,1)>img_profile(1:end-2,ii,1))&(img_profile(2:end-1,ii,1)>img_profile(3:end,ii,1))&(img_profile(2:end-1,ii,1)>threshold),1)+1;
            if isempty(rpeak) % in case of platteau
                rpeak=globalind;
            end
            % find first threshold crossing before peak
            search_min=find(img_profile(1:rpeak,ii,1)>threshold,1);
            if isempty(search_min)
                search_min=ceil(rpeak/2);
            end
            % remove pixels outside
            max_profile(1:search_min,ii)=0;
            max_profile(search_min+round(search_range*profile_samples/search_radius):end,ii)=0;
        end
       
        % optimize path
        cell_contour = round(find_path(max_profile));
        
end

%% Convert profile contour to ROI contour image coordinates

% get corresponding pixels from img_profile(:,:,[2 3]) and reshape to original size
ROI_coordinates = [img_profile(sub2ind(size(img_profile),cell_contour,1:angle_samples,2*ones(1,angle_samples)))',...
    img_profile(sub2ind(size(img_profile),cell_contour,1:angle_samples,3*ones(1,angle_samples)))'];
ROI_coordinates = bsxfun(@rdivide,ROI_coordinates,expand_factor([2 1]));
    
end

% -- helper function to smooth path
function max_path=find_path(reward)
% Use "closeness to threshold" as reward function and find path that
% maximizes it. Nearby pixels are used to find local maximum to smooth the
% path
% start with empty maximum path
max_path=zeros(1,size(reward,2));
% ADD moving maximum twice to find optimal path
for i=2:size(reward,2)
    for j=2:(size(reward,1)-1)
        reward(j,i)=reward(j,i) + max(reward((j-1):(j+1),i-1));
    end
end
reward(:,1)=reward(:,end);
for i=2:size(reward,2)
    for j=2:(size(reward,1)-1)
        reward(j,i)=reward(j,i) + max(reward((j-1):(j+1),i-1));
    end
end
% go through maximum path by starting at a point and moving to the maximum
% of the 3 nearest neighbors
[~,current_point]=max(reward(:,end));
max_path(end) = current_point;
for sample=size(reward,2)-1:-1:1
    [~,next_max]=max(reward(current_point-1:current_point+1,sample));
    current_point = current_point - 2 + next_max;
    max_path(sample)=current_point;
end
end
