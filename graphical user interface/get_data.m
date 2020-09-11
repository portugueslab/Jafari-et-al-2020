% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%         CODE TO CALCULAT AVERAGE ROI AND NEUROPILE SIGNAL
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load file
if ~exist('PathName','var')
    PathName = pwd;
end
[FileName,PathName] = uigetfile([PathName,filesep,'*_extracted.mat'],'MultiSelect','off','Select file to open');
load([PathName,filesep,FileName])

% get data
data_array = NaN(length(data),4);
for roi_ii = 1:length(data)
       
    data_array(roi_ii,1) = data(roi_ii).label;
    % add data only of neurons
    if ismember(data(roi_ii).type,{'Glial Cell','Neuron'})
    data_array(roi_ii,2) = mean(data(roi_ii).signal(1,:));
    data_array(roi_ii,3) = mean(data(roi_ii).neuropile(1,:));
    data_array(roi_ii,4) = data_array(roi_ii,2)/data_array(roi_ii,3);
    end
end

% make table
close all
h_fig = figure('Name',FileName);
t = uitable(h_fig,'Data',data_array,'ColumnName',{'ID','Signal','Neuropile','S/N'});
ext = get(t,'Extent');pos =get(t,'Position');
set(t,'Position',[pos(1) pos(2) ext(3) ext(4)]);