function  varargout = do_cluster_analysis(data,correlation_cutoff)
% Check the synchronization of the recording by analyzing the correlation
% between the ROIs.

%% get dissimilarity cutoff
if nargin == 1
    correlation_cutoff = 0.25;
end
dissimilarity_cutoff = 1 - correlation_cutoff;

%% process data
[processed_signal,roi_labels] = extract_data(data);

%% calculate hierarchical structure from correlation
correlation_matrix = corrcoef(processed_signal);
dissimilarity = 1 - correlation_matrix;
hierarchical_clusters = linkage(squareform(dissimilarity),'complete');

%% separate valid clusters and get ID
separated_clusters = cell(0,1);
joint_cluster_ind = size(hierarchical_clusters,1)+1;
for ii = 1:size(hierarchical_clusters,1)
    % break if cutoff reached
    if hierarchical_clusters(ii,3) > dissimilarity_cutoff
        break
    end
    joint_cluster_ind = joint_cluster_ind + 1;
    
    % find corresponding cluster
    corresponding_cluster = find(cellfun(@(x)any(ismember(x,hierarchical_clusters(ii,1:2))),separated_clusters));
    if isempty(corresponding_cluster)
        % make new
        separated_clusters{end+1} = [hierarchical_clusters(ii,1:2) joint_cluster_ind];
    elseif length(corresponding_cluster)>1
        % join two clusters
        separated_clusters{corresponding_cluster(1)} = [separated_clusters{corresponding_cluster(1)} separated_clusters{corresponding_cluster(2)} joint_cluster_ind];
        separated_clusters(corresponding_cluster(2)) = [];
    else
        % append to previous
        separated_clusters{corresponding_cluster} = [separated_clusters{corresponding_cluster} hierarchical_clusters(ii,1:2) joint_cluster_ind];
    end
end
% remove joint cluster indices and convert to real labels
for ii = 1:length(separated_clusters)
    separated_clusters{ii}(separated_clusters{ii}>size(hierarchical_clusters,1)+1) = [];
    separated_clusters{ii} = sort(roi_labels(separated_clusters{ii}))';
end
% sort by size
[~,idx] = sort(cellfun(@length,separated_clusters),'descend');
separated_clusters = separated_clusters(idx);

%% plot dendrogram
h_fig = findobj('type','figure','name','Dendrogram');
if isempty(h_fig)
    figure('Name','Dendrogram','Color','w');
else
    figure(h_fig)
    clf
end
dendrogram(hierarchical_clusters,0,'ColorThreshold',dissimilarity_cutoff);

% label axes
xlabel('ROI label')
ylabel('Pairwise correlation')

title(sprintf('%d/%d cells in %d clusters (threshold corr > %0.2f)', length(cat(2,separated_clusters{:})), size(hierarchical_clusters,1)+1, length(separated_clusters), 1-dissimilarity_cutoff))

% change labels to match correlation and roi_label
h = gca;

old_tick = get(h,'XTickLabel');
if iscell(old_tick)
    old_tick = cellfun(@(x)str2num(x),old_tick);
else
    old_tick = str2num(old_tick);
end
new_tick = roi_labels(old_tick);
set(h,'XTickLabel',new_tick);

old_tick = get(h,'YTickLabel');if iscell(old_tick)
    old_tick = cellfun(@(x)str2num(x),old_tick);
else
    old_tick = str2num(old_tick);
end
new_tick = 1-old_tick;
set(h,'YTickLabel',new_tick);

%% print report

fprintf('\n=== REPORT:\n\n')
fprintf('Correlation threshold = %0.2f\n',1-dissimilarity_cutoff)
fprintf('%d/%d cells in %d clusters\n\n', length(cat(2,separated_clusters{:})), size(hierarchical_clusters,1)+1, length(separated_clusters))
for ii = 1:length(separated_clusters)
    g=sprintf('%d ', separated_clusters{ii});
    fprintf('Cluster %i (%i members) : %s \n',ii,length(separated_clusters{ii}),g)
end

if nargout==1
    varargout{1} = separated_clusters;    
end

end

%% function to process the data

function [processed_signal,roi_labels] = extract_data(data)
% extract signal from data
processed_signal = zeros(size(data(1).processed,2),0);
roi_labels = zeros(0,1);
for ii = 1:length(data)
    if strcmp(data(ii).type,'Neuron')
        roi_labels = [roi_labels; data(ii).label];
        processed_signal = [processed_signal, data(ii).processed(1,:)'];
    end
end
end
