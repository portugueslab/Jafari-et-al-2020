function table_data = calculate_synchronization(data,correlation_cutoff)
% Check the synchronization of the recording by analyzing the correlation
% between the ROIs.

% parameters
if nargin == 1
    correlation_cutoff = 0.25;
end

% calculations
[processed_signal,roi_labels] = extract_data(data);
table_data = do_table(processed_signal,roi_labels);
hierarchical_clusters = cluster_data(processed_signal);
separated_clusters = separate_clusters(hierarchical_clusters, 1-correlation_cutoff);
do_plot(hierarchical_clusters,separated_clusters,1-correlation_cutoff,roi_labels)

end

%%

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

function table_data = do_table(processed_signal,roi_labels)
% calculate pairwise correlations
correlation_matrix = corrcoef(processed_signal);
correlation_matrix(eye(size(correlation_matrix))==1) = 0;

average_correlation = (sum(correlation_matrix,2))/(length(correlation_matrix)-1);
max_correlation = max(correlation_matrix,[],2);

correlation_matrix(eye(size(correlation_matrix))==1) = NaN;

median_correlation = zeros(length(correlation_matrix),1);
for ii = 1:length(correlation_matrix)
    tmp = correlation_matrix(ii,:);
    tmp(isnan(tmp)) = [];
    median_correlation(ii) = median(tmp);
end

% make a table data
table_data = NaN(max(roi_labels),4);
table_data(:,1) = 1:max(roi_labels);
table_data(roi_labels,2) = average_correlation;
table_data(roi_labels,3) = max_correlation;
table_data(roi_labels,4) = median_correlation;

open table_data
% display
h_result = findobj('type','figure','name','Correlations');
if isempty(h_result)
    h_result = figure('Name','Correlations','Color','w');
else
    figure(h_result)
    clf
end
t = uitable(h_result,'Data',table_data,...
    'ColumnName',{'ID','Average','Max','Median'},...
    'ColumnWidth',{120});
ext = get(t,'Extent');pos =get(t,'Position');
set(t,'Position',[pos(1) pos(2) ext(3) ext(4)]);

end


function hierarchical_clusters = cluster_data(processed_signal)
% calculate pairwise correlations
correlation_matrix = corrcoef(processed_signal);
% convert correlation into dissimilarity
dissimilarity = 1 - correlation_matrix;
% calculate hierarchical structure
hierarchical_clusters = linkage(squareform(dissimilarity),'complete');
end


function separated_clusters = separate_clusters(hierarchical_clusters,cutoff)

separated_clusters = cell(0,1);
joint_cluster_ind = size(hierarchical_clusters,1)+1;
for ii = 1:7%size(hierarchical_clusters,1)
    % break if cutoff reached
    if hierarchical_clusters(ii,3) > cutoff
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

% remove joint cluster indices
for ii = 1:length(separated_clusters)
    separated_clusters{ii}(separated_clusters{ii}>size(hierarchical_clusters,1)+1) = [];
    separated_clusters{ii} = sort(separated_clusters{ii});
end

end


function do_plot(hierarchical_clusters,separated_clusters,dissimilarity_cutoff,roi_labels)

h_result = findobj('type','figure','name','Dendrogram');
if isempty(h_result)
    figure('Name','Dendrogram','Color','w');
else
    figure(h_result)
    clf
end
dendrogram(hierarchical_clusters,0,'ColorThreshold',dissimilarity_cutoff);

% labels
xlabel('Cell ID')
ylabel('Pairwise correlation')

title(sprintf('%d/%d cells in %d clusters (threshold corr > %0.2f)', length(cat(2,separated_clusters{:})), size(hierarchical_clusters,1)+1, length(separated_clusters), 1-dissimilarity_cutoff))

% change labels
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

end
