function export_dFF_to_excel(roi_list)


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


%% open file

[FileName,PathName] = uigetfile([pwd,filesep,'*_extracted.mat'],'MultiSelect','off','Select file to open');
source_file = [PathName,filesep,FileName];

[~,FileNameNoExtension,~] = fileparts(source_file);
destination_file = [PathName,filesep,FileNameNoExtension,'_signal',file_sufix];

if exist(destination_file,'file')
    delete(destination_file);
end

%% get data and export

load(source_file,'data');

num_roi = length(data);
max_label = max(cat(1,data.label));
num_bins = length(data(1).processed);

t = linspace(0,num_bins/data(1).frame_rate,num_bins)';

data_array = cell(num_bins+1,max_label+1);
data_array{1} = 'time / ID';
data_array(2:end,1) = num2cell(t);

data_array(1,2:end) = num2cell(1:max_label);

if nargin == 0
   roi_list = 1:num_roi; 
end
for roi_ii = roi_list(:)'
    if ismember(data(roi_ii).type,{'Glial Cell','Neuron'})
        data_array(2:end,data(roi_ii).label+1) = ...
            num2cell(data(roi_ii).processed(1,:)');
    end
end

expfun(destination_file,data_array,'signal');

end
