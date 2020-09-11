function cm = getColorMap(type,num_bins)

if nargin < 2
    num_bins = 12;
end
if nargin < 1
    type = 'bluewhitered' ;
end

try
    eval(['cm = ',type,'(num_bins);'])
catch
    error('Colormap not defined')
end

if size(cm,1) ~= num_bins
    cm = interp1(1:size(cm,1),cm,linspace(1,size(cm,1),num_bins));
end

end

function cm = default(~)
% color map from http://geog.uoregon.edu/datagraphics/color_scales.htm
cm = ...
    [41   10  216;
    38   77  255;
    63  160  255;
    114  217  255;
    170  247  255;
    224  255  255;
    255  255  191;
    255  224  153;
    255  173  114;
    247  109  94;
    216   38  50;
    165    0  33]/255;
end

function cm = bluered(~)
% color map from http://colorbrewer.org.
cm =...
    [103	0       31;
    178     24      43;
    214     96      77;
    244     165     130;
    253     219     199;
    209     229     240;
    146     197     222;
    67	    147     195;
    33      102     172;
    5       48      97]/255;

cm = flipud(cm);

end

function cm = bluewhitered(~)
cm = ...
    [0        0         0.5000;
    0         0.1667    0.6667;
    0         0.3333    0.8333;
    0         0.5000    1.0000;
    0.3333    0.6667    1.0000;
    0.6667    0.8333    1.0000;
    1.0000    1.0000    1.0000;
    1.0000    1.0000    1.0000;
    1.0000    0.6000    0.6000;
    1.0000    0.2000    0.2000;
    0.9000    0         0;
    0.7000    0         0;
    0.5000    0         0];
end

