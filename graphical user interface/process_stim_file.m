function signal = signal_stim_file(filePath,num_frames,frame_rate)
% Function that opens a stimulation protocol file and converts it in the
% necessesary output format:
% signal = array(N,3) 
%     rows: stimulation set
%    1.col: start time in seconds
%    2.col: end time in seconds
%    3.col: stimulus ID (e.g. orientation)
%%

total_time = num_frames/frame_rate;

% fake data with 5 repetitions
reps = 5;
tmp = linspace(0,total_time,3*reps+1);
signal = [tmp(2:3:3*reps)',tmp(3:3:3*reps)',ones(reps,1)];


end