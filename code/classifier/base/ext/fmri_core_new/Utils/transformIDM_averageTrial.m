% transformIDM_averageTrial(info,data,meta) 
%
% transformIDM_averageTrial(info,data,meta) 
% Returns a copy of info,data,meta in which each trial is the average of the 
% snapshots in this trial.
%
% Example: 
% [info2,data2,meta2] = transformIDM_averageTrial(info,data,meta) 
%
% History
% - 6/13/05 Wei Wang Created file.

function [rinfo,rdata,rmeta] = transformIDM_averageTrial(info,data,meta)
    rinfo = info;
    rmeta = meta;
    rdata = data;
    for j=1:length(info)
        rinfo(j).len = 1;
        rdata{j} = mean(data{j},1);
    end
    rmeta.nsnapshots = length(info);