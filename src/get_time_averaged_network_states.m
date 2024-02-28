function [states,partition] = get_time_averaged_network_states(ts,roi2net,activity_thr)
%% get time-averaged network states eg Ceballos 2024 https://www.biorxiv.org/content/10.1101/2024.01.24.577068v1.full
% ts = time by roi matrix - NOTE - regional time-series (ie each column)
% should be zscored prior to input to be consistent with reference.
% roi2net = roi by num_nets matrix where each column is a binary vector
% with 1's indicating regions belonging to that network
% activity_thr = activity threshold by which the max network must reach in
% order for that timepoint to recieve an assignment. If the threshol is not
% met, the timepoint will be assigned NaN. Default = 0.5 (based on ref).

[rois,num_nets] = size(roi2net);
states = NaN(rois,num_nets);
[trs,rois2] = size(ts);
if rois~=rois2
    ts = transpose(ts);
    [trs,rois2] = size(ts);
    if rois~=rois2
        error('dimensions of ts and roi2nets not consistent')
    end
end

% For each TR, calculate the network with the most (positively) active ROIs
time_labels = zeros(trs, 1);
for tr = 1:trs
    ts_tr = transpose(ts(tr,:));
    ts_tr = repmat(ts_tr,1,num_nets);
    
    % Extract values for the given network
    % Mean across all ROIs in network with positive activity
    ts_tr_net = ts_tr.*roi2net;
    ts_tr_net(ts_tr_net<=0)=NaN;
    net_means = nanmean(ts_tr_net);
    
    
    [~, idx] = max(net_means);
    if net_means(idx) >= activity_thr
        time_labels(tr) = idx;
    else
        time_labels(tr) = NaN;
    end
end

% Calculate mean activity across dominant time points for each state
for i = 1:num_nets
    states(:,i) = mean(ts(time_labels==i, :), 1);
end

partition=time_labels;

end