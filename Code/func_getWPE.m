function [WPE_pattern, WPE_para, WPE_key_idx ] = func_getWPE(currentProfile_timeaxis, currentProfile, thresh_freq_change, thresh_duration_change )

% Whistle Pattern Encoding function
% returns encoded whistle pattern into ADS_ symbols based on defined threshold

    
WPE_code_Fs = currentProfile(1); % start frequency
WPE_pattern = '';
WPE_idx = [];

changeGrad = zeros(1, length(currentProfile));
prevPattern = '';
for n = 1: length(currentProfile)-1
    gap = currentProfile(n+1) - currentProfile(n);

    if gap > thresh_freq_change
        changeGrad(n+1) = 1;
        pattern_now = 'A';
    elseif abs(gap)<= thresh_freq_change
        changeGrad(n+1) = 0;
        pattern_now = 'S';
    elseif gap < -thresh_freq_change
        changeGrad(n+1) = -1;
        pattern_now = 'D';
    end
    % check the pattern
    if ~strcmp(pattern_now,prevPattern)      
        WPE_pattern = [WPE_pattern pattern_now];
        prevPattern = pattern_now;
        WPE_idx = [WPE_idx n];
    end
end

% the start and stop point of each segment
WPE_idx_start = [];
WPE_idx_stop = [];

if length(WPE_idx) ==1 
    WPE_idx_start = [1];
    WPE_idx_stop = [length(currentProfile_timeaxis)];
else
    for n = 1:length(WPE_idx)
        if n ==1
            WPE_idx_start = [ 1 ];
            WPE_idx_stop = [WPE_idx(n+1)-1];
        elseif n == length(WPE_idx)
            WPE_idx_start = [WPE_idx_start; WPE_idx(n)];
            WPE_idx_stop = [WPE_idx_stop; length(currentProfile)];
        else
            WPE_idx_start = [WPE_idx_start; WPE_idx(n)];
            WPE_idx_stop = [WPE_idx_stop; WPE_idx(n+1)-1];
        end
    end
end

WPE_para = zeros(length(WPE_pattern), 6);
% syntax [startP_x startP_y stopP_x stopP_y duration f_change]
for n = 1: length(WPE_idx_start)
    WPE_para(n, :) = [currentProfile_timeaxis(WPE_idx_start(n))  currentProfile(WPE_idx_start(n))*1e-3...
                        currentProfile_timeaxis(WPE_idx_stop(n))  currentProfile(WPE_idx_stop(n))*1e-3...
                        currentProfile_timeaxis(WPE_idx_stop(n))-currentProfile_timeaxis(WPE_idx_start(n)) ...
                        (currentProfile(WPE_idx_stop(n))-currentProfile(WPE_idx_start(n)))*1e-3 ];
end

    % remove the S class meet the condition
    idx_low_change = find([WPE_para(:,5)]<=thresh_duration_change);
    idx_S_class = find(WPE_pattern =='S');

    if ~isempty(idx_low_change) && ~isempty(idx_S_class)
        % find the intersect
        idx_remove = intersect(idx_low_change, idx_S_class);   
        idx_exception = [];
        % remain the case when both side of the S patter has the same
        % pattern like ASA, DSD to avoid a reulst of AA or DD
        for iii = 1: length(idx_remove)
            if idx_remove(iii)>1 && idx_remove(iii)< length(WPE_pattern)
                prior_pattern = WPE_pattern(idx_remove(iii)-1);
                following_pattern = WPE_pattern(idx_remove(iii)+1);
                if strcmp(prior_pattern, following_pattern)
                    idx_exception = [idx_exception iii];
                end
            end
        end

        if ~isempty(idx_exception)
            idx_remove(idx_exception) = [];
        end

        WPE_para(idx_remove,: )=[];
        WPE_idx_start(idx_remove) = [];
        WPE_idx_stop(idx_remove) = [];
        WPE_pattern(idx_remove)=[];
    end


WPE_key_idx.start = WPE_idx_start;
WPE_key_idx.stop = WPE_idx_stop;