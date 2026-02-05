function [WPE_pattern, WPE_para, WPE_key_idx] = func_getWPE(t_axis, freq_profile, thresh_freq, thresh_dur)
% FUNC_GETWPE Whistle Pattern Encoding
% Encodes a frequency profile into Ascent (A), Descent (D), and Static (S) segments.

    num_samples = length(freq_profile);
    if num_samples < 2
        WPE_pattern = ''; WPE_para = []; WPE_key_idx = struct('start',[],'stop',[]);
        return;
    end

    % 1. Determine local patterns based on frequency gradient
    gaps = diff(freq_profile);
    local_patterns = repmat('S', 1, num_samples - 1);
    local_patterns(gaps >  thresh_freq) = 'A';
    local_patterns(gaps < -thresh_freq) = 'D';

    % 2. Identify segment boundaries (where pattern changes)
    % Find indices where the pattern differs from the previous sample
    change_idx = [1, find(local_patterns(1:end-1) ~= local_patterns(2:end)) + 1];
    
    WPE_pattern = local_patterns(change_idx);
    WPE_idx_start = change_idx';
    WPE_idx_stop  = [change_idx(2:end) - 1, num_samples]';

    % 3. Calculate parameters for each segment
    % Format: [start_t, start_f_kHz, stop_t, stop_f_kHz, duration, delta_f_kHz]
    num_segs = length(WPE_pattern);
    WPE_para = zeros(num_segs, 6);
    
    for n = 1:num_segs
        i1 = WPE_idx_start(n);
        i2 = WPE_idx_stop(n);
        
        dur = t_axis(i2) - t_axis(i1);
        df  = (freq_profile(i2) - freq_profile(i1)) * 1e-3;
        
        WPE_para(n, :) = [t_axis(i1), freq_profile(i1)*1e-3, ...
                          t_axis(i2), freq_profile(i2)*1e-3, ...
                          dur, df];
    end

    % 4. Filter out short 'S' segments (unless they bridge identical patterns)
    % e.g., Keep 'S' in "ASA", but remove in "ASD"
    idx_S = find(WPE_pattern == 'S');
    is_too_short = WPE_para(:, 5) <= thresh_dur;
    
    to_remove = [];
    for i = 1:length(idx_S)
        curr_idx = idx_S(i);
        
        if is_too_short(curr_idx)
            % Check neighbors: If it's a bridge between same patterns, keep it
            if curr_idx > 1 && curr_idx < length(WPE_pattern)
                if WPE_pattern(curr_idx-1) == WPE_pattern(curr_idx+1)
                    continue; % Skip removal (exception case)
                end
            end
            to_remove = [to_remove, curr_idx]; %#ok<AGROW>
        end
    end

    % Apply removal
    WPE_pattern(to_remove)    = [];
    WPE_para(to_remove, :)    = [];
    WPE_idx_start(to_remove)  = [];
    WPE_idx_stop(to_remove)   = [];

    % 5. Format Output
    WPE_key_idx.start = WPE_idx_start;
    WPE_key_idx.stop  = WPE_idx_stop;
end