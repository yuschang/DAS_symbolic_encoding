% ===== User Input (Example) =====
    fmin   = 6.0;   % kHz
    fmax   = 20.0;  % kHz
    fstart = 8.0;   % kHz
    fend   = 16.0;  % kHz
    dur    = 1.2;   % s

    fs_plot = 1500;                  % Sampling frequency for visualization (Hz)
    plateau_ratio = 0.30;            % Ratio of the static plateau at the end of ADS
    dx = 0.03;                       % ★ Resampling interval (seconds) e.g., 0.01 s

    % ===== Encoding Thresholds (Hz / s) =====
    thresh_freq_change_Hz     = 150;   % Frequency change threshold between adjacent samples
    thresh_duration_change_s  = 0.030; % Threshold to remove "S" segments that are too short

    % ===== Generate Whistles 1 and 2 =====
    [t1, w1_kHz] = make_whistle_modulated_equal_minmax(fmin, fmax, fstart, fend, dur, fs_plot, 2);
    [t2, w2_kHz] = make_whistle_ADS(fmin, fmax, fstart, fend, dur, fs_plot, plateau_ratio);

    % ===== (Figure 1, 2) Verify Original Contours: Display Start/End only =====
    figure('Name','Whistle 1: smooth modulated (bounded)','Color','w');
    plot(t1, w1_kHz, 'LineWidth', 2); hold on;
    scatter(0,   fstart, 50, '^', 'filled', 'DisplayName','Start');
    scatter(dur, fend,   50, 'v', 'filled', 'DisplayName','End');
    xlim([-0.2 1.4])
    ylim([0 24])
    xlabel('Time (s)'); ylabel('Frequency (kHz)'); grid on; box on;
    title('Whistle 1 (smooth modulated)'); legend('Contour','Start','End','Location','best');
    grid off

    figure('Name','Whistle 2: ADS (ascent → descent → static)','Color','w');
    plot(t2, w2_kHz, 'LineWidth', 2); hold on;
    scatter(0,   fstart, 50, '^', 'filled', 'DisplayName','Start');
    scatter(dur, fend,   50, 'v', 'filled', 'DisplayName','End');
    xlim([-0.2 1.4])
    ylim([0 24])    
    xlabel('Time (s)'); ylabel('Frequency (kHz)'); grid on; box on;
    title('Whistle 2 (ADS)'); legend('Contour','Start','End','Location','best');
    grid off

    % ===== ★ Resampling → Input to func_getWPE (x: s, y: Hz) =====
    [x1, y1_Hz] = resample_curve(t1, w1_kHz*1e3, dx);  % kHz to Hz
    [x2, y2_Hz] = resample_curve(t2, w2_kHz*1e3, dx);

    [pat1, para1, key1] = func_getWPE(x1, y1_Hz, thresh_freq_change_Hz, thresh_duration_change_s);
    [pat2, para2, key2] = func_getWPE(x2, y2_Hz, thresh_freq_change_Hz, thresh_duration_change_s);

    % ===== (Figure 3) Whistle 1 ADS Encoding Visualization (Segments after Resampling) =====
    figure('Name','Figure 3: Whistle 1 ADS encoding (resampled)','Color','w');
    hold on; grid on; box on;
    plot(x1, y1_Hz*1e-3, 'Color', [0.8 0.8 0.8], 'LineWidth', 1.5, 'DisplayName','Contour'); % Hz to kHz
    cmap = struct('A',[0.85 0.2 0.2], 'D',[0.2 0.2 0.85], 'S',[0.2 0.6 0.2]);
    for i = 1:numel(key1.start)
        seg_t = x1(key1.start(i):key1.stop(i));
        seg_f = y1_Hz(key1.start(i):key1.stop(i))*1e-3; % kHz
        plot(seg_t, seg_f, 'LineWidth', 3, 'Color', cmap.(pat1(i)));
    end
    scatter(0,   fstart, 50, '^', 'filled', 'k', 'DisplayName','Start');
    scatter(dur, fend,   50, 'v', 'filled', 'k', 'DisplayName','End');
    xlim([-0.2 1.4])
    ylim([0 24])
    xlabel('Time (s)'); ylabel('Frequency (kHz)');
    title(sprintf('Figure 3: Whistle 1 ADS (pattern = %s)', pat1));
    % legend('Location','best');
    grid off

    % ===== (Figure 4) Whistle 2 ADS Encoding Visualization (Segments after Resampling) =====
    figure('Name','Figure 4: Whistle 2 ADS encoding (resampled)','Color','w');
    hold on; grid on; box on;
    plot(x2, y2_Hz*1e-3, 'Color', [0.8 0.8 0.8], 'LineWidth', 1.5, 'DisplayName','Contour');
    for i = 1:numel(key2.start)
        seg_t = x2(key2.start(i):key2.stop(i));
        seg_f = y2_Hz(key2.start(i):key2.stop(i))*1e-3; % kHz
        plot(seg_t, seg_f, 'LineWidth', 3, 'Color', cmap.(pat2(i)));
    end
    scatter(0,   fstart, 50, '^', 'filled', 'k', 'DisplayName','Start');
    scatter(dur, fend,   50, 'v', 'filled', 'k', 'DisplayName','End');
    xlim([-0.2 1.4])
    ylim([0 24])
    xlabel('Time (s)'); ylabel('Frequency (kHz)');
    title(sprintf('Figure 4: Whistle 2 ADS (pattern = %s)', pat2));
    % legend('Location','best');
    grid off

    % ===== Console Report (Maintain previous requirements) =====
    report_ads_summary('Whistle 1', pat1, para1);
    report_ads_summary('Whistle 2', pat2, para2);


% ---------- Resampling Utility (Interpolation at equal dx intervals; x: s, y: Hz) ----------
function [new_x, new_y] = resample_curve(x, y, dx)
    % Sort and remove duplicates (for interpolation stability)
    [ux, ia] = unique(x(:).', 'stable');  % Maintain row vector
    uy = y(ia);
    % Create evenly spaced x
    new_x = ux(1):dx:ux(end);
    % Linear interpolation (do not generate values outside the original range)
    new_y = interp1(ux, uy, new_x, 'linear');
end

% ===== Generator/Report Utilities =====

function [t, w_kHz] = make_whistle_modulated_equal_minmax(fmin, fmax, fstart, fend, dur, fs_plot, harmonicFactor)
    if fmin > fmax, tmp=fmin; fmin=fmax; fmax=tmp; end
    fstart = min(max(fstart, fmin), fmax);
    fend   = min(max(fend,   fmin), fmax);
    if dur <= 0, error('Duration must be positive.'); end
    minSE = min(fstart, fend); maxSE = max(fstart, fend);

    N = max(300, round(fs_plot*dur));
    t = linspace(0, dur, N);
    tau = t/dur;

    b = fstart + (fend - fstart)*tau;
    phi = sin(harmonicFactor*pi*2 * tau);

    phi_pos = max(phi, 0);  phi_neg = max(-phi, 0);
    idx_pos = phi_pos > 0;  idx_neg = phi_neg > 0;

    if any(idx_pos), a_pos_max = min((maxSE - b(idx_pos)) ./ phi_pos(idx_pos)); else, a_pos_max = 0; end
    if any(idx_neg), a_neg_max = min((b(idx_neg) - minSE) ./ phi_neg(idx_neg)); else, a_neg_max = 0; end
    a_pos = 0.999 * a_pos_max;  a_neg = 0.999 * a_neg_max;

    w_kHz = b + a_pos*phi_pos - a_neg*phi_neg;
end

function [t, w_kHz] = make_whistle_ADS(fmin, fmax, fstart, fend, dur, fs_plot, plateau_ratio)
    if fmin > fmax, tmp=fmin; fmin=fmax; fmax=tmp; end
    fstart = min(max(fstart, fmin), fmax);
    fend   = min(max(fend,   fmin), fmax);
    if dur <= 0, error('Duration must be positive.'); end
    plateau_ratio = min(max(plateau_ratio, 0.05), 0.8);

    N = max(300, round(fs_plot*dur));
    t = linspace(0, dur, N);

    p = plateau_ratio * dur;
    r = dur - p; if r <= 0, error('plateau_ratio is too large (r > 0 required)'); end

    m = 0.5*(fstart + fend);
    c = m + (p/r)*(fstart - fend);
    c = min(max(c, min([fmin, fmax, fstart, fend])), max([fmin, fmax, fstart, fend]));

    r1 = r/2; r2 = r/2;
    w_kHz = zeros(size(t));

    idx1 = (t >= 0) & (t <= r1);
    if any(idx1)
        tau1 = (t(idx1))/r1;
        z1   = 0.5*(1 - cos(pi * tau1));
        w_kHz(idx1) = fstart + (c - fstart) * z1;
    end

    idx2 = (t > r1) & (t <= r1 + r2);
    if any(idx2)
        tau2 = (t(idx2) - r1)/r2;
        z2   = 0.5*(1 - cos(pi * tau2));
        w_kHz(idx2) = c + (fend - c) * z2;
    end

    idx3 = (t > r1 + r2) & (t <= dur);
    w_kHz(idx3) = fend;

    w_kHz = min(max(w_kHz, fmin), fmax);
end

function report_ads_summary(tag, pattern, para)
    Aidx = find(pattern=='A'); Didx = find(pattern=='D'); Sidx = find(pattern=='S');
    if ~isempty(Aidx)
        start_f_A = para(Aidx(1), 2); end_f_A = para(Aidx(end),4); dur_A = sum(para(Aidx,5));
    else
        start_f_A = NaN; end_f_A = NaN; dur_A = 0;
    end
    if ~isempty(Didx)
        start_f_D = para(Didx(1), 2); end_f_D = para(Didx(end),4); dur_D = sum(para(Didx,5));
    else
        start_f_D = NaN; end_f_D = NaN; dur_D = 0;
    end
    dur_S = sum(para(Sidx,5));
    fprintf('--- %s ADS summary (resampled) ---\n', tag);
    fprintf('start f A (kHz): %.3f | end f A (kHz): %.3f | dur A (s): %.3f\n', start_f_A, end_f_A, dur_A);
    fprintf('start f D (kHz): %.3f | end f D (kHz): %.3f | dur D (s): %.3f\n', start_f_D, end_f_D, dur_D);
    fprintf('dur S (s): %.3f\n', dur_S);
end