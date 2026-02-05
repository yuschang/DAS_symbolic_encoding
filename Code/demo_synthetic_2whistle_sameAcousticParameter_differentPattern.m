    % ===== 사용자 입력 (예시) =====
    fmin   = 6.0;   % kHz
    fmax   = 20.0;  % kHz
    fstart = 8.0;   % kHz
    fend   = 16.0;  % kHz
    dur    = 1.2;   % s

    fs_plot = 1500;                  % 시각화용 시간 샘플링(Hz)
    plateau_ratio = 0.30;            % ADS 끝 Static 비율
    dx = 0.03;                       % ★ 리샘플링 간격(초) 예: 0.01 s

    % ===== 인코딩 임계값 (Hz / s) =====
    thresh_freq_change_Hz     = 150;   % 인접 샘플 간 주파수 변화 임계값
    thresh_duration_change_s  = 0.030; % 너무 짧은 S 구간 제거 임계값

    % ===== 휘슬 1, 2 생성 =====
    [t1, w1_kHz] = make_whistle_modulated_equal_minmax(fmin, fmax, fstart, fend, dur, fs_plot, 2);
    [t2, w2_kHz] = make_whistle_ADS(fmin, fmax, fstart, fend, dur, fs_plot, plateau_ratio);

    % ===== (그림 1, 2) 원형 윤곽 확인: start/end만 표시 =====
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

    % ===== ★ 리샘플링 → func_getWPE에 투입 (x: s, y: Hz) =====
    [x1, y1_Hz] = resample_curve(t1, w1_kHz*1e3, dx);  % kHz→Hz
    [x2, y2_Hz] = resample_curve(t2, w2_kHz*1e3, dx);

    [pat1, para1, key1] = func_getWPE(x1, y1_Hz, thresh_freq_change_Hz, thresh_duration_change_s);
    [pat2, para2, key2] = func_getWPE(x2, y2_Hz, thresh_freq_change_Hz, thresh_duration_change_s);

    % ===== (그림 3) 휘슬 1 ADS 인코딩 시각화(리샘플 후 세그먼트 표시) =====
    figure('Name','Figure 3: Whistle 1 ADS encoding (resampled)','Color','w');
    hold on; grid on; box on;
    plot(x1, y1_Hz*1e-3, 'Color', [0.8 0.8 0.8], 'LineWidth', 1.5, 'DisplayName','Contour'); % Hz→kHz
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

    % ===== (그림 4) 휘슬 2 ADS 인코딩 시각화(리샘플 후 세그먼트 표시) =====
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

    % ===== 콘솔 리포트(기존 요건 유지) =====
    report_ads_summary('Whistle 1', pat1, para1);
    report_ads_summary('Whistle 2', pat2, para2);




% ---------- 리샘플링 유틸 (등간격 dx로 보간; x: s, y: Hz) ----------
function [new_x, new_y] = resample_curve(x, y, dx)
    % 정렬 및 중복 제거(보간 안정화)
    [ux, ia] = unique(x(:).', 'stable');  % row vector 유지
    uy = y(ia);
    % 등간격 x 생성
    new_x = ux(1):dx:ux(end);
    % 선형 보간(구간 밖은 생성하지 않음)
    new_y = interp1(ux, uy, new_x, 'linear');
end

% ===== 이하: 이전에 사용하던 생성기/리포트 유틸 그대로 =====

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
    r = dur - p; if r <= 0, error('plateau_ratio가 너무 큼(r>0 필요)'); end

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
