function result = cab_receiver_qam16(rx_signal, n_symbols, mcs, gi_samples, ...
                                     phase_method, use_sym_asm, use_phase_track)
% CAB_RECEIVER_QAM16  CAB receiver with 16-QAM tag demodulation.
%
%   Combines:
%     - CAB's WLS zero-subcarrier estimation (phase + amplitude from pilots)
%     - PilotScatter's amplitude tracking (Eq. 13) for amplitude recovery
%     - PilotScatter's differential algorithm (Eq. 7) to eliminate CPE
%     - 16-QAM demapping for tag data (4 bits/symbol)
%
%   Inputs:
%     rx_signal        - Received time-domain signal
%     n_symbols        - Number of payload OFDM symbols
%     mcs              - MCS for ambient data (0-3)
%     gi_samples       - Guard interval samples (8 or 16)
%     phase_method     - Phase estimation: 'wls' (recommended)
%     use_sym_asm      - Enable symbol assembly (default: true)
%     use_phase_track  - Enable phase tracking for ambient (default: true)
%
%   Output: struct with fields:
%     tag_bits, tag_symbols_est, ambient_bits, Z_raw, Z_diff

    if nargin < 3, mcs = 0; end
    if nargin < 4, gi_samples = 8; end
    if nargin < 5, phase_method = 'wls'; end
    if nargin < 6, use_sym_asm = true; end
    if nargin < 7, use_phase_track = true; end

    C = cab_constants();

    %% ====== WiFi Front-End (same as CAB) ======
    front = cab_wlan_toolbox('rx_frontend', rx_signal, mcs, gi_samples);
    rx = front.rx_corrected;
    H = front.channel_est_full;
    phy = local_phy_from_wlan(front);
    gi_samples = phy.gi_samples;

    preamble_len = phy.preamble_len;

    %% ====== Symbol Extraction ======
    freq_symbols = extract_and_equalize(rx, H, preamble_len, n_symbols, ...
                                        gi_samples, zeros(1,n_symbols), phy);

    % Remove known pilot modulation -> complex pilot estimates
    pilot_corrected = remove_pilot_mod(freq_symbols, n_symbols, phy);

    % WLS phase estimation (also get beta1 for SFO)
    [beta0, beta1] = estimate_zero_sc_phase(pilot_corrected, phase_method, phy);

    % Symbol assembly (SFO compensation)
    if use_sym_asm
        offsets = compute_symbol_offsets(beta1, n_symbols, C);
        if any(offsets ~= 0)
            freq_symbols = extract_and_equalize(rx, H, preamble_len, n_symbols, ...
                                                gi_samples, offsets, phy);
            pilot_corrected = remove_pilot_mod(freq_symbols, n_symbols, phy);
            [beta0, beta1] = estimate_zero_sc_phase(pilot_corrected, phase_method, phy);
        end
    end

    %% ====== 16-QAM Tag Demodulation (CAB + PilotScatter fusion) ======
    %
    % Key insight: pilot_corrected(k,:) are complex values containing both
    % amplitude and phase information from the tag modulation.
    %
    % Z_k = mean of pilot_corrected(k,:) is the complex zero-subcarrier
    % estimate that encodes:  Z_k ~ tag_cumulative_gain(k) * CPE_k
    %
    % PilotScatter's differential algorithm (Eq. 7):
    %   D_k = Z_k / Z_{k-1}  cancels CPE since CPE changes slowly
    %   D_k ~ tag_gains(k) = tag_symbols(k) / tag_symbols(k-1)
    %
    % PilotScatter's amplitude tracking (Eq. 13):
    %   |Z_k| / reference_amplitude gives the amplitude level

    % Step 1: Compute complex zero-subcarrier estimate Z_k per symbol
    %  (weighted mean of corrected pilots, using magnitude as weight)
    Z_raw = zeros(1, n_symbols);
    for k = 1:n_symbols
        pc = pilot_corrected(k, :);
        weights = abs(pc);
        Z_raw(k) = sum(weights .* pc) / sum(weights);
    end

    % Step 2: Differential decoding (PilotScatter Eq. 7)
    %  D_k = Z_k / Z_{k-1}
    %  This cancels the slowly-varying CPE: CPE_k / CPE_{k-1} ~ 1
    Z_diff = zeros(1, n_symbols);
    Z_diff(1) = Z_raw(1);  % first symbol: reference is 1 (no previous)
    for k = 2:n_symbols
        if abs(Z_raw(k-1)) > 1e-10
            Z_diff(k) = Z_raw(k) / Z_raw(k-1);
        else
            Z_diff(k) = Z_raw(k);
        end
    end

    % Step 3: Reconstruct absolute QAM symbols from differential gains
    %  cumulative(k) = prod(D_1..D_k) = tag_symbols(k)
    tag_symbols_est = zeros(1, n_symbols);
    tag_symbols_est(1) = Z_diff(1);
    for k = 2:n_symbols
        tag_symbols_est(k) = tag_symbols_est(k-1) * Z_diff(k);
    end

    % Step 4: Amplitude normalization (PilotScatter amplitude tracking)
    %  Use first few symbols' amplitude as reference to normalize.
    %  The 16-QAM constellation has unit average power = 1.
    ref_amp = mean(abs(tag_symbols_est(1:min(5, n_symbols))));
    qam_rms = 1.0;  % 16-QAM normalized to unit power
    if ref_amp > 1e-10
        % Compute expected RMS of 16-QAM constellation
        pts = cab_modulation('tag16qam_constellation');
        qam_rms = sqrt(mean(abs(pts).^2));
        tag_symbols_est = tag_symbols_est * (qam_rms / ref_amp);
    end

    % Step 5: Demap to 16-QAM
    tag_bits = cab_modulation('tag16qam_demap', tag_symbols_est);

    %% ====== Ambient-Data Demodulation (same as CAB) ======
    % Pre-filter: remove zero-subcarrier phase from all subcarriers
    filtered = freq_symbols;
    for k = 1:n_symbols
        filtered(k, :) = filtered(k, :) .* exp(-1i * beta0(k));
    end

    if use_phase_track
        filtered = phase_tracking(filtered, n_symbols, phy);
    end

    qam_orders = [2, 4, 16, 64];
    qam_order  = qam_orders(mcs + 1);
    ambient_bits = [];
    for k = 1:n_symbols
        data_syms = filtered(k, phy.data_fft_bins);
        b = cab_modulation('qam_demap', qam_order, data_syms);
        ambient_bits = [ambient_bits, b]; %#ok<AGROW>
    end

    %% ====== Output ======
    result.tag_bits         = tag_bits;
    result.tag_symbols_est  = tag_symbols_est;
    result.ambient_bits     = ambient_bits;
    result.Z_raw            = Z_raw;
    result.Z_diff           = Z_diff;
    result.beta0            = beta0;
    result.beta1            = beta1;
end

%% ============================================================
%  Internal helper functions (same as cab_receiver.m)
%  ============================================================

function freq_sym = extract_and_equalize(rx, H, preamble_len, n_symbols, ...
                                          gi_samples, offsets, phy)
    symbol_len = gi_samples + phy.n_fft;
    freq_sym = zeros(n_symbols, phy.n_fft);
    for k = 1:n_symbols
        ofs = offsets(k);
        start = preamble_len + (k-1)*symbol_len + ofs;
        data_start = start + gi_samples + 1;
        data_end   = data_start + phy.n_fft - 1;
        if data_end > length(rx), break; end
        time_samples = rx(data_start:data_end);
        freq = fft(time_samples, phy.n_fft);
        for b = 1:phy.n_fft
            if abs(H(b)) > 1e-10
                freq(b) = freq(b) / H(b);
            else
                freq(b) = 0;
            end
        end
        freq_sym(k, :) = freq;
    end
end

function pilot_cor = remove_pilot_mod(freq_symbols, n_symbols, phy)
    pilot_cor = zeros(n_symbols, phy.n_pilot);
    for k = 1:n_symbols
        pol_idx  = mod(k-1, length(phy.pilot_polarity)) + 1;
        polarity = phy.pilot_polarity(pol_idx);
        known_pilots = phy.pilot_base * polarity;
        raw_pilots   = freq_symbols(k, phy.pilot_fft_bins);
        pilot_cor(k, :) = raw_pilots ./ known_pilots;
    end
end

function [beta0, beta1] = estimate_zero_sc_phase(pilot_corrected, method, phy)
    n_sym = size(pilot_corrected, 1);
    x = double(phy.pilot_sc_idx(:));
    beta0 = zeros(1, n_sym);
    beta1 = zeros(1, n_sym);
    for k = 1:n_sym
        A   = abs(pilot_corrected(k, :))';
        Phi = angle(pilot_corrected(k, :))';
        for i = 2:length(Phi)
            while Phi(i) - Phi(1) > pi,  Phi(i) = Phi(i) - 2*pi; end
            while Phi(i) - Phi(1) < -pi, Phi(i) = Phi(i) + 2*pi; end
        end
        switch method
            case 'wls'
                X = [ones(length(x),1), x];
                W = diag(A + 1e-12);
                b_vec = (X'*W*X) \ (X'*W*Phi);
                beta0(k) = b_vec(1);
                beta1(k) = b_vec(2);
            otherwise
                X = [ones(length(x),1), x];
                b_vec = X \ Phi;
                beta0(k) = b_vec(1);
                beta1(k) = b_vec(2);
        end
    end
end

function offsets = compute_symbol_offsets(beta1, n_symbols, C)
    offsets = zeros(1, n_symbols);
    cumul = 0;
    for k = 1:n_symbols
        if abs(beta1(k)) > C.SFO_THRESHOLD
            if beta1(k) > 0, cumul = cumul - 1;
            else,            cumul = cumul + 1;
            end
        end
        offsets(k) = cumul;
    end
end

function filtered = phase_tracking(filtered_in, n_symbols, phy)
    filtered = filtered_in;
    window = 5;
    wrap_fn = @(a) mod(a + pi, 2*pi) - pi;
    for k = 1:n_symbols
        pol_idx = mod(k-1, length(phy.pilot_polarity)) + 1;
        polarity = phy.pilot_polarity(pol_idx);
        known_pilots = phy.pilot_base * polarity;
        pilots_k = filtered(k, phy.pilot_fft_bins) ./ known_pilots;
        lo = max(1, k - floor(window/2));
        hi = min(n_symbols, k + floor(window/2));
        pilot_phases = zeros(1, hi-lo+1);
        idx = 1;
        for m = lo:hi
            pol_m = phy.pilot_polarity(mod(m-1, length(phy.pilot_polarity)) + 1);
            known_m = phy.pilot_base * pol_m;
            p_m = filtered(m, phy.pilot_fft_bins) ./ known_m;
            pilot_phases(idx) = angle(mean(p_m));
            idx = idx + 1;
        end
        avg_phase = mean(pilot_phases);
        x_p = double(phy.pilot_sc_idx(:));
        phi_p = angle(pilots_k(:));
        phi_p = wrap_fn(phi_p - avg_phase);
        p = polyfit(x_p, phi_p, 1);
        slope = p(1);
        for d_idx = 1:phy.n_data
            sc = phy.data_sc_idx(d_idx);
            correction = avg_phase + slope * sc;
            bin = phy.data_fft_bins(d_idx);
            filtered(k, bin) = filtered(k, bin) * exp(-1i * correction);
        end
    end
end

function phy = local_phy_from_wlan(front)
    phy.n_fft = front.n_fft;
    phy.n_data = front.n_data;
    phy.n_pilot = front.n_pilot;
    phy.gi_samples = front.gi_samples;
    phy.preamble_len = front.preamble_len;
    phy.data_fft_bins = front.data_fft_bins;
    phy.pilot_fft_bins = front.pilot_fft_bins;
    phy.data_sc_idx = front.data_sc_idx;
    phy.pilot_sc_idx = front.pilot_sc_idx;
    phy.pilot_base = front.pilot_base;
    phy.pilot_polarity = front.pilot_polarity;
end
