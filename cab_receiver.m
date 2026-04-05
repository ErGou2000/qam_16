function result = cab_receiver(rx_signal, n_symbols, mcs, gi_samples, psk_order, ...
                              phase_method, use_cpe_sep, use_sym_asm, use_phase_track)
% CAB_RECEIVER  Full CAB receiver: demodulate tag and ambient data.
%
%   result = cab_receiver(rx_signal, n_symbols, mcs, gi_samples, psk_order,
%                         phase_method, use_cpe_sep, use_sym_asm, use_phase_track)
%
%   Inputs:
%     rx_signal        - Received time-domain signal (row vector)
%     n_symbols        - Number of payload OFDM symbols
%     mcs              - Modulation coding scheme (0-3)
%     gi_samples       - Guard interval samples (8 or 16)
%     psk_order        - Tag PSK order (2, 4, 8, 16)
%     phase_method     - Phase estimation: 'mean','wmean','ls','wls','spline'
%     use_cpe_sep      - Enable CPE separation (default: true)
%     use_sym_asm      - Enable symbol assembly (default: true)
%     use_phase_track  - Enable phase tracking for ambient (default: true)
%
%   Output: struct with fields:
%     tag_bits, tag_phases, ambient_bits, cpe_estimates, beta0, beta1

    if nargin < 3,  mcs = 0; end
    if nargin < 4,  gi_samples = 8; end
    if nargin < 5,  psk_order = 4; end
    if nargin < 6,  phase_method = 'wls'; end
    if nargin < 7,  use_cpe_sep = true; end
    if nargin < 8,  use_sym_asm = true; end
    if nargin < 9,  use_phase_track = true; end

    C = cab_constants();
    preamble_len = 400;
    symbol_len = gi_samples + C.N_FFT;

    %% ====== WiFi Front-End ======

    % Coarse CFO from L-STF (16-sample repetition)
    lstf = rx_signal(1:160);
    corr_coarse = sum(lstf(17:160) .* conj(lstf(1:144)));
    cfo_coarse  = angle(corr_coarse) / (2*pi * 16 / C.SAMPLING_RATE);
    rx = correct_cfo(rx_signal, cfo_coarse, C);

    % Fine CFO from L-LTF (two 64-sample repetitions)
    ltf1_start = 160 + 32 + 1;  % 1-based
    ltf1 = rx(ltf1_start : ltf1_start+63);
    ltf2 = rx(ltf1_start+64 : ltf1_start+127);
    corr_fine = sum(ltf2 .* conj(ltf1));
    cfo_fine  = angle(corr_fine) / (2*pi * 64 / C.SAMPLING_RATE);
    rx = correct_cfo(rx, cfo_fine, C);

    % Channel estimation from L-LTF
    H = channel_estimate(rx, C);

    %% ====== Initial Symbol Extraction ======
    freq_symbols = extract_and_equalize(rx, H, preamble_len, n_symbols, ...
                                        gi_samples, zeros(1,n_symbols), C);

    %% ====== Tag-Data Demodulation (Section 3.2) ======

    % Remove known pilot modulation
    pilot_corrected = remove_pilot_mod(freq_symbols, n_symbols, C);

    % Phase estimation for zero-subcarriers
    [beta0, beta1] = estimate_zero_sc_phase(pilot_corrected, phase_method, C);

    % Symbol assembly (SFO compensation)
    if use_sym_asm
        offsets = compute_symbol_offsets(beta1, n_symbols, C);
        if any(offsets ~= 0)
            freq_symbols = extract_and_equalize(rx, H, preamble_len, n_symbols, ...
                                                gi_samples, offsets, C);
            pilot_corrected = remove_pilot_mod(freq_symbols, n_symbols, C);
            [beta0, beta1] = estimate_zero_sc_phase(pilot_corrected, phase_method, C);
        end
    end

    % CPE separation
    if use_cpe_sep
        [tag_phases, cpe_est] = cpe_separation(beta0, psk_order);
    else
        tag_phases = zeros(1, n_symbols);
        cpe_est    = zeros(1, n_symbols);
        for i = 1:n_symbols
            [tag_phases(i), ~] = cab_modulation('snap_psk_angle', beta0(i), psk_order);
        end
    end

    % Tag phases to bits
    tag_bits = tag_phases_to_bits(tag_phases, psk_order);

    %% ====== Ambient-Data Demodulation (Section 3.3) ======

    % Pre-filter: remove zero-subcarrier phase from all subcarriers
    filtered = freq_symbols;
    for k = 1:n_symbols
        filtered(k, :) = filtered(k, :) .* exp(-1i * beta0(k));
    end

    % Phase tracking
    if use_phase_track
        filtered = phase_tracking(filtered, n_symbols, C);
    end

    % Demap ambient data
    qam_orders = [2, 4, 16, 64];
    qam_order  = qam_orders(mcs + 1);
    ambient_bits = [];
    for k = 1:n_symbols
        data_syms = filtered(k, C.DATA_FFT_BINS);
        b = cab_modulation('qam_demap', qam_order, data_syms);
        ambient_bits = [ambient_bits, b]; %#ok<AGROW>
    end

    %% ====== Output ======
    result.tag_bits      = tag_bits;
    result.tag_phases    = tag_phases;
    result.ambient_bits  = ambient_bits;
    result.cpe_estimates = cpe_est;
    result.beta0         = beta0;
    result.beta1         = beta1;
    result.pilot_corrected = pilot_corrected;
end

%% ============================================================
%  Internal helper functions
%  ============================================================

function out = correct_cfo(sig, cfo_hz, C)
    n = 0:length(sig)-1;
    out = sig .* exp(-1i * 2*pi * cfo_hz * n / C.SAMPLING_RATE);
end

function H = channel_estimate(rx, C)
% Estimate channel from two L-LTF repetitions.
    ltf1_start = 160 + 32 + 1;  % 1-based
    ltf1_time = rx(ltf1_start     : ltf1_start+63);
    ltf2_time = rx(ltf1_start+64  : ltf1_start+127);
    ltf1_freq = fft(ltf1_time, C.N_FFT);
    ltf2_freq = fft(ltf2_time, C.N_FFT);
    H = zeros(1, C.N_FFT);
    for b = 1:C.N_FFT
        if C.L_LTF_FREQ(b) ~= 0
            H(b) = 0.5 * (ltf1_freq(b) + ltf2_freq(b)) / C.L_LTF_FREQ(b);
        end
    end
end

function freq_sym = extract_and_equalize(rx, H, preamble_len, n_symbols, ...
                                          gi_samples, offsets, C)
% Extract payload symbols, remove CP, FFT, equalize.
    symbol_len = gi_samples + C.N_FFT;
    freq_sym = zeros(n_symbols, C.N_FFT);
    for k = 1:n_symbols
        ofs = offsets(k);
        start = preamble_len + (k-1)*symbol_len + ofs;  % 0-based sample position
        data_start = start + gi_samples + 1;             % 1-based MATLAB index
        data_end   = data_start + C.N_FFT - 1;
        if data_end > length(rx), break; end
        time_samples = rx(data_start:data_end);
        freq = fft(time_samples, C.N_FFT);
        % Zero-forcing equalization
        for b = 1:C.N_FFT
            if abs(H(b)) > 1e-10
                freq(b) = freq(b) / H(b);
            else
                freq(b) = 0;
            end
        end
        freq_sym(k, :) = freq;
    end
end

function pilot_cor = remove_pilot_mod(freq_symbols, n_symbols, C)
% Remove known pilot signs from pilot subcarriers.
    pilot_cor = zeros(n_symbols, C.N_PILOT);
    for k = 1:n_symbols
        pol_idx  = mod(k-1, length(C.PILOT_POLARITY)) + 1;
        polarity = C.PILOT_POLARITY(pol_idx);
        known_pilots = C.PILOT_BASE * polarity;
        raw_pilots   = freq_symbols(k, C.PILOT_FFT_BINS);
        pilot_cor(k, :) = raw_pilots ./ known_pilots;
    end
end

function [beta0, beta1] = estimate_zero_sc_phase(pilot_corrected, method, C)
% Estimate zero-subcarrier phase via various methods (Eq. 3).
    n_sym = size(pilot_corrected, 1);
    x = double(C.PILOT_SC_IDX(:));   % column vector
    beta0 = zeros(1, n_sym);
    beta1 = zeros(1, n_sym);

    for k = 1:n_sym
        A   = abs(pilot_corrected(k, :))';    % column
        Phi = angle(pilot_corrected(k, :))';   % column

        % Phase unwrapping relative to first pilot
        for i = 2:length(Phi)
            while Phi(i) - Phi(1) > pi,  Phi(i) = Phi(i) - 2*pi; end
            while Phi(i) - Phi(1) < -pi, Phi(i) = Phi(i) + 2*pi; end
        end

        switch method
            case 'mean'
                beta0(k) = mean(Phi);
                beta1(k) = 0;
            case 'wmean'
                beta0(k) = sum(A .* Phi) / sum(A);
                beta1(k) = 0;
            case 'ls'
                X = [ones(length(x),1), x];
                b_vec = X \ Phi;
                beta0(k) = b_vec(1);
                beta1(k) = b_vec(2);
            case 'wls'
                X = [ones(length(x),1), x];
                W = diag(A + 1e-12);
                XtWX   = X' * W * X;
                XtWPhi = X' * W * Phi;
                b_vec  = XtWX \ XtWPhi;
                beta0(k) = b_vec(1);
                beta1(k) = b_vec(2);
            case 'spline'
                beta0(k) = interp1(x, Phi, 0, 'spline');
                % Numerical derivative at 0
                dx = 0.1;
                beta1(k) = (interp1(x,Phi,dx,'spline') - interp1(x,Phi,-dx,'spline'))/(2*dx);
            otherwise
                error('Unknown method: %s', method);
        end
    end
end

function [tag_phases, cpe] = cpe_separation(beta0, psk_order)
% Iterative CPE separation (Section 3.2).
    n_sym = length(beta0);
    tag_phases = zeros(1, n_sym);
    cpe = zeros(1, n_sym);

    for i = 1:n_sym
        % Step 1: assign CPE from previous
        if i == 1
            cpe(i) = 0;
        else
            cpe(i) = cpe(i-1);
        end
        % Step 2: estimate tag phase, snap to nearest PSK point
        raw_tag = beta0(i) - cpe(i);
        [snapped, ~] = cab_modulation('snap_psk_angle', raw_tag, psk_order);
        tag_phases(i) = snapped;
        % Step 3: update CPE
        cpe(i) = beta0(i) - tag_phases(i);
    end
end

function offsets = compute_symbol_offsets(beta1, n_symbols, C)
% Compute per-symbol sample offsets for SFO compensation.
    offsets = zeros(1, n_symbols);
    cumul = 0;
    for k = 1:n_symbols
        if abs(beta1(k)) > C.SFO_THRESHOLD
            if beta1(k) > 0
                cumul = cumul - 1;
            else
                cumul = cumul + 1;
            end
        end
        offsets(k) = cumul;
    end
end

function bits = tag_phases_to_bits(tag_phases, psk_order)
% Convert demodulated tag phases to bits.
    n_sym = length(tag_phases);
    syms = exp(1i * tag_phases);
    bits = cab_modulation('psk_demap', psk_order, syms);
end

function filtered = phase_tracking(filtered_in, n_symbols, C)
% Customized phase tracking for ambient data (Section 3.3).
    filtered = filtered_in;
    window = 5;
    wrap = @(a) mod(a + pi, 2*pi) - pi;

    for k = 1:n_symbols
        % Current symbol pilots (remove known polarity)
        pol_idx = mod(k-1, length(C.PILOT_POLARITY)) + 1;
        polarity = C.PILOT_POLARITY(pol_idx);
        known_pilots = C.PILOT_BASE * polarity;
        pilots_k = filtered(k, C.PILOT_FFT_BINS) ./ known_pilots;

        % Sliding window pilot averaging
        lo = max(1, k - floor(window/2));
        hi = min(n_symbols, k + floor(window/2));
        pilot_phases = zeros(1, hi-lo+1);
        idx = 1;
        for m = lo:hi
            pol_m = C.PILOT_POLARITY(mod(m-1, length(C.PILOT_POLARITY)) + 1);
            known_m = C.PILOT_BASE * pol_m;
            p_m = filtered(m, C.PILOT_FFT_BINS) ./ known_m;
            pilot_phases(idx) = angle(mean(p_m));
            idx = idx + 1;
        end
        avg_phase = mean(pilot_phases);

        % Linear fit across pilot positions for SFO tracking
        x_p = double(C.PILOT_SC_IDX(:));
        phi_p = angle(pilots_k(:));
        phi_p = wrap(phi_p - avg_phase);
        p = polyfit(x_p, phi_p, 1);
        slope = p(1);

        % Correct data subcarriers
        for d_idx = 1:C.N_DATA
            sc = C.DATA_SC_IDX(d_idx);
            correction = avg_phase + slope * sc;
            bin = C.DATA_FFT_BINS(d_idx);
            filtered(k, bin) = filtered(k, bin) * exp(-1i * correction);
        end
    end
end
