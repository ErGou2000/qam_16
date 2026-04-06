%% CAB_RUN_DEMO  Main demo script for CAB simulation (MATLAB version).
%
%  Runs:
%    1) Sanity tests (no-channel and moderate-channel)
%    2) Figure 6 reproduction: constellation diagram pipeline
%
%  Usage: Run this script directly in MATLAB.

clear; clc; close all;
fprintf('==========================================================\n');
fprintf('  CAB Simulation (MATLAB) - Content-Agnostic Backscatter\n');
fprintf('==========================================================\n\n');

%% ==================== Sanity Tests ====================
fprintf('--- Sanity Tests (high SNR, no serious impairments) ---\n');
all_passed = true;

for psk_order = [2, 4]
    for mcs = [0, 1]
        rng(42);
        tx = cab_generate_packet(30, mcs, true, 42);
        tag_result = cab_tag_modulate(tx.signal, tx.n_symbols, psk_order, ...
                                      tx.preamble_len, tx.gi_samples, [], 0, 123);
        % No channel, just tiny noise
        ch = cab_channel(tag_result.signal, 60, 0, 0, false, 4, 0);

        result = cab_receiver(ch.signal, tx.n_symbols, mcs, tx.gi_samples, ...
                              psk_order, 'wls', true, true, true);

        tag_ber = compute_ber(tag_result.tag_bits, result.tag_bits);
        amb_ber = compute_ber(tx.data_bits, result.ambient_bits);

        if tag_ber < 0.01 && amb_ber < 0.05
            status = 'PASS';
        else
            status = 'FAIL';
            all_passed = false;
        end
        fprintf('  PSK=%2d  MCS=%d  Tag BER=%.4f  Ambient BER=%.4f  [%s]\n', ...
                psk_order, mcs, tag_ber, amb_ber, status);
    end
end

fprintf('\n--- With moderate channel (SNR=25dB, multipath) ---\n');
for psk_order = [2, 4]
    rng(42);
    tx = cab_generate_packet(50, 0, true, 42);
    tag_result = cab_tag_modulate(tx.signal, tx.n_symbols, psk_order, ...
                                  tx.preamble_len, tx.gi_samples, [], 0, 123);
    ch = cab_channel(tag_result.signal, 25, 50, 0, true, 4, 77);

    result = cab_receiver(ch.signal, tx.n_symbols, 0, tx.gi_samples, ...
                          psk_order, 'wls', true, true, true);

    tag_ber = compute_ber(tag_result.tag_bits, result.tag_bits);
    amb_ber = compute_ber(tx.data_bits, result.ambient_bits);

    if tag_ber < 0.1
        status = 'PASS';
    else
        status = 'FAIL';
        all_passed = false;
    end
    fprintf('  PSK=%2d  Tag BER=%.4f  Ambient BER=%.4f  [%s]\n', ...
            psk_order, tag_ber, amb_ber, status);
end

fprintf('\n==========================================================\n');
if all_passed
    fprintf('  Sanity tests: ALL PASSED\n');
else
    fprintf('  Sanity tests: SOME FAILED\n');
end
fprintf('==========================================================\n\n');

%% ==================== Figure 6: Constellation Diagrams ====================
fprintf('--- Generating Figure 6: Constellation Pipeline ---\n');

n_symbols = 100;
psk_order = 4;   % QPSK tag
mcs = 1;         % QPSK ambient

rng(42);
tx = cab_generate_packet(n_symbols, mcs, true, 42);
tag_result = cab_tag_modulate(tx.signal, n_symbols, psk_order, ...
                              tx.preamble_len, tx.gi_samples, [], 0, 123);
ch = cab_channel(tag_result.signal, 25, 200, 2.0, true, 4, 77);
rx_signal = ch.signal;

C = cab_constants();
front = cab_wlan_toolbox('rx_frontend', rx_signal, mcs, tx.gi_samples);
rx = front.rx_corrected;
H = front.channel_est_full;
pilot_fft_bins = front.pilot_fft_bins;
pilot_sc_idx = front.pilot_sc_idx;
pilot_base = front.pilot_base;
pilot_polarity = front.pilot_polarity;
n_fft = front.n_fft;

% Extract symbols
gi_samp = tx.gi_samples;
sym_len = gi_samp + n_fft;
freq_sym = zeros(n_symbols, n_fft);
for k = 1:n_symbols
    ds = tx.preamble_len + (k-1)*sym_len + gi_samp + 1;
    de = ds + n_fft - 1;
    if de > length(rx), break; end
    f = fft(rx(ds:de), n_fft);
    for b = 1:n_fft
        if abs(H(b)) > 1e-10
            f(b) = f(b) / H(b);
        else
            f(b) = 0;
        end
    end
    freq_sym(k,:) = f;
end

% Remove pilot modulation
pilot_cor = zeros(n_symbols, numel(pilot_fft_bins));
for k = 1:n_symbols
    pol_idx = mod(k-1, length(pilot_polarity)) + 1;
    pol = pilot_polarity(pol_idx);
    known = pilot_base * pol;
    pilot_cor(k,:) = freq_sym(k, pilot_fft_bins) ./ known;
end

% Helper: estimate phase with unwrapping
est_phase = @(pc, method) local_est_phase(pc, method, pilot_sc_idx);

% (a) Simple averaging
beta0_mean = est_phase(pilot_cor, 'mean');
pts_a = exp(1i * beta0_mean);

% (b) WLS
beta0_wls = est_phase(pilot_cor, 'wls');
pts_b = exp(1i * beta0_wls);

% (c) CPE separation
tag_ph_c = zeros(1, n_symbols);
cpe_c    = zeros(1, n_symbols);
for i = 1:n_symbols
    if i == 1, cpe_c(i) = 0; else, cpe_c(i) = cpe_c(i-1); end
    raw = beta0_wls(i) - cpe_c(i);
    [snp, ~] = cab_modulation('snap_psk_angle', raw, psk_order);
    tag_ph_c(i) = snp;
    cpe_c(i) = beta0_wls(i) - snp;
end
pts_c = exp(1i * tag_ph_c);

% (d) Symbol assembly: use full receiver
result_full = cab_receiver(rx_signal, n_symbols, mcs, gi_samp, psk_order, ...
                           'wls', true, true, true);
pts_d = exp(1i * result_full.tag_phases);

% --- Plot ---
figure('Position', [100, 100, 1400, 350]);
titles = {'(a) Simple averaging', '(b) WLS', '(c) CPE separation', '(d) Symbol assembly'};
all_pts = {pts_a, pts_b, pts_c, pts_d};

for sp = 1:4
    subplot(1, 4, sp);
    pts = all_pts{sp};
    plot(real(pts), imag(pts), 'b.', 'MarkerSize', 8);
    xlim([-1.5, 1.5]); ylim([-1.5, 1.5]);
    xlabel('In-phase'); ylabel('Quadrature');
    title(titles{sp});
    axis square; grid on;
end
sgtitle('CAB Tag-Data Demodulation Pipeline (Figure 6)', 'FontSize', 14);

fprintf('Done! Figure 6 displayed.\n');
fprintf('==========================================================\n');

%% ==================== Helper Functions ====================

function ber = compute_ber(tx_bits, rx_bits)
    n = min(length(tx_bits), length(rx_bits));
    if n == 0, ber = 0; return; end
    ber = sum(tx_bits(1:n) ~= rx_bits(1:n)) / n;
end

function beta0 = local_est_phase(pilot_cor, method, pilot_sc_idx)
% Estimate zero-subcarrier phase.
    n_sym = size(pilot_cor, 1);
    x = double(pilot_sc_idx(:));
    beta0 = zeros(1, n_sym);
    for k = 1:n_sym
        A   = abs(pilot_cor(k,:))';
        Phi = angle(pilot_cor(k,:))';
        for i = 2:length(Phi)
            while Phi(i)-Phi(1) > pi,  Phi(i) = Phi(i)-2*pi; end
            while Phi(i)-Phi(1) < -pi, Phi(i) = Phi(i)+2*pi; end
        end
        switch method
            case 'mean'
                beta0(k) = mean(Phi);
            case 'wls'
                X = [ones(length(x),1), x];
                W = diag(A + 1e-12);
                b_vec = (X'*W*X) \ (X'*W*Phi);
                beta0(k) = b_vec(1);
        end
    end
end
