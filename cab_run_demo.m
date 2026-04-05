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

% --- WiFi front-end ---
% Coarse CFO
lstf = rx_signal(1:160);
corr_c = sum(lstf(17:160) .* conj(lstf(1:144)));
cfo_c  = angle(corr_c) / (2*pi*16/C.SAMPLING_RATE);
rx = rx_signal .* exp(-1i * 2*pi*cfo_c*(0:length(rx_signal)-1)/C.SAMPLING_RATE);

% Fine CFO
ltf1s = 160+32+1;
ltf1 = rx(ltf1s:ltf1s+63);
ltf2 = rx(ltf1s+64:ltf1s+127);
corr_f = sum(ltf2 .* conj(ltf1));
cfo_f  = angle(corr_f) / (2*pi*64/C.SAMPLING_RATE);
rx = rx .* exp(-1i * 2*pi*cfo_f*(0:length(rx)-1)/C.SAMPLING_RATE);

% Channel estimate
ltf1_time = rx(ltf1s:ltf1s+63);
ltf2_time = rx(ltf1s+64:ltf1s+127);
H = zeros(1, C.N_FFT);
L1F = fft(ltf1_time, C.N_FFT);
L2F = fft(ltf2_time, C.N_FFT);
for b = 1:C.N_FFT
    if C.L_LTF_FREQ(b) ~= 0
        H(b) = 0.5*(L1F(b)+L2F(b)) / C.L_LTF_FREQ(b);
    end
end

% Extract symbols
gi_samp = tx.gi_samples;
sym_len = gi_samp + C.N_FFT;
freq_sym = zeros(n_symbols, C.N_FFT);
for k = 1:n_symbols
    ds = tx.preamble_len + (k-1)*sym_len + gi_samp + 1;
    de = ds + C.N_FFT - 1;
    if de > length(rx), break; end
    f = fft(rx(ds:de), C.N_FFT);
    for b = 1:C.N_FFT
        if abs(H(b)) > 1e-10
            f(b) = f(b) / H(b);
        else
            f(b) = 0;
        end
    end
    freq_sym(k,:) = f;
end

% Remove pilot modulation
pilot_cor = zeros(n_symbols, C.N_PILOT);
for k = 1:n_symbols
    pol_idx = mod(k-1, length(C.PILOT_POLARITY)) + 1;
    pol = C.PILOT_POLARITY(pol_idx);
    known = C.PILOT_BASE * pol;
    pilot_cor(k,:) = freq_sym(k, C.PILOT_FFT_BINS) ./ known;
end

% Helper: estimate phase with unwrapping
est_phase = @(pc, method) local_est_phase(pc, method, C);

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

function beta0 = local_est_phase(pilot_cor, method, C)
% Estimate zero-subcarrier phase.
    n_sym = size(pilot_cor, 1);
    x = double(C.PILOT_SC_IDX(:));
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
