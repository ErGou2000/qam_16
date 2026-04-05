%% CAB_QAM16_DEMO  16-QAM tag modulation/demodulation demo.
%
%  Combines CAB's zero-subcarrier approach with PilotScatter's amplitude
%  tracking and differential demodulation to achieve square 16-QAM for
%  backscatter tag data (4 bits per OFDM symbol).
%
%  This script:
%    1) Sanity test: no channel -> BER = 0
%    2) Moderate channel test: multipath + AWGN + CFO
%    3) Constellation diagram: received vs ideal 16-QAM
%    4) BER vs SNR comparison: 16-QAM vs QPSK tag

clear; clc; close all;
fprintf('============================================================\n');
fprintf('  CAB + PilotScatter: 16-QAM Tag Modulation Demo\n');
fprintf('============================================================\n\n');

%% ==================== 1. Sanity Test (no channel) ====================
fprintf('--- Sanity Test: High SNR, no multipath ---\n');

n_symbols = 50;
mcs = 0;  % BPSK ambient

rng(42);
tx = cab_generate_packet(n_symbols, mcs, true, 42);

% 16-QAM tag modulation
tag16 = cab_tag_modulate_qam16(tx.signal, n_symbols, tx.preamble_len, ...
                                tx.gi_samples, [], 0, 123);

% No channel, tiny noise
ch = cab_channel(tag16.signal, 60, 0, 0, false, 4, 0);

% 16-QAM receiver
r16 = cab_receiver_qam16(ch.signal, n_symbols, mcs, tx.gi_samples);

tag_ber = compute_ber(tag16.tag_bits, r16.tag_bits);
amb_ber = compute_ber(tx.data_bits, r16.ambient_bits);
fprintf('  16-QAM Tag BER = %.6f  |  Ambient BER = %.4f\n', tag_ber, amb_ber);
if tag_ber < 0.01
    fprintf('  [PASS]\n');
else
    fprintf('  [FAIL]\n');
end

%% ==================== 2. Moderate Channel Test ====================
fprintf('\n--- Moderate Channel: SNR=25dB, multipath, CFO=100Hz ---\n');

snr_vals = [30, 25, 20];
for snr = snr_vals
    rng(42);
    tx = cab_generate_packet(n_symbols, mcs, true, 42);
    tag16 = cab_tag_modulate_qam16(tx.signal, n_symbols, tx.preamble_len, ...
                                    tx.gi_samples, [], 0, 123);
    ch = cab_channel(tag16.signal, snr, 100, 0, true, 4, 77);
    r16 = cab_receiver_qam16(ch.signal, n_symbols, mcs, tx.gi_samples);

    tag_ber = compute_ber(tag16.tag_bits, r16.tag_bits);
    amb_ber = compute_ber(tx.data_bits, r16.ambient_bits);
    fprintf('  SNR=%2ddB: 16-QAM Tag BER=%.4f  Ambient BER=%.4f\n', ...
            snr, tag_ber, amb_ber);
end

%% ==================== 3. Constellation Diagram ====================
fprintf('\n--- Generating 16-QAM Constellation Diagram ---\n');

n_symbols = 200;
rng(42);
tx = cab_generate_packet(n_symbols, mcs, true, 42);
tag16 = cab_tag_modulate_qam16(tx.signal, n_symbols, tx.preamble_len, ...
                                tx.gi_samples, [], 0, 123);
ch = cab_channel(tag16.signal, 30, 50, 0, true, 4, 77);
r16 = cab_receiver_qam16(ch.signal, n_symbols, mcs, tx.gi_samples);

% Get ideal constellation
pts = cab_modulation('tag16qam_constellation');

figure('Position', [100, 100, 1200, 500]);

% (a) Received tag symbols (after differential decoding + normalization)
subplot(1, 3, 1);
plot(real(r16.tag_symbols_est), imag(r16.tag_symbols_est), 'b.', 'MarkerSize', 10);
hold on;
plot(real(pts), imag(pts), 'r+', 'MarkerSize', 12, 'LineWidth', 2);
hold off;
xlim([-1.5, 1.5]); ylim([-1.5, 1.5]);
xlabel('In-phase'); ylabel('Quadrature');
title('(a) Demodulated 16-QAM (SNR=30dB)');
legend('Received', 'Ideal', 'Location', 'southeast');
axis square; grid on;

% (b) Differential values D_k
subplot(1, 3, 2);
plot(real(r16.Z_diff), imag(r16.Z_diff), 'g.', 'MarkerSize', 8);
xlim([-3, 3]); ylim([-3, 3]);
xlabel('In-phase'); ylabel('Quadrature');
title('(b) Differential D_k = Z_k / Z_{k-1}');
axis square; grid on;

% (c) Ideal 16-QAM reference
subplot(1, 3, 3);
plot(real(pts), imag(pts), 'r+', 'MarkerSize', 14, 'LineWidth', 2);
hold on;
for i = 1:length(pts)
    text(real(pts(i))+0.05, imag(pts(i))+0.05, dec2bin(i-1,4), 'FontSize', 7);
end
hold off;
xlim([-1.5, 1.5]); ylim([-1.5, 1.5]);
xlabel('In-phase'); ylabel('Quadrature');
title('(c) Ideal 16-QAM constellation');
axis square; grid on;

sgtitle('CAB+PilotScatter: Tag 16-QAM Modulation/Demodulation', 'FontSize', 14);

tag_ber = compute_ber(tag16.tag_bits, r16.tag_bits);
fprintf('  Constellation plot done. BER at SNR=30dB: %.4f\n', tag_ber);

%% ==================== 4. BER vs SNR: 16-QAM vs QPSK ====================
fprintf('\n--- BER vs SNR: 16-QAM tag vs QPSK tag ---\n');

n_symbols = 80;
snr_range = 5:5:40;
n_trials  = 30;

ber_qam16 = zeros(size(snr_range));
ber_qpsk  = zeros(size(snr_range));

for s_idx = 1:length(snr_range)
    snr = snr_range(s_idx);
    bers_16 = zeros(1, n_trials);
    bers_4  = zeros(1, n_trials);
    for trial = 1:n_trials
        sd = trial * 1000;

        % --- 16-QAM tag ---
        tx = cab_generate_packet(n_symbols, 0, true, sd);
        tag16 = cab_tag_modulate_qam16(tx.signal, n_symbols, ...
                    tx.preamble_len, tx.gi_samples, [], 0, sd+1);
        ch = cab_channel(tag16.signal, snr, 50, 0, true, 4, sd+2);
        r16 = cab_receiver_qam16(ch.signal, n_symbols, 0, tx.gi_samples);
        bers_16(trial) = compute_ber(tag16.tag_bits, r16.tag_bits);

        % --- QPSK tag (baseline, using original CAB) ---
        tx2 = cab_generate_packet(n_symbols, 0, true, sd);
        tag4 = cab_tag_modulate(tx2.signal, n_symbols, 4, ...
                    tx2.preamble_len, tx2.gi_samples, [], 0, sd+1);
        ch2 = cab_channel(tag4.signal, snr, 50, 0, true, 4, sd+2);
        r4  = cab_receiver(ch2.signal, n_symbols, 0, tx2.gi_samples, 4);
        bers_4(trial) = compute_ber(tag4.tag_bits, r4.tag_bits);
    end
    ber_qam16(s_idx) = max(mean(bers_16), 1e-5);
    ber_qpsk(s_idx)  = max(mean(bers_4),  1e-5);
    fprintf('  SNR=%2ddB: 16-QAM BER=%.4f  QPSK BER=%.4f\n', ...
            snr, ber_qam16(s_idx), ber_qpsk(s_idx));
end

figure;
semilogy(snr_range, ber_qam16, 'rs-', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
semilogy(snr_range, ber_qpsk,  'bo-', 'LineWidth', 1.5, 'MarkerSize', 8); hold off;
xlabel('SNR (dB)'); ylabel('Tag BER');
title('Tag BER: 16-QAM (4 bits/sym) vs QPSK (2 bits/sym)');
legend('16-QAM (CAB+PilotScatter)', 'QPSK (CAB only)', 'Location', 'southwest');
grid on;

% Throughput comparison
fprintf('\n--- Throughput Comparison ---\n');
sym_rate = 1 / 3.6e-6;  % symbols/sec (with 0.4us GI)
fprintf('  Symbol rate: %.0f sym/s\n', sym_rate);
fprintf('  QPSK:   %.1f kbps (2 bits/sym)\n', sym_rate * 2 / 1000);
fprintf('  16-QAM: %.1f kbps (4 bits/sym) -> 2x throughput gain\n', sym_rate * 4 / 1000);

fprintf('\n============================================================\n');
fprintf('  Demo complete.\n');
fprintf('============================================================\n');

%% ==================== Helper ====================
function ber = compute_ber(tx_bits, rx_bits)
    n = min(length(tx_bits), length(rx_bits));
    if n == 0, ber = 0; return; end
    ber = sum(tx_bits(1:n) ~= rx_bits(1:n)) / n;
end
