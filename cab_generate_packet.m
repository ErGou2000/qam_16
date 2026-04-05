function tx = cab_generate_packet(n_symbols, mcs, gi_short, seed)
% CAB_GENERATE_PACKET  Generate an 802.11n-like OFDM WiFi packet.
%
%   tx = cab_generate_packet(n_symbols, mcs, gi_short, seed)
%
%   Inputs:
%     n_symbols  - Number of payload OFDM symbols (default: 50)
%     mcs        - Modulation coding scheme: 0=BPSK,1=QPSK,2=16QAM,3=64QAM
%     gi_short   - true for 0.4us GI, false for 0.8us GI
%     seed       - Random seed
%
%   Output: struct with fields:
%     signal, data_bits, data_symbols, pilot_symbols,
%     n_symbols, mcs, gi_samples, preamble_len

    if nargin < 1, n_symbols = 50; end
    if nargin < 2, mcs = 0; end
    if nargin < 3, gi_short = true; end
    if nargin < 4, seed = 42; end

    C = cab_constants();
    rng(seed);

    if gi_short
        gi_samples = C.N_GI_SHORT;
    else
        gi_samples = C.N_GI_LONG;
    end

    % QAM order from MCS
    qam_orders = [2, 4, 16, 64];
    qam_order = qam_orders(mcs + 1);
    bps = log2(qam_order);

    %% Generate preamble
    % L-STF: 10 short repetitions = 160 samples
    lstf_time = ifft(C.L_STF_FREQ, C.N_FFT);
    short_sym = lstf_time(1:16);
    lstf = repmat(short_sym, 1, 10);  % 160 samples

    % L-LTF: GI2(32) + LTF(64) + LTF(64) = 160 samples
    lltf_time = ifft(C.L_LTF_FREQ, C.N_FFT);
    gi2 = lltf_time(33:64);  % last 32 samples
    lltf = [gi2, lltf_time, lltf_time];  % 160 samples

    % L-SIG: one OFDM symbol with BPSK (placeholder)
    rng_state = rng; rng(42);
    sig_freq = zeros(1, C.N_FFT);
    sig_data = 2*randi([0,1], 1, C.N_DATA) - 1;  % BPSK: +/-1
    sig_freq(C.DATA_FFT_BINS) = sig_data;
    sig_freq(C.PILOT_FFT_BINS) = C.PILOT_BASE;
    sig_time = ifft(sig_freq, C.N_FFT);
    sig_cp = sig_time(end-C.N_GI_LONG+1 : end);
    lsig = [sig_cp, sig_time];  % 80 samples
    rng(rng_state);  % restore RNG

    preamble = [lstf, lltf, lsig];
    preamble_len = length(preamble);  % 400

    %% Generate payload
    total_data_bits = n_symbols * C.N_DATA * bps;
    data_bits = randi([0, 1], 1, total_data_bits);

    data_symbols_all = zeros(n_symbols, C.N_DATA);
    pilot_symbols_all = zeros(n_symbols, C.N_PILOT);
    payload = [];

    for k = 1:n_symbols
        % Map data bits to QAM
        bit_start = (k-1) * C.N_DATA * bps + 1;
        bit_end   = k * C.N_DATA * bps;
        sym_bits  = data_bits(bit_start:bit_end);
        data_syms = cab_modulation('qam_map', qam_order, sym_bits);
        data_symbols_all(k, :) = data_syms;

        % Pilot values with polarity
        pol_idx = mod(k-1, length(C.PILOT_POLARITY)) + 1;
        polarity = C.PILOT_POLARITY(pol_idx);
        pilots = C.PILOT_BASE * polarity;
        pilot_symbols_all(k, :) = pilots;

        % Build frequency-domain OFDM symbol
        freq = zeros(1, C.N_FFT);
        freq(C.DATA_FFT_BINS)  = data_syms;
        freq(C.PILOT_FFT_BINS) = pilots;

        % IFFT and add cyclic prefix
        time_domain = ifft(freq, C.N_FFT);
        cp = time_domain(end-gi_samples+1 : end);
        payload = [payload, cp, time_domain]; %#ok<AGROW>
    end

    signal = [preamble, payload];

    tx.signal        = signal;
    tx.data_bits     = data_bits;
    tx.data_symbols  = data_symbols_all;
    tx.pilot_symbols = pilot_symbols_all;
    tx.n_symbols     = n_symbols;
    tx.mcs           = mcs;
    tx.gi_samples    = gi_samples;
    tx.preamble_len  = preamble_len;
end
