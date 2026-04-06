function tx = cab_generate_packet(n_symbols, mcs, gi_short, seed)
% CAB_GENERATE_PACKET  Generate an HT mixed-format 802.11n-like packet.
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

    rng(seed);
    if gi_short
        gi_samples = 8;
    else
        gi_samples = 16;
    end

    phy = cab_wlan_toolbox('tx_context', mcs, gi_samples);
    gi_samples = phy.gi_samples;
    n_fft = phy.n_fft;
    data_fft_bins = phy.data_fft_bins;
    pilot_fft_bins = phy.pilot_fft_bins;
    pilot_base = phy.pilot_base;
    pilot_polarity = phy.pilot_polarity;
    preamble = phy.preamble;
    preamble_len = phy.preamble_len;
    n_data = phy.n_data;

    % QAM order from MCS
    qam_orders = [2, 4, 16, 64];
    qam_order = qam_orders(mcs + 1);
    bps = log2(qam_order);

    %% Generate payload
    total_data_bits = n_symbols * n_data * bps;
    data_bits = randi([0, 1], 1, total_data_bits);

    data_symbols_all = zeros(n_symbols, n_data);
    pilot_symbols_all = zeros(n_symbols, numel(pilot_fft_bins));
    payload = [];

    for k = 1:n_symbols
        % Map data bits to QAM
        bit_start = (k-1) * n_data * bps + 1;
        bit_end   = k * n_data * bps;
        sym_bits  = data_bits(bit_start:bit_end);
        data_syms = cab_modulation('qam_map', qam_order, sym_bits);
        data_symbols_all(k, :) = data_syms;

        % Pilot values with polarity
        pol_idx = mod(k-1, length(pilot_polarity)) + 1;
        polarity = pilot_polarity(pol_idx);
        pilots = pilot_base * polarity;
        pilot_symbols_all(k, :) = pilots;

        % Build frequency-domain OFDM symbol
        freq = zeros(1, n_fft);
        freq(data_fft_bins) = data_syms;
        freq(pilot_fft_bins) = pilots;

        % IFFT and add cyclic prefix
        time_domain = ifft(freq, n_fft);
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
