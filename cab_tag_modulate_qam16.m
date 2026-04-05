function tag = cab_tag_modulate_qam16(signal, n_symbols, preamble_samples, ...
                                      gi_samples, tag_bits, freq_shift_hz, seed)
% CAB_TAG_MODULATE_QAM16  Apply per-symbol 16-QAM tag modulation.
%
%   Unlike PSK (phase-only), 16-QAM tag modulation applies a complex gain
%   A_k * exp(j*phi_k) per OFDM symbol, modifying both amplitude and phase.
%   This combines CAB's zero-subcarrier approach with PilotScatter's
%   amplitude modulation via variable reflection coefficients.
%
%   Inputs:
%     signal           - Time-domain WiFi signal (row vector)
%     n_symbols        - Number of payload OFDM symbols
%     preamble_samples - Preamble length in samples (default: 400)
%     gi_samples       - Guard interval samples (default: 8)
%     tag_bits         - Tag data bits ([] for random), 4 bits/symbol
%     freq_shift_hz    - Frequency shift in Hz (default: 0)
%     seed             - Random seed
%
%   Output: struct with fields:
%     signal      - Modulated time-domain signal
%     tag_bits    - Tag data bits (4 per symbol)
%     tag_symbols - Complex 16-QAM symbols applied per symbol
%     tag_gains   - Complex gain per symbol (for differential: g_k)

    C = cab_constants();

    if nargin < 3, preamble_samples = 400; end
    if nargin < 4, gi_samples = 8; end
    if nargin < 5, tag_bits = []; end
    if nargin < 6, freq_shift_hz = 0; end
    if nargin < 7, seed = 123; end

    rng(seed);
    bps = 4;  % 16-QAM: 4 bits per symbol
    symbol_len = gi_samples + C.N_FFT;

    % Generate or use provided tag bits
    n_tag_bits = n_symbols * bps;
    if isempty(tag_bits)
        tag_bits = randi([0, 1], 1, n_tag_bits);
    else
        tag_bits = tag_bits(1:min(end, n_tag_bits));
    end

    % Map tag bits to 16-QAM symbols
    tag_symbols = cab_modulation('tag16qam_map', tag_bits);

    % --- Differential encoding (PilotScatter Advanced Onion) ---
    % The tag applies a differential complex gain so that the receiver
    % can use D_k = Z_k / Z_{k-1} to cancel CPE.
    % Reference symbol s_0 = 1 (identity gain for first symbol).
    % For symbol k: actual gain g_k is such that the cumulative product
    % of gains equals tag_symbols(k).
    %
    % cumulative(k) = prod(g_1..g_k) = tag_symbols(k)
    % g_1 = tag_symbols(1) / 1 = tag_symbols(1)
    % g_k = tag_symbols(k) / tag_symbols(k-1)  for k >= 2
    tag_gains = zeros(1, n_symbols);
    tag_gains(1) = tag_symbols(1);  % relative to reference = 1
    for k = 2:n_symbols
        tag_gains(k) = tag_symbols(k) / tag_symbols(k-1);
    end

    % Apply per-symbol complex gain (preamble untouched)
    out = signal;
    cumulative_gain = 1.0;
    for k = 1:n_symbols
        cumulative_gain = cumulative_gain * tag_gains(k);
        s = preamble_samples + (k-1) * symbol_len + 1;
        e = s + symbol_len - 1;
        if e > length(out), break; end
        out(s:e) = out(s:e) .* cumulative_gain;
    end

    % Optional frequency shift
    if freq_shift_hz ~= 0
        t = (0:length(out)-1) / C.SAMPLING_RATE;
        out = out .* exp(1i * 2 * pi * freq_shift_hz * t);
    end

    tag.signal      = out;
    tag.tag_bits    = tag_bits;
    tag.tag_symbols = tag_symbols;
    tag.tag_gains   = tag_gains;
end
