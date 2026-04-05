function tag = cab_tag_modulate(signal, n_symbols, psk_order, preamble_samples, ...
                               gi_samples, tag_bits, freq_shift_hz, seed)
% CAB_TAG_MODULATE  Apply per-symbol PSK tag modulation to WiFi signal.
%
%   tag = cab_tag_modulate(signal, n_symbols, psk_order, preamble_samples,
%                          gi_samples, tag_bits, freq_shift_hz, seed)
%
%   Inputs:
%     signal           - Time-domain WiFi signal (row vector)
%     n_symbols        - Number of payload OFDM symbols
%     psk_order        - PSK order: 2, 4, 8, or 16
%     preamble_samples - Preamble length in samples (default: 400)
%     gi_samples       - Guard interval samples (default: 8)
%     tag_bits         - Tag data bits ([] for random)
%     freq_shift_hz    - Frequency shift in Hz (default: 0)
%     seed             - Random seed
%
%   Output: struct with fields: signal, tag_bits, tag_phases

    C = cab_constants();

    if nargin < 3, psk_order = 4; end
    if nargin < 4, preamble_samples = 400; end
    if nargin < 5, gi_samples = 8; end
    if nargin < 6, tag_bits = []; end
    if nargin < 7, freq_shift_hz = 0; end
    if nargin < 8, seed = 123; end

    rng(seed);
    bps = log2(psk_order);
    symbol_len = gi_samples + C.N_FFT;

    % Generate or use provided tag bits
    n_tag_bits = n_symbols * bps;
    if isempty(tag_bits)
        tag_bits = randi([0, 1], 1, n_tag_bits);
    else
        tag_bits = tag_bits(1:min(end, n_tag_bits));
    end

    % Map tag bits to PSK symbols and get phases
    tag_symbols = cab_modulation('psk_map', psk_order, tag_bits);
    tag_phases  = angle(tag_symbols);

    % Apply per-symbol phase multiplication (preamble untouched)
    out = signal;
    for k = 1:n_symbols
        s = preamble_samples + (k-1) * symbol_len + 1;  % 1-based
        e = s + symbol_len - 1;
        if e > length(out), break; end
        out(s:e) = out(s:e) .* exp(1i * tag_phases(k));
    end

    % Optional frequency shift
    if freq_shift_hz ~= 0
        t = (0:length(out)-1) / C.SAMPLING_RATE;
        out = out .* exp(1i * 2 * pi * freq_shift_hz * t);
    end

    tag.signal     = out;
    tag.tag_bits   = tag_bits;
    tag.tag_phases = tag_phases;
end
