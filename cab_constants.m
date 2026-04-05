function C = cab_constants()
% CAB_CONSTANTS  Return struct with all IEEE 802.11n OFDM parameters.
%
%   C = cab_constants()
%
%   Returns a struct C containing:
%     N_FFT, N_DATA, N_PILOT, SAMPLING_RATE, N_GI_SHORT, N_GI_LONG,
%     PILOT_SC_IDX, PILOT_BASE, PILOT_FFT_BINS, DATA_SC_IDX, DATA_FFT_BINS,
%     PILOT_POLARITY, L_STF_FREQ, L_LTF_FREQ, SFO_THRESHOLD

    C.N_FFT  = 64;
    C.N_DATA = 48;
    C.N_PILOT = 4;
    C.N_USED = 52;
    C.SAMPLING_RATE = 20e6;
    C.T_FFT = C.N_FFT / C.SAMPLING_RATE;  % 3.2 us
    C.N_GI_SHORT = 8;   % 0.4 us
    C.N_GI_LONG  = 16;  % 0.8 us

    % Pilot subcarrier indices (centered around DC)
    C.PILOT_SC_IDX  = [-21, -7, 7, 21];
    C.PILOT_BASE    = [1, 1, 1, -1];

    % Data subcarrier indices (48 subcarriers, excluding pilots/DC/guard)
    all_used = [-26:-1, 1:26];  % 52 values, no DC
    pilot_set = C.PILOT_SC_IDX;
    mask = true(1, length(all_used));
    for p = pilot_set
        mask = mask & (all_used ~= p);
    end
    C.DATA_SC_IDX = all_used(mask);  % 48 data subcarrier indices

    % Convert subcarrier indices to 1-based FFT bin indices
    % sc_to_bin: centered index k -> MATLAB FFT bin = mod(k, N_FFT) + 1
    C.PILOT_FFT_BINS = mod(C.PILOT_SC_IDX, C.N_FFT) + 1;
    C.DATA_FFT_BINS  = mod(C.DATA_SC_IDX,  C.N_FFT) + 1;

    % SFO threshold (Section 3.2 footnote 2)
    C.SFO_THRESHOLD = 0.9 * 2 * pi / C.N_FFT;

    % Pilot polarity sequence (127 elements, Table 20-17 in 802.11n)
    C.PILOT_POLARITY = [ ...
         1, 1, 1, 1,-1,-1,-1, 1,-1,-1,-1,-1, 1, 1,-1, 1, ...
        -1,-1, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 1,-1, 1, ...
         1, 1,-1, 1, 1,-1,-1, 1, 1, 1,-1, 1,-1,-1,-1, 1, ...
        -1, 1,-1,-1, 1,-1,-1, 1, 1, 1, 1, 1,-1,-1, 1, 1, ...
        -1,-1, 1,-1, 1,-1, 1, 1,-1,-1,-1, 1, 1,-1,-1,-1, ...
        -1, 1,-1,-1, 1,-1, 1, 1, 1, 1,-1, 1,-1, 1,-1, 1, ...
        -1,-1,-1,-1,-1, 1,-1, 1, 1,-1, 1,-1, 1, 1, 1,-1, ...
        -1, 1,-1,-1,-1, 1, 1, 1,-1,-1,-1,-1,-1,-1,-1];

    % L-STF frequency-domain sequence (64 bins, 1-based)
    C.L_STF_FREQ = zeros(1, C.N_FFT);
    lstf_map_sc  = [-24,-20,-16,-12, -8, -4,  4,  8, 12, 16, 20, 24];
    lstf_map_val = [-1-1i, 1-1i, 1+1i, -1+1i, 1+1i, -1-1i, ...
                    -1-1i, 1-1i, 1+1i, -1+1i, 1+1i, -1-1i];
    for idx = 1:length(lstf_map_sc)
        bin = mod(lstf_map_sc(idx), C.N_FFT) + 1;
        C.L_STF_FREQ(bin) = lstf_map_val(idx) * sqrt(13/6);
    end

    % L-LTF frequency-domain sequence (64 bins, 1-based)
    % Subcarriers -26..+26 (53 values), DC=0
    lltf_vals = [ ...
        0, 1,-1,-1, 1, 1,-1, 1,-1, 1,-1,-1,-1,-1,-1, 1, ...
        1,-1,-1, 1,-1, 1,-1, 1, 1, 1, 1, ...  % sc -26..-1 (first is index for -26)
        0, ...                                   % DC
        1, 1,-1,-1, 1, 1,-1, 1,-1, 1, 1, 1, 1, 1, 1,-1, ...
       -1, 1, 1,-1, 1,-1, 1, 1, 1, 1];         % sc 1..26
    C.L_LTF_FREQ = zeros(1, C.N_FFT);
    for i = 1:length(lltf_vals)
        sc = i - 27;  % subcarrier index: -26..+26
        bin = mod(sc, C.N_FFT) + 1;
        C.L_LTF_FREQ(bin) = lltf_vals(i);
    end

end
