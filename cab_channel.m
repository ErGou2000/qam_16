function ch = cab_channel(signal, snr_db, cfo_hz, sfo_ppm, multipath_en, n_taps, seed)
% CAB_CHANNEL  Apply wireless channel impairments.
%
%   ch = cab_channel(signal, snr_db, cfo_hz, sfo_ppm, multipath_en, n_taps, seed)
%
%   Order: multipath -> CFO -> SFO -> AWGN
%
%   Inputs:
%     signal       - Input signal (row vector)
%     snr_db       - SNR in dB (default: 30)
%     cfo_hz       - Carrier frequency offset in Hz (default: 0)
%     sfo_ppm      - Sampling frequency offset in ppm (default: 0)
%     multipath_en - Enable multipath fading (default: true)
%     n_taps       - Number of multipath taps (default: 4)
%     seed         - Random seed
%
%   Output: struct with fields: signal, channel_taps

    C = cab_constants();

    if nargin < 2, snr_db = 30; end
    if nargin < 3, cfo_hz = 0; end
    if nargin < 4, sfo_ppm = 0; end
    if nargin < 5, multipath_en = true; end
    if nargin < 6, n_taps = 4; end
    if nargin < 7, seed = 77; end

    rng(seed);
    out = signal;
    channel_taps = [];

    %% Multipath
    if multipath_en
        delay_profiles = {'Model-A', 'Model-B', 'Model-C', 'Model-D', 'Model-E', 'Model-F'};
        profile_idx = min(max(round(n_taps), 1), numel(delay_profiles));
        tgn = wlanTGnChannel( ...
            'SampleRate', C.SAMPLING_RATE, ...
            'DelayProfile', delay_profiles{profile_idx}, ...
            'LargeScaleFadingEffect', 'None', ...
            'NormalizeChannelOutputs', true, ...
            'RandomStream', 'mt19937ar with seed', ...
            'Seed', seed);
        out = tgn(out(:));
        out = reshape(out, 1, []);
    end

    %% CFO
    if abs(cfo_hz) > 0
        n = 0:length(out)-1;
        out = out .* exp(1i * 2 * pi * cfo_hz * n / C.SAMPLING_RATE);
    end

    %% SFO
    if abs(sfo_ppm) > 1e-10
        n_orig = length(out);
        ratio = 1.0 + sfo_ppm * 1e-6;
        n_new = round(n_orig * ratio);
        t_orig = 0:n_orig-1;
        t_new  = linspace(0, n_orig-1, n_new);

        re_interp = interp1(t_orig, real(out), t_new, 'spline', 0);
        im_interp = interp1(t_orig, imag(out), t_new, 'spline', 0);
        out_sfo = re_interp + 1i * im_interp;

        if length(out_sfo) >= n_orig
            out = out_sfo(1:n_orig);
        else
            out = [out_sfo, zeros(1, n_orig - length(out_sfo))];
        end
    end

    %% AWGN
    sig_power   = mean(abs(out).^2);
    noise_power = sig_power / (10^(snr_db/10));
    noise = sqrt(noise_power/2) * (randn(1,length(out)) + 1i*randn(1,length(out)));
    out = out + noise;

    ch.signal       = out;
    ch.channel_taps = channel_taps;
end
