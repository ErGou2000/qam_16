function varargout = cab_wlan_toolbox(action, varargin)
% CAB_WLAN_TOOLBOX  Bridge CAB processing to MATLAB WLAN Toolbox HT helpers.
%
%   ctx = cab_wlan_toolbox('tx_context', mcs, gi_samples)
%   front = cab_wlan_toolbox('rx_frontend', rx_signal, mcs, gi_samples)

    switch action
        case 'tx_context'
            varargout{1} = local_tx_context(varargin{:});
        case 'rx_frontend'
            varargout{1} = local_rx_frontend(varargin{:});
        otherwise
            error('Unknown action: %s', action);
    end
end

function ctx = local_tx_context(mcs, gi_samples)
    cfg = local_create_cfg(mcs, gi_samples);
    data_idx = wlanFieldIndices(cfg, 'HT-Data');
    bits = zeros(8 * cfg.PSDULength, 1, 'int8');
    waveform = wlanWaveformGenerator(bits, cfg);
    ofdm_info = wlanHTOFDMInfo('HT-Data', cfg);
    hltf_info = wlanHTOFDMInfo('HT-LTF', cfg);

    active_fft_bins = double(ofdm_info.ActiveFFTIndices(:)).';
    active_sc_idx = double(ofdm_info.ActiveFrequencyIndices(:)).';
    data_active_idx = double(ofdm_info.DataIndices(:)).';
    pilot_active_idx = double(ofdm_info.PilotIndices(:)).';

    ctx.cfg = cfg;
    ctx.sample_rate = ofdm_info.SampleRate;
    ctx.n_fft = ofdm_info.FFTLength;
    ctx.gi_samples = ofdm_info.CPLength;
    ctx.preamble = waveform(1:data_idx(1)-1, 1).';
    ctx.preamble_len = data_idx(1) - 1;
    ctx.data_field_indices = double(data_idx(:)).';
    ctx.active_fft_bins = active_fft_bins;
    ctx.active_sc_idx = active_sc_idx;
    ctx.data_active_idx = data_active_idx;
    ctx.pilot_active_idx = pilot_active_idx;
    ctx.data_fft_bins = active_fft_bins(data_active_idx);
    ctx.pilot_fft_bins = active_fft_bins(pilot_active_idx);
    ctx.data_sc_idx = active_sc_idx(data_active_idx);
    ctx.pilot_sc_idx = active_sc_idx(pilot_active_idx);
    ctx.n_data = numel(ctx.data_fft_bins);
    ctx.n_pilot = numel(ctx.pilot_fft_bins);
    ctx.pilot_base = [1, 1, 1, -1];
    ctx.pilot_polarity = cab_constants().PILOT_POLARITY;
    ctx.htltf_active_fft_bins = double(hltf_info.ActiveFFTIndices(:)).';
end

function front = local_rx_frontend(rx_signal, mcs, gi_samples)
    ctx = local_tx_context(mcs, gi_samples);

    lstf_idx = wlanFieldIndices(ctx.cfg, 'L-STF');
    lltf_idx = wlanFieldIndices(ctx.cfg, 'L-LTF');

    lstf = rx_signal(lstf_idx(1):min(lstf_idx(2), length(rx_signal)));
    coarse_cfo = wlanCoarseCFOEstimate(lstf(:), ctx.cfg.ChannelBandwidth);
    rx = local_correct_cfo(rx_signal, coarse_cfo, ctx.sample_rate);

    lltf = rx(lltf_idx(1):min(lltf_idx(2), length(rx)));
    fine_cfo = wlanFineCFOEstimate(lltf(:), ctx.cfg.ChannelBandwidth);
    rx = local_correct_cfo(rx, fine_cfo, ctx.sample_rate);

    htltf_idx = wlanFieldIndices(ctx.cfg, 'HT-LTF');
    htltf = rx(htltf_idx(1):min(htltf_idx(2), length(rx)));
    demod_htltf = wlanHTLTFDemodulate(htltf(:), ctx.cfg);
    ch_est_active = wlanHTLTFChannelEstimate(demod_htltf, ctx.cfg);
    ch_est_active = reshape(ch_est_active(:, 1, 1), [], 1);

    h_full = zeros(1, ctx.n_fft);
    h_full(ctx.htltf_active_fft_bins) = ch_est_active;

    front = ctx;
    front.rx_corrected = reshape(rx, 1, []);
    front.coarse_cfo = coarse_cfo;
    front.fine_cfo = fine_cfo;
    front.channel_est_active = ch_est_active(:).';
    front.channel_est_full = h_full;
end

function cfg = local_create_cfg(mcs, gi_samples)
    cfg = wlanHTConfig;
    if isprop(cfg, 'ChannelBandwidth')
        cfg.ChannelBandwidth = 'CBW20';
    end
    if isprop(cfg, 'MCS')
        cfg.MCS = local_map_mcs(mcs);
    end
    if isprop(cfg, 'PSDULength')
        cfg.PSDULength = 1;
    end
    if isprop(cfg, 'GuardInterval')
        if gi_samples <= 8
            cfg.GuardInterval = 'Short';
        else
            cfg.GuardInterval = 'Long';
        end
    end
    if isprop(cfg, 'NumTransmitAntennas')
        cfg.NumTransmitAntennas = 1;
    end
    if isprop(cfg, 'NumSpaceTimeStreams')
        cfg.NumSpaceTimeStreams = 1;
    end
end

function mcs_ht = local_map_mcs(mcs_cab)
    % Preserve the original "modulation order only" abstraction while
    % selecting HT MCS values with matching constellations for one stream.
    mcs_map = [0, 1, 3, 5];
    idx = min(max(round(mcs_cab), 0), numel(mcs_map) - 1) + 1;
    mcs_ht = mcs_map(idx);
end

function out = local_correct_cfo(sig, cfo_hz, sample_rate)
    n = 0:length(sig)-1;
    out = sig .* exp(-1i * 2 * pi * cfo_hz * n / sample_rate);
end
