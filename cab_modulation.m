function varargout = cab_modulation(action, varargin)
% CAB_MODULATION  QAM and PSK modulation/demodulation functions.
%
%   symbols = cab_modulation('qam_map',  order, bits)
%   bits    = cab_modulation('qam_demap', order, symbols)
%   symbols = cab_modulation('psk_map',  order, bits)
%   bits    = cab_modulation('psk_demap', order, symbols)
%   [snapped_angle, idx] = cab_modulation('snap_psk_angle', angle_rad, psk_order)
%   symbols = cab_modulation('tag16qam_map', bits)     % 16-QAM for tag
%   bits    = cab_modulation('tag16qam_demap', symbols) % 16-QAM demap
%   snapped = cab_modulation('tag16qam_snap', symbols)  % snap to nearest
%   points  = cab_modulation('tag16qam_constellation')  % get constellation
%
%   order: 2=BPSK, 4=QPSK, 8=8PSK, 16=16QAM/16PSK, 64=64QAM

    switch action
        case 'qam_map'
            varargout{1} = qam_map(varargin{1}, varargin{2});
        case 'qam_demap'
            varargout{1} = qam_demap(varargin{1}, varargin{2});
        case 'psk_map'
            varargout{1} = psk_map(varargin{1}, varargin{2});
        case 'psk_demap'
            varargout{1} = psk_demap(varargin{1}, varargin{2});
        case 'snap_psk_angle'
            [varargout{1}, varargout{2}] = snap_psk_angle(varargin{1}, varargin{2});
        case 'tag16qam_map'
            varargout{1} = tag16qam_map(varargin{1});
        case 'tag16qam_demap'
            varargout{1} = tag16qam_demap(varargin{1});
        case 'tag16qam_snap'
            varargout{1} = tag16qam_snap(varargin{1});
        case 'tag16qam_constellation'
            varargout{1} = tag16qam_constellation();
        otherwise
            error('Unknown action: %s', action);
    end
end

%% ---- QAM ----

function points = qam_constellation(order)
    switch order
        case 2   % BPSK
            points = [-1, 1];
        case 4   % QPSK
            points = [-1-1i, -1+1i, 1-1i, 1+1i] / sqrt(2);
        case 16  % 16-QAM
            vals = [-3, -1, 1, 3];
            points = zeros(1, 16);
            idx = 1;
            for q = vals
                for ii = vals
                    points(idx) = complex(ii, q) / sqrt(10);
                    idx = idx + 1;
                end
            end
        case 64  % 64-QAM
            vals = [-7,-5,-3,-1,1,3,5,7];
            points = zeros(1, 64);
            idx = 1;
            for q = vals
                for ii = vals
                    points(idx) = complex(ii, q) / sqrt(42);
                    idx = idx + 1;
                end
            end
        otherwise
            error('Unsupported QAM order: %d', order);
    end
end

function symbols = qam_map(order, bits)
    pts = qam_constellation(order);
    bps = log2(order);
    n_sym = floor(length(bits) / bps);
    symbols = zeros(1, n_sym);
    for i = 1:n_sym
        b = bits((i-1)*bps+1 : i*bps);
        dec_val = bin2dec_vec(b);
        gray_dec = gray_decode(dec_val);
        symbols(i) = pts(gray_dec + 1);  % +1 for MATLAB indexing
    end
end

function bits = qam_demap(order, symbols)
    pts = qam_constellation(order);
    bps = log2(order);
    n_sym = length(symbols);
    bits = zeros(1, n_sym * bps);
    for i = 1:n_sym
        dists = abs(symbols(i) - pts);
        [~, idx] = min(dists);
        idx = idx - 1;  % 0-based constellation index
        gray_idx = gray_encode(idx);
        b = dec2bin_vec(gray_idx, bps);
        bits((i-1)*bps+1 : i*bps) = b;
    end
end

%% ---- PSK ----

function points = psk_constellation(order)
    points = exp(1i * 2 * pi * (0:order-1) / order);
end

function symbols = psk_map(order, bits)
    pts = psk_constellation(order);
    bps = log2(order);
    n_sym = floor(length(bits) / bps);
    symbols = zeros(1, n_sym);
    for i = 1:n_sym
        b = bits((i-1)*bps+1 : i*bps);
        dec_val = bin2dec_vec(b);
        gray_dec = gray_decode(dec_val);
        symbols(i) = pts(gray_dec + 1);
    end
end

function bits = psk_demap(order, symbols)
    pts = psk_constellation(order);
    bps = log2(order);
    n_sym = length(symbols);
    bits = zeros(1, n_sym * bps);
    for i = 1:n_sym
        dists = abs(symbols(i) - pts);
        [~, idx] = min(dists);
        idx = idx - 1;  % 0-based
        gray_idx = gray_encode(idx);
        b = dec2bin_vec(gray_idx, bps);
        bits((i-1)*bps+1 : i*bps) = b;
    end
end

function [snapped, idx] = snap_psk_angle(angle_rad, psk_order)
% Snap angle to nearest M-PSK constellation angle.
    step = 2 * pi / psk_order;
    idx = mod(round(angle_rad / step), psk_order);
    snapped = idx * step;
end

%% ---- Gray code helpers ----

function g = gray_encode(n)
    g = bitxor(n, bitshift(n, -1));
end

function n = gray_decode(g)
    n = g;
    mask = bitshift(g, -1);
    while mask > 0
        n = bitxor(n, mask);
        mask = bitshift(mask, -1);
    end
end

function d = bin2dec_vec(bits)
% Convert binary vector [MSB ... LSB] to decimal.
    d = 0;
    for i = 1:length(bits)
        d = d * 2 + bits(i);
    end
end

function b = dec2bin_vec(d, n)
% Convert decimal to binary vector of length n.
    b = zeros(1, n);
    for i = n:-1:1
        b(i) = mod(d, 2);
        d = floor(d / 2);
    end
end

%% ---- Tag 16-QAM (square, for backscatter tag modulation) ----
% Constellation: {-3,-1,+1,+3}^2 / sqrt(10), Gray coded
% Each symbol = complex gain applied by tag (amplitude + phase)

function pts = tag16qam_constellation()
% Return 16 constellation points with Gray coding order.
% Gray-coded index -> constellation point
    vals = [-3, -1, 1, 3];
    pts = zeros(1, 16);
    idx = 1;
    for q = vals
        for ii = vals
            pts(idx) = complex(ii, q) / sqrt(10);
            idx = idx + 1;
        end
    end
end

function symbols = tag16qam_map(bits)
% Map bits to tag 16-QAM symbols. 4 bits per symbol.
    pts = tag16qam_constellation();
    bps = 4;
    n_sym = floor(length(bits) / bps);
    symbols = zeros(1, n_sym);
    for i = 1:n_sym
        b = bits((i-1)*bps+1 : i*bps);
        dec_val = bin2dec_vec(b);
        gray_dec = gray_decode(dec_val);
        symbols(i) = pts(gray_dec + 1);
    end
end

function bits = tag16qam_demap(symbols)
% Demap complex symbols to bits using tag 16-QAM constellation.
    pts = tag16qam_constellation();
    bps = 4;
    n_sym = length(symbols);
    bits = zeros(1, n_sym * bps);
    for i = 1:n_sym
        dists = abs(symbols(i) - pts);
        [~, idx] = min(dists);
        idx = idx - 1;  % 0-based
        gray_idx = gray_encode(idx);
        b = dec2bin_vec(gray_idx, bps);
        bits((i-1)*bps+1 : i*bps) = b;
    end
end

function snapped = tag16qam_snap(symbols)
% Snap complex symbols to nearest tag 16-QAM constellation point.
    pts = tag16qam_constellation();
    snapped = zeros(size(symbols));
    for i = 1:length(symbols)
        dists = abs(symbols(i) - pts);
        [~, idx] = min(dists);
        snapped(i) = pts(idx);
    end
end
