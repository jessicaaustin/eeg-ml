function [out] = gen_epoch_features(varargin)
% Start Cache
[success, INTERN_cache_desc, varargin] = cache_enter(varargin);
if (success) == 1, out = INTERN_cache_desc; return; end;

data = varargin{1};
bands = varargin{2};
sampling_rate = varargin{3};

for t = 1:length(data.start_time)
    % TODO reconsider this below line
    for e = 1:numel(data.epochs{t})
        epoch = data.epochs{t}(e);
        if isempty(epoch.signal)
            continue
        end
        % GENERATE OSCILLATION BANDS
        
        features = gen_oscillation_features(bands, epoch, sampling_rate);
        data.epochs{t}(e).features = features;
    end
end

out = data;
% End Cache
out = cache_exit(INTERN_cache_desc, out);
end

function [features] = gen_oscillation_features(bands, channel, sampling_rate)
        % gen powerspec
        x = channel.signal;

        fs = sampling_rate;        % Sample frequency (Hz)
        nfft = pow2(nextpow2(fs)); % Transform length

        f = fs / 2 * linspace(0, 1, nfft/2+1);
        F = fft(x', nfft, 2);
        f = f(2:end);
        F = F(:,2:floor(nfft/2)+1);

        y = F .* conj(F) / nfft;
        %y = 2 * abs(F) / fs;
        
        powerspec = y;

        % gen oscillation features
        for i = 1:length(bands) - 1,
            low = bands(i);
            high = bands(i+1)-1;
            features.(sprintf('band_%d_%d', low, high)) = nansum(powerspec(:, low:high), 2) ./ ...
                nansum(powerspec(:, 1:sampling_rate/2), 2);
        end % i
end