function [out] = smooth(varargin)

% Start Cache
[success, INTERN_cache_desc, varargin] = cache_enter(varargin);
if (success) == 1, out = INTERN_cache_desc; return; end;

data = varargin{1};
denoise_N = varargin{2};

for d = 1:length(data)
    [thr, sorh, keepapp] = ddencmp('den', 'wv', data{d}.rawwave.signal); % find the threshholding value
    data{d}.rawwave.signal = wdencmp('gbl', data{d}.rawwave.signal, 'db3', denoise_N, thr, sorh, keepapp); % de-noising using the threshholding value
end

out = data;

% End Cache
out = cache_exit(INTERN_cache_desc, out);
end
