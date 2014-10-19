function applyColorLineStyleIndependently(h)

c = num2cell(get(0,'DefaultAxesColorOrder'), 2);
l = cellstr(get(0,'DefaultAxesLineStyleOrder'));
set(h, {'Color'}, c(rem((1:numel(h))-1,numel(c))+1), ...
    {'LineStyle'}, l(rem((1:numel(h))-1,numel(l))+1))
