function [headers, vals] = describe(expt),

    [headers, vals] = describe_fields(expt, {'feature_headers'});
    [fheaders, fvals] = describe_features(expt.feature_headers);
    headers = sprintf('%s\t%s', fheaders, headers);
    vals = sprintf('%s\t%s', fvals, vals);
    
end

function [header, vals] = describe_fields(expt, exceptions),
    headers = {};
    val_array = {};
    fields = fieldnames(expt);
    for i = 1:numel(fields),
        field = fields{i};
        if (ismember(field, exceptions)),
            continue;
        end
        val_array{end + 1} = val2str(eval(sprintf('expt.%s', field)));
        headers{end + 1} = field;
    end
    vals = [sprintf('%s\t',val_array{1:end-1}),val_array{end}];
    header = [sprintf('%s\t',headers{1:end-1}),headers{end}];
end

function [text] = val2str(val)
    text = '';
    if (ischar(val)),
        text = sprintf('%s', val);
    elseif (isfloat(val)),
        if (mod(val, 1)),
            text = sprintf('%.4f', val);
        else
            text = sprintf('%d', val);
        end
    elseif (iscell(val)),
        if length(val) == 0,
            return;
        end
        text = val2str(val{1});
        for i = 2:length(val),
            text = sprintf('%s, %s', text, val2str(val{i}));
        end
    else
        text = sprintf('variable could not be printed');
    end
end

function [headers, vals] = describe_features(features),
    [headers, vals] = describe_fields(features, {'segments'});
    
    starts_array = {};
    ends_array = {};
    for s = 1:length(features.segments),
        starts_array{end + 1} = features.segments(s).start_time;
        ends_array{end + 1} = features.segments(s).end_time;
    end
    
    starts = starts_array{1};
    ends = ends_array{1};
    for i = 2:length(starts_array),
        starts = sprintf('%s, %s', starts, starts_array{i});
        ends = sprintf('%s, %s', ends, ends_array{i});
    end
    headers = sprintf('%s\t%s\t%s', headers, 'feature_starts', 'feature_ends');
    vals = sprintf('%s\t%s\t%s', vals, starts, ends);
end