function [out] = balance_data(varargin)
    I = varargin{1};
    data = varargin{2};
    balance = varargin{3};
    % TODO don't read conds this way
    H = data.task_H;
    list = H; enum = 1; enum_list;
    balance_I = balance_class(data.task_M(I, COND), find(data.task_M(I, COND)), balance);
    out = I(balance_I);
end