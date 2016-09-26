function [varargout] = rep2struct(varargin)

% [s.target] = rep2struct(dat)
% replicate the value dat into each element of structure s in field target.
% if dat has same nb or elements as s, each element of dat goes into one
% element of s. if dat is more dimensional and doesn't have the same number
% of elements as s, and has same size along dimension 1, then pass each
% slice into s.target.

if numel(varargin) == 1
    dat = varargin{1};
    if numel(dat) == nargout
        for i = 1:nargout
            varargout{i} = dat(i);
        end
    elseif size(dat,1) == nargout
        for i = 1:nargout
            varargout{i} = dat(i,:);
        end
    else
        for i = 1:nargout
            varargout{i} = dat;
        end
    end
elseif numel(varargin) == nargout
    for i = 1:nargout
        varargout{i} = varargin{i};
    end
else
    error('Wrong number of arguments');
end

