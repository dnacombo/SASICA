function h = vline(x,varargin)

% h = vline(x,varargin)
% add vertical line(s) on the current axes at x
% all varargin arguments are passed to plot...

x = x(:);
ho = ishold;
hold on
h = plot([x x]',repmat(ylim,numel(x),1)',varargin{:});
if not(ho)
    hold off
end
if nargout == 0
    clear h
end
