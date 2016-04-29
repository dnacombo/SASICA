function h = hline(y,varargin)

% h = hline(y,varargin)
% add horizontal line(s) on the current axes at y
% all varargin arguments are passed to plot...

y = y(:);
ho = ishold;
hold on
h = plot(repmat(xlim,numel(y),1)',[y y]',varargin{:});
if not(ho)
    hold off
end
if nargout == 0
    clear h
end


