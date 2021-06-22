function pos = getpointerpos(h,evt)

if not(exist('h','var')) || ~strcmp(get(h,'type'),'axes')
    h = gca;
end

prev = get(0,'units');
set(0,'units','pixels');
pointp = get(0,'PointerLocation');
set(0,'units',prev);

prev = get(gcf,'units');
set(gcf,'units','pixels');
fp = get(gcf,'Position');
set(gcf,'units',prev);

prev = get(h,'units');
set(h,'units','pixels');
axp = get(h,'Position');
set(h,'units',prev);

pointp = pointp - fp(1:2) - axp(1:2);

xl = xlim;yl = ylim;
pos = (pointp + [1 2]  )./ [axp(3) axp(4)] .* [diff(xl) diff(yl)] + [xl(1) yl(1)];

% set(h,'xlimmode','manual','ylimmode','manual')
% hold on
% scatter(pos(1),pos(2))

