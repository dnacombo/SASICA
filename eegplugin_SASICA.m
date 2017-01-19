function vers = eegplugin_SASICA(fig,try_strings,catch_strings)

vers = 'SASICA_1.3.4';

if nargin == 0
    return
end

toolsmenu = findobj(fig,'Tag','tools');
abovemenupos = get(findobj(toolsmenu,'Label','Remove components'),'position');


uimenu( toolsmenu,'position',abovemenupos+1,'separator','on', 'label', 'SASICA', 'callback',...
    [ 'SASICA(EEG);' ]); 