function vers = eegplugin_SASICA(fig,try_strings,catch_strings)

vers = 'SASICA_1.3.7';

if nargin == 0
    return
end

toolsmenu = findobj(fig,'Tag','tools');

u = uimenu( toolsmenu,'separator','on', 'label', 'SASICA', 'callback',[ 'SASICA;' ], 'userdata', 'startup:off;study:on;chanloc:on;ica:on'); 
