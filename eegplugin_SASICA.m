function vers = eegplugin_SASICA(fig,try_strings,catch_strings)

vers = 'SASICA 1.3.8';

if nargin == 0
    return
end

toolsmenu = findobj(fig,'Tag','tools');

u = uimenu( toolsmenu,'separator','on', 'label', 'SASICA', 'callback',[ 'SASICA;' ], 'userdata', 'startup:off;study:on;chanloc:on;ica:on'); 
