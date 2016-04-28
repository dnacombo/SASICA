function setctxt(hfig,EEG,cfg)
COLREJ = '[1 0.6 0.6]';
COLACC = '[0.75 1 0.75]';
buttons = findobj(hfig,'-regexp','tag','^comp\d{1,3}$');
buttonnums = regexp(get(buttons,'tag'),'comp(\d{1,3})','tokens');
if numel(buttonnums)>1
    buttonnums = cellfun(@(x)(str2num(x{1}{1})),buttonnums);
else
    buttonnums = str2num(buttonnums{1}{1});
end
for i = 1:numel(buttonnums)
    hcmenu = uicontextmenu;

    if ~isempty(EEG.reject.gcompreject)
        status = EEG.reject.gcompreject(buttonnums(i));
    else
        status = 0;
    end;

    hcb1 = ['EEG.reject.gcompreject(' num2str(buttonnums(i)) ') = ~EEG.reject.gcompreject(' num2str(buttonnums(i)) ');'...
        'set(gco,''backgroundcolor'',fastif(EEG.reject.gcompreject(' num2str(buttonnums(i)) '), ' COLREJ ',' COLACC '));'...
        'set(findobj(''tag'',''ctxt' num2str(buttonnums(i)) '''), ''Label'',fastif(EEG.reject.gcompreject(' num2str(buttonnums(i)) '),''ACCEPT'',''REJECT''));' ];
    uimenu(hcmenu, 'Label', fastif(status,'ACCEPT','REJECT'), 'Callback', hcb1,'tag',['ctxt' num2str(buttonnums(i))]);

    mycb = strrep(get(buttons(i),'Callback'),'''','''''');
    mycb = regexprep(mycb,'pop_prop','eeg_SASICA(EEG,''pop_prop');
    mycb = [mycb ''');'];
    set(buttons(i),'CallBack',mycb)
    set(buttons(i),'uicontextmenu',hcmenu)
end

