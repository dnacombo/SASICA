function rm_frompath(what)

p = [path pathsep];

ps = regexp(p,['(.*?)' pathsep],'tokens');
for i = 1:numel(ps)
    ps{i} = ps{i}{1};
end
ps = ps(regexpcell(ps,what,'inv'));
p = [];
for i= 1:numel(ps)
    p = [p pathsep ps{i}];
end
path(p);


