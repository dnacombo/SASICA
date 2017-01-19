function s = vararg2struct(v,tag)

% s = vararg2struct(v,tag)
% 
% translate a sequence of varargin 'name', value pairs into a structure.
% substructure fields can be defined by using underscores (if tag is
% provided, another character can be used)
% 
% ex:   v = {'name','toto','size',55,'hair_style','cool','hair_color','blue'}
%       s = vararg2struct(v)
% s = 
%     name: 'toto'
%     size: 55
%     hair: [1x1 struct]
% s.hair
% ans = 
%     style: 'cool'
%     color: 'blue'
% 

if not(exist('tag','var'))
    tag = '_';
end
s = struct;
f = regexp(v(1:2:end),['[^'  regexptranslate('escape',tag) ']*'],'match');
for i_f = 1:numel(f)
    str = 's';
    for i_ff = 1:numel(f{i_f})
        str = [str '.' f{i_f}{i_ff}];
    end
    str = [str ' = v{i_f*2};'];
    eval(str);
end

