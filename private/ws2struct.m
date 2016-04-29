function [sortie] = ws2struct(varargin)

% [s] = ws2struct(varargin)
%
% Description : This function returns a structure containing variables
% of the current workspace.
% __________________________________
% Inputs :
%   re (string optional) :  a regular expression matching the variables to
%                           be returned.
% Outputs :
%   s (structure array) :   a structure containing all variables of the
%                           calling workspace. If re input is specified,
%                           only variables matching re are returned.
% _____________________________________
% See also : struct2ws ; regexp
%
% Maximilien Chaumon v1.0 02/2007


if not(isempty(varargin))
    re = varargin{1};
else
    re = '.*';
end

vars = evalin('caller','who');
vmatch = regexp(vars,re);
varsmatch = [];
for i = 1:length(vmatch)
    if isempty(vmatch{i}) || not(vmatch{i} == 1)
        continue
    end
    varsmatch{end+1} = vars{i};
end

for i = 1:length(varsmatch)
    dat{i} = evalin('caller',varsmatch{i});
end

sortie = cell2struct(dat,varsmatch,2);

