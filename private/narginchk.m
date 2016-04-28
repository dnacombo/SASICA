function narginchk(min,max)

n = evalin('caller','nargin');
if  n < min || n > max
    error('number of arguments')
end
