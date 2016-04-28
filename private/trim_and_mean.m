% trim_and_mean() - Computes average value from vector 'vettore'
% after removing the top .1% of the values
% (to be outlier resistant)
%
% Usage:
%   >> valore=trim_and_mean(vettore);
%
% Inputs:
%   vettore    - row vector
%
% Outputs:
%   valore     - result
%
% Copyright (C) 2009-2014 Andrea Mognon (1) and Marco Buiatti (2),
% (1) Center for Mind/Brain Sciences, University of Trento, Italy
% (2) INSERM U992 - Cognitive Neuroimaging Unit, Gif sur Yvette, France
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


function valore=trim_and_mean(vettore)


dim=floor(.01*size(vettore,2)); % = 1% of vector length

tmp=sort(vettore);
valore= mean (tmp(1:(length(vettore)-dim)));





