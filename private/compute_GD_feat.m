% compute_GD_feat() - Computes Generic Discontinuity spatial feature
%
% Usage:
%   >> res = compute_GD_feat(topografie,canali,num_componenti);
%
% Inputs:
%   topografie - topographies vector
%   canali     - EEG.chanlocs struct
%   num_componenti  - number of components
%
% Outputs:
%   res       - GDSF values

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


function res = compute_GD_feat(topografie,canali,num_componenti)

% Computes GDSF, discontinuity spatial feature
% topografie is the topography weights matrix
% canali is the structure EEG.chanlocs
% num_componenti is the number of ICs
% res is GDSF values

xpos=[canali.X];ypos=[canali.Y];zpos=[canali.Z];
pos=[xpos',ypos',zpos'];

res=zeros(1,num_componenti);

for ic=1:num_componenti

    % consider the vector topografie(ic,:)

    aux=[];

    for el=1:length(canali)-1

        P=pos(el,:); %position of current electrode
        d=pos-repmat(P,length(canali),1);
        %d=pos-repmat(P,62,1);
        dist=sqrt(sum((d.*d),2));

        [y,I]=sort(dist);
        repchas=I(2:11); % list of 10 nearest channels to el
        weightchas=exp(-y(2:11)); % respective weights, computed wrt distance

        aux=[aux abs(topografie(ic,el)-mean(weightchas.*topografie(ic,repchas)'))];
        % difference between el and the average of 10 neighbors
        % weighted according to weightchas
    end

    res(ic)=max(aux);

end


