function [i_te, i_tr] = findconvertedlabels(pos_3d, chanlocs)
% IN  pos_3d  - 3d-positions of training channel labels
%     chanlocs - EEG.chanlocs structure of data to be classified

%compute spherical coordinates theta and phi for the training channel
%label
[theta, phi, r] = cart2sph(pos_3d(1,:),pos_3d(2,:), pos_3d(3,:));
theta = theta - pi/2;
theta(theta < -pi) = theta(theta < -pi) + 2*pi;
theta = theta*180/pi;
phi = phi * 180/pi;
theta(find(pos_3d(1,:) == 0 & pos_3d(2,:) == 0)) = 0; %exception for Cz


clab_common = {};
i_te = [];
i_tr = [];

%For each channel in EEG.chanlocs, try to find matching channel in
%training data
for chan = 1:length(chanlocs)
    if not(isempty(chanlocs(chan).sph_phi))
        idx = find((theta <= chanlocs(chan).sph_theta + 6) ...
            & (theta >= chanlocs(chan).sph_theta - 6) ...
            & (phi <= chanlocs(chan).sph_phi + 6) ...
            & (phi >= chanlocs(chan).sph_phi - 6));
        if not(isempty(idx))
            i_tr = [i_tr, idx(1)];
            i_te = [i_te, chan];
        end
    end
end
