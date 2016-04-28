function [M100, idx_clab_desired] = get_M100_ADE(clab_desired)
% [M100, idx_clab_desired] = get_M100_ADEC(clab_desired)
%
% IN  clab_desired - channel setup for which M100 should be calculated
% OUT M100
%     idx_clab_desired
% M100 is the matrix such that  feature = norm(M100*ica_pattern(idx_clab_desired), 'fro')
%
% (c) Stefan Haufe

lambda = 100;

load inv_matrix_icbm152; %L (forward matrix 115 x 2124 x 3), clab (channel labels)

[cl_ ia idx_clab_desired] = intersect(clab, clab_desired);
F = L(ia, :, :); %forward matrix for desired channels labels
[n_channels m foo] = size(F);  %m = 2124, number of dipole locations
F = reshape(F, n_channels, 3*m);

%H - matrix that centralizes the pattern, i.e. mean(H*pattern) = 0
H = eye(n_channels) -  ones(n_channels, n_channels)./ n_channels;
%W - inverse of the depth compensation matrix Lambda
W = sloreta_invweights(L);

L = H*F*W;

%We have inv(L'L +lambda eye(size(L'*L))* L' = L'*inv(L*L' + lambda
%eye(size(L*L')), which is easier to calculate as number of dimensions is
%much smaller

%calulate the inverse of L*L' + lambda * eye(size(L*L')
[U D] = eig(L*L');
d = diag(D);
di = d+lambda;
di = 1./di;
di(d < 1e-10) = 0;
inv1 = U*diag(di)*U';  %inv1 = inv(L*L' + lambda *eye(size(L*L'))

%get M100
M100 = L'*inv1*H;



