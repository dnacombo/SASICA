function alpha = msk2alpha(msk,m,M)

% replace in msk all 0 by m and all 1 by M
if not(islogical(msk)) || isequal(unique(msk),[0,1])
    error('Cannot deal with non binary msk')
end
alpha = msk;
alpha(alpha==0) = m;
alpha(alpha==1) = M;

