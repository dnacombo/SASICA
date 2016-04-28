function thresh = readauto(thresh,dat,comp)
% if thresh starts with 'auto'
% compute auto threshold as mean(dat) +/- N std(dat)
% with N read in the string thresh = 'auto N'
% if not, use thresh as a value
if isstr(thresh) && strncmp(thresh,'auto',4)
    if numel(thresh) > 4
        threshsigma = str2num(thresh(5:end));
    else
        threshsigma = 2;
    end
    thresh = eval(['mean(dat,2)' comp 'threshsigma * std(dat,[],2)']);
end



