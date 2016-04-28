%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prob=Bayes(med,var,point)
if var==0
    prob=1;
else
    prob=((1/(sqrt(2*pi*var)))*exp((-1)*((point-med)^2)/(2*var)));
end





