%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EM - ADJUST package
%
% Performs automatic threshold on the digital numbers
% of the input vector 'vec'; based on Expectation - Maximization algorithm

% Reference paper:
% Bruzzone, L., Prieto, D.F., 2000. Automatic analysis of the difference image
% for unsupervised change detection.
% IEEE Trans. Geosci. Remote Sensing 38, 1171:1182

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Usage:
%   >> [last,med1,med2,var1,var2,prior1,prior2]=EM(vec);
%
% Input: vec (row vector, to be thresholded)
%
% Outputs: last (threshold value)
%          med1,med2 (mean values of the Gaussian-distributed classes 1,2)
%          var1,var2 (variance of the Gaussian-distributed classes 1,2)
%          prior1,prior2 (prior probabilities of the Gaussian-distributed classes 1,2)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


function [last,med1,med2,var1,var2,prior1,prior2]=EM(vec)

if size(vec,2)>1
    len=size(vec,2); %number of elements
else
    vec=vec';
    len=size(vec,2);
end

c_FA=1; % False Alarm cost
c_MA=1; % Missed Alarm cost

med=mean(vec);
standard=std(vec);
mediana=(max(vec)+min(vec))/2;

alpha1=0.01*(max(vec)-mediana); % initialization parameter/ righthand side
alpha2=0.01*(mediana-min(vec)); % initialization parameter/ lefthand side

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPECTATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

train1=[]; % Expectation of class 1
train2=[];
train=[]; % Expectation of 'unlabeled' samples

for i=1:(len)
    if (vec(i)<(mediana-alpha2))
        train2=[train2 vec(i)];
    elseif (vec(i)>(mediana+alpha1))
        train1=[train1 vec(i)];
    else
        train=[train vec(i)];
    end
end

n1=length(train1);
n2=length(train2);

med1=mean(train1);
med2=mean(train2);
prior1=n1/(n1+n2);
prior2=n2/(n1+n2);
var1=var(train1);
var2=var(train2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAXIMIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count=0;
dif_med_1=1; % difference between current and previous mean
dif_med_2=1;
dif_var_1=1; % difference between current and previous variance
dif_var_2=1;
dif_prior_1=1; % difference between current and previous prior
dif_prior_2=1;
stop=0.0001;

while((dif_med_1>stop)&&(dif_med_2>stop)&&(dif_var_1>stop)&&(dif_var_2>stop)&&(dif_prior_1>stop)&&(dif_prior_2>stop))

    count=count+1;

    med1_old=med1;
    med2_old=med2;
    var1_old=var1;
    var2_old=var2;
    prior1_old=prior1;
    prior2_old=prior2;
    prior1_i=[];
    prior2_i=[];

    % FOLLOWING FORMULATION IS ACCORDING TO REFERENCE PAPER:

    for i=1:len
        prior1_i=[prior1_i prior1_old*Bayes(med1_old,var1_old,vec(i))/...
            (prior1_old*Bayes(med1_old,var1_old,vec(i))+prior2_old*Bayes(med2_old,var2_old,vec(i)))];
        prior2_i=[prior2_i prior2_old*Bayes(med2_old,var2_old,vec(i))/...
            (prior1_old*Bayes(med1_old,var1_old,vec(i))+prior2_old*Bayes(med2_old,var2_old,vec(i)))];
    end


    prior1=sum(prior1_i)/len;
    prior2=sum(prior2_i)/len;
    med1=sum(prior1_i.*vec)/(prior1*len);
    med2=sum(prior2_i.*vec)/(prior2*len);
    var1=sum(prior1_i.*((vec-med1_old).^2))/(prior1*len);
    var2=sum(prior2_i.*((vec-med2_old).^2))/(prior2*len);

    dif_med_1=abs(med1-med1_old);
    dif_med_2=abs(med2-med2_old);
    dif_var_1=abs(var1-var1_old);
    dif_var_2=abs(var2-var2_old);
    dif_prior_1=abs(prior1-prior1_old);
    dif_prior_2=abs(prior2-prior2_old);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THRESHOLDING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=c_MA/c_FA;
a=(var1-var2)/2;
b= ((var2*med1)-(var1*med2));
c=(log((k*prior1*sqrt(var2))/(prior2*sqrt(var1)))*(var2*var1))+(((((med2)^2)*var1)-(((med1)^2)*var2))/2);
rad=(b^2)-(4*a*c);
if rad<0
    disp('Negative Discriminant!');
    return;
end

soglia1=(-b+sqrt(rad))/(2*a);
soglia2=(-b-sqrt(rad))/(2*a);

if ((soglia1<med2)||(soglia1>med1))
    last=soglia2;
else
    last=soglia1;
end

if isnan(last) % TO PREVENT CRASHES
    last=mediana;
end

return



