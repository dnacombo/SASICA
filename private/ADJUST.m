%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% BELOW IS ADJUST CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ADJUST() - Automatic EEG artifact Detector
% with Joint Use of Spatial and Temporal features
%
% Usage:
%   >> [art, horiz, vert, blink, disc,...
%         soglia_DV, diff_var, soglia_K, med2_K, meanK, soglia_SED, med2_SED, SED, soglia_SAD, med2_SAD, SAD, ...
%         soglia_GDSF, med2_GDSF, GDSF, soglia_V, med2_V, nuovaV, soglia_D, maxdin]=ADJUST (EEG,out)
%
% Inputs:
%   EEG        - current dataset structure or structure array (has to be epoched)
%   out        - (string) report file name
%
% Outputs:
%   art        - List of artifacted ICs
%   horiz      - List of HEM ICs
%   vert       - List of VEM ICs
%   blink      - List of EB ICs
%   disc       - List of GD ICs
%   soglia_DV  - SVD threshold
%   diff_var   - SVD feature values
%   soglia_K   - TK threshold
%   meanK      - TK feature values
%   soglia_SED - SED threshold
%   SED        - SED feature values
%   soglia_SAD - SAD threshold
%   SAD        - SAD feature values
%   soglia_GDSF- GDSF threshold
%   GDSF       - GDSF feature values
%   soglia_V   - MEV threshold
%   nuovaV     - MEV feature values
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ADJUST
% Automatic EEG artifact Detector based on the Joint Use of Spatial and Temporal features
%
% Developed 2007-2014
% Andrea Mognon (1) and Marco Buiatti (2),
% (1) Center for Mind/Brain Sciences, University of Trento, Italy
% (2) INSERM U992 - Cognitive Neuroimaging Unit, Gif sur Yvette, France
%
% Last update: 02/05/2014 by Marco Buiatti
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reference paper:
% Mognon A, Bruzzone L, Jovicich J, Buiatti M,
% ADJUST: An Automatic EEG artifact Detector based on the Joint Use of Spatial and Temporal features.
% Psychophysiology 48 (2), 229-240 (2011).
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% VERSIONS LOG
%
% 02/05/14: Modified text in Report.txt (MB).
%
% 30/03/14: Removed 'message to the user' (redundant). (MB)
%
% 22/03/14: kurtosis is replaced by kurt for compatibility if signal processing
%           toolbox is missing (MB).
%
% V2 (07 OCTOBER 2010) - by Andrea Mognon
% Added input 'nchannels' to compute_SAD and compute_SED_NOnorm;
% this is useful to differentiate the number of ICs (n) and the number of
% sensors (nchannels);
% bug reported by Guido Hesselman on October, 1 2010.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% function [art, horiz, vert, blink, disc,...
%         soglia_DV, diff_var, soglia_K, meanK, soglia_SED, SED, soglia_SAD, SAD, ...
%         soglia_GDSF, GDSF, soglia_V, nuovaV, soglia_D, maxdin]=ADJUST (EEG,out)
function [art, horiz, vert, blink, disc,...
    soglia_DV, diff_var, soglia_K, med2_K, meanK, soglia_SED, med2_SED, SED, soglia_SAD, med2_SAD, SAD, ...
    soglia_GDSF, med2_GDSF, GDSF, soglia_V, med2_V, nuovaV, soglia_D, maxdin]=ADJUST (EEG)


%% Settings

% ----------------------------------------------------
% |  Change experimental settings in this section    |
% ----------------------------------------------------

% ----------------------------------------------------
% |  Initial message to user:                        |
% ----------------------------------------------------
%
% disp(' ')
% disp('Detects Horizontal and Vertical eye movements,')
% disp('Blinks and Discontinuities in dataset:')
% disp([EEG.filename])
% disp(' ')

% ----------------------------------------------------
% |  Collect useful data from EEG structure          |
% ----------------------------------------------------

%number of ICs=size(EEG.icawinv,1);

%number of time points=size(EEG.data,2);

if length(size(EEG.data))==3

    num_epoch=size(EEG.data,3);

else

    num_epoch=1;

end

% Check the presence of ICA activations
EEG.icaact = eeg_getica(EEG);

topografie=EEG.icawinv'; %computes IC topographies

% Topographies and time courses normalization
%
% disp(' ');
% disp('Normalizing topographies...')
% disp('Scaling time courses...')

for i=1:size(EEG.icawinv,2) % number of ICs

    ScalingFactor=norm(topografie(i,:));

    topografie(i,:)=topografie(i,:)/ScalingFactor;

    if length(size(EEG.data))==3
        EEG.icaact(i,:,:)=ScalingFactor*EEG.icaact(i,:,:);
    else
        EEG.icaact(i,:)=ScalingFactor*EEG.icaact(i,:);
    end

end
%
% disp('Done.')
% disp(' ')

% Variables memorizing artifacted ICs indexes

blink=[];

horiz=[];

vert=[];

disc=[];

%% Check EEG channel position information
nopos_channels=[];
for el=1:length(EEG.chanlocs)
    if(any(isempty(EEG.chanlocs(1,el).X)&isempty(EEG.chanlocs(1,el).Y)&isempty(EEG.chanlocs(1,el).Z)&isempty(EEG.chanlocs(1,el).theta)&isempty(EEG.chanlocs(1,el).radius)))
        nopos_channels=[nopos_channels el];
    end;
end

if ~isempty(nopos_channels)
    disp(['Warning : Channels ' num2str(nopos_channels) ' have incomplete location information. They will NOT be used to compute ADJUST spatial features']);
    disp(' ');
end;

pos_channels=setdiff(EEG.icachansind,nopos_channels); %pos_channels=setdiff(1:length(EEG.chanlocs),nopos_channels); NF edit

%% Feature extraction

disp(' ')
disp('Features Extraction:')

%GDSF - General Discontinuity Spatial Feature

disp('GDSF - General Discontinuity Spatial Feature...')

GDSF = compute_GD_feat(topografie,EEG.chanlocs(1,pos_channels),size(EEG.icawinv,2));


%SED - Spatial Eye Difference

disp('SED - Spatial Eye Difference...')

[SED,medie_left,medie_right]=computeSED_NOnorm(topografie,EEG.chanlocs(1,pos_channels),size(EEG.icawinv,2));


%SAD - Spatial Average Difference

disp('SAD - Spatial Average Difference...')

[SAD,var_front,var_back,mean_front,mean_back]=computeSAD(topografie,EEG.chanlocs(1,pos_channels),size(EEG.icawinv,2));


%SVD - Spatial Variance Difference between front zone and back zone

diff_var=var_front-var_back;

%epoch dynamic range, variance and kurtosis

K=zeros(num_epoch,size(EEG.icawinv,2)); %kurtosis
Kloc=K;

Vmax=zeros(num_epoch,size(EEG.icawinv,2)); %variance

% disp('Computing variance and kurtosis of all epochs...')

for i=1:size(EEG.icawinv,2) % number of ICs

    for j=1:num_epoch
        Vmax(j,i)=var(EEG.icaact(i,:,j));
        %         Kloc(j,i)=kurtosis(EEG.icaact(i,:,j));
        K(j,i)=kurt(EEG.icaact(i,:,j));
    end
end

% check that kurt and kurtosis give the same values:
% [a,b]=max(abs(Kloc(:)-K(:)))

%TK - Temporal Kurtosis

disp('Temporal Kurtosis...')

meanK=zeros(1,size(EEG.icawinv,2));

for i=1:size(EEG.icawinv,2)
    if num_epoch>100
        meanK(1,i)=trim_and_mean(K(:,i));
    else meanK(1,i)=mean(K(:,i));
    end

end


%MEV - Maximum Epoch Variance

disp('Maximum epoch variance...')

maxvar=zeros(1,size(EEG.icawinv,2));
meanvar=zeros(1,size(EEG.icawinv,2));


for i=1:size(EEG.icawinv,2)
    if num_epoch>100
        maxvar(1,i)=trim_and_max(Vmax(:,i)');
        meanvar(1,i)=trim_and_mean(Vmax(:,i)');
    else
        maxvar(1,i)=max(Vmax(:,i));
        meanvar(1,i)=mean(Vmax(:,i));
    end
end

% MEV in reviewed formulation:

nuovaV=maxvar./meanvar;



%% Thresholds computation

disp('Computing EM thresholds...')

% soglia_K=EM(meanK);
%
% soglia_SED=EM(SED);
%
% soglia_SAD=EM(SAD);
%
% soglia_GDSF=EM(GDSF);
%
% soglia_V=EM(nuovaV);
[soglia_K,med1_K,med2_K]=EM(meanK);

[soglia_SED,med1_SED,med2_SED]=EM(SED);

[soglia_SAD,med1_SAD,med2_SAD]=EM(SAD);

[soglia_GDSF,med1_GDSF,med2_GDSF]=EM(GDSF);

[soglia_V,med1_V,med2_V]=EM(nuovaV);



%% Horizontal eye movements (HEM)


horiz=intersect(intersect(find(SED>=soglia_SED),find(medie_left.*medie_right<0)),...
    (find(nuovaV>=soglia_V)));




%% Vertical eye movements (VEM)



vert=intersect(intersect(find(SAD>=soglia_SAD),find(medie_left.*medie_right>0)),...
    intersect(find(diff_var>0),find(nuovaV>=soglia_V)));



%% Eye Blink (EB)


blink=intersect ( intersect( find(SAD>=soglia_SAD),find(medie_left.*medie_right>0) ) ,...
    intersect ( find(meanK>=soglia_K),find(diff_var>0) ));



%% Generic Discontinuities (GD)



disc=intersect(find(GDSF>=soglia_GDSF),find(nuovaV>=soglia_V));


%compute output variable
art = nonzeros( union (union(blink,horiz) , union(vert,disc)) )'; %artifact ICs

% these three are old outputs which are no more necessary in latest ADJUST version.
soglia_D=0;
soglia_DV=0;
maxdin=zeros(1,size(EEG.icawinv,2));

return

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


% computeSED_NOnorm() - Computes Spatial Eye Difference feature
% without normalization
%
% Usage:
%   >> [out,medie_left,medie_right]=computeSED_NOnorm(topog,chanlocs,n);
%
% Inputs:
%   topog      - topographies vector
%   chanlocs   - EEG.chanlocs struct
%   n          - number of ICs
%   nchannels  - number of channels
%
% Outputs:
%   out        - SED values
%   medie_left - Left Eye area average values
%   medie_right- Right Eye area average values
%
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

function [out,medie_left,medie_right]=computeSED_NOnorm(topog,chanlocs,n)

nchannels=length(chanlocs);

%% Define scalp zones

% Find electrodes in Left Eye area (LE)
dimleft=0; %number of LE electrodes
index1=zeros(1,nchannels); %indexes of LE electrodes

for k=1:nchannels
    if (-61<chanlocs(1,k).theta) && (chanlocs(1,k).theta<-35) && (chanlocs(1,k).radius>0.30) %electrodes are in LE
        dimleft=dimleft+1; %count electrodes
        index1(1,dimleft)=k;
    end
end

% Find electrodes in Right Eye area (RE)
dimright=0; %number of RE electrodes
index2=zeros(1,nchannels); %indexes of RE electrodes
for g=1:nchannels
    if (34<chanlocs(1,g).theta) && (chanlocs(1,g).theta<61) && (chanlocs(1,g).radius>0.30) %electrodes are in RE
        dimright=dimright+1; %count electrodes
        index2(1,dimright)=g;
    end
end

% Find electrodes in Posterior Area (PA)
dimback=0;
index3=zeros(1,nchannels);
for h=1:nchannels
    if (abs(chanlocs(1,h).theta)>110)
        dimback=dimback+1;
        index3(1,dimback)=h;
    end
end

if dimleft*dimright*dimback==0
    disp('ERROR: no channels included in some scalp areas.')
    disp('Check channels distribution and/or change scalp areas definitions in computeSAD.m and computeSED_NOnorm.m')
    disp('ADJUST session aborted.')
    return
end

%% Outputs

out=zeros(1,n); %memorizes SED
medie_left=zeros(1,n); %memorizes LE mean value
medie_right=zeros(1,n); %memorizes RE mean value

%% Output computation

for i=1:n  % for each topography
    %create LE electrodes vector
    left=zeros(1,dimleft);
    for h=1:dimleft
        left(1,h)=topog(i,index1(1,h));
    end

    %create RE electrodes vector
    right=zeros(1,dimright);
    for h=1:dimright
        right(1,h)=topog(i,index2(1,h));
    end

    %create PA electrodes vector
    back=zeros(1,dimback);
    for h=1:dimback
        back(1,h)=topog(i,index3(1,h));
    end



    %compute features
    out1=abs(mean(left)-mean(right));
    out2=var(back);
    out(1,i)=out1; % SED not notmalized
    medie_left(1,i)=mean(left);
    medie_right(1,i)=mean(right);


end


% computeSAD() - Computes Spatial Average Difference feature
%
% Usage:
%   >> [rapp,var_front,var_back,mean_front,mean_back]=computeSAD(topog,chanlocs,n);
%
% Inputs:
%   topog      - topographies vector
%   chanlocs   - EEG.chanlocs struct
%   n          - number of ICs
%   nchannels  - number of channels
%
% Outputs:
%   rapp       - SAD values
%   var_front  - Frontal Area variance values
%   var_back   - Posterior Area variance values
%   mean_front - Frontal Area average values
%   mean_back  - Posterior Area average values
%
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


function [rapp,var_front,var_back,mean_front,mean_back]=computeSAD(topog,chanlocs,n)

nchannels=length(chanlocs);

%% Define scalp zones

% Find electrodes in Frontal Area (FA)
dimfront=0; %number of FA electrodes
index1=zeros(1,nchannels); %indexes of FA electrodes

for k=1:nchannels
    if (abs(chanlocs(1,k).theta)<60) && (chanlocs(1,k).radius>0.40) %electrodes are in FA
        dimfront=dimfront+1; %count electrodes
        index1(1,dimfront)=k;
    end
end

% Find electrodes in Posterior Area (PA)
dimback=0;
index3=zeros(1,nchannels);
for h=1:nchannels
    if (abs(chanlocs(1,h).theta)>110)
        dimback=dimback+1;
        index3(1,dimback)=h;
    end
end

if dimfront*dimback==0
    disp('ERROR: no channels included in some scalp areas.')
    disp('Check channels distribution and/or change scalp areas definitions in computeSAD.m and computeSED_NOnorm.m')
    disp('ADJUST session aborted.')
    return
end

%% Outputs

rapp=zeros(1,n); % SAD
mean_front=zeros(1,n); % FA electrodes mean value
mean_back=zeros(1,n); % PA electrodes mean value
var_front=zeros(1,n); % FA electrodes variance value
var_back=zeros(1,n); % PA electrodes variance value

%% Output computation

for i=1:n % for each topography

    %create FA electrodes vector
    front=zeros(1,dimfront);
    for h=1:dimfront
        front(1,h)=topog(i,index1(1,h));
    end

    %create PA electrodes vector
    back=zeros(1,dimback);
    for h=1:dimback
        back(1,h)=topog(i,index3(1,h));
    end



    %compute features

    rapp(1,i)=abs(mean(front))-abs(mean(back)); % SAD
    mean_front(1,i)=mean(front);
    mean_back(1,i)=mean(back);
    var_back(1,i)=var(back);
    var_front(1,i)=var(front);

end


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
    [last,med1,med2,var1,var2,prior1,prior2] = rep2struct(NaN);
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prob=Bayes(med,var,point)
if var==0
    prob=1;
else
    prob=((1/(sqrt(2*pi*var)))*exp((-1)*((point-med)^2)/(2*var)));
end





% trim_and_max() - Computes maximum value from vector 'vettore'
% after removing the top 1% of the values
% (to be outlier resistant)
%
% Usage:
%   >> valore=trim_and_max(vettore);
%
% Inputs:
%   vettore    - row vector
%
% Outputs:
%   valore     - result
%
%
% Author: Andrea Mognon, Center for Mind/Brain Sciences, University of
% Trento, 2009

% Motivation taken from the following comment to our paper:
% "On page 11 the authors motivate the use of the max5 function when computing
% Maximum Epoch Variance because the simple maximum would be too sensitive
% to spurious outliers. This is a good concern, however the max5 function would
% still be sensitive to spurious outliers for very large data sets. In other words, if
% the data set is large enough, one will be very likely to record more than five
% outliers. The authors should use a trimmed max function that computes the
% simple maximum after the top say .1% of the values have been removed from
% consideration. This rejection criteria scales appropriately with the size of the data
% set."

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


function valore=trim_and_max(vettore)


dim=floor(.01*size(vettore,2)); % = 1% of vector length

tmp=sort(vettore);
valore= tmp(length(vettore)-dim);



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





