%***************************
%*  Segment explorer for partitioning purposes
%   Desgined for segment training
%*  2013,4,11
%*  References :
%*  Mittal10, 
%*  
%***************************
%
%******************************


mkdir tmp_data;

clear all; clc; close all;
addpath src
addpath ../ReadData3D_version1k/mha
inDir = 'Training_seg/';

%% Flags
visualDebug=false;
%% Training mode (decision criters)
% 0-> binary 0=healthy,soft|1=calc,mix
% 1-> binary 0=healthy,soft,mix|1=calc
% 2-> binary 0 healthy|1=calc,mix,narr>50
% 3-> multi 0=healthy|1=calc,mix|2=narr>50
% 4-> binary 0=healthy|1=grade narrowing>50
trainingMode=0;

yTrain=[];
numseg=0;
ssCounter=0;
ssList=[];
hsCounter=0;
hsList=[];
L=5;
tic
%% Segment lecture
MHD= dir(fullfile(inDir,'*.mhd'));
for j =1:numel(MHD),
    
    %% Files lecture
    refFilename=fullfile(inDir,[MHD(j).name(1:end-4) '.txt'])
    reference=load(refFilename);
    [dims]=size(reference)
    numseg=numseg+1;
    %key=input('key input');
    
    sickSegFlag=0;
    segmentyTrain=[];
    for z0=round(L/2):1:dims(1)-round(L/2)+1     % loop through CPR

        if trainingMode==0
            if reference(z0,6)>1.0
                segmentyTrain=[segmentyTrain; 1];
                sickSegFlag=1;
            else
                segmentyTrain=[segmentyTrain; 0];
            end
        elseif trainingMode==1
            if reference(z0,6)==2.0
                segmentyTrain=[segmentyTrain; 1];
                sickSegFlag=1;
            else
                segmentyTrain=[segmentyTrain; 0];
            end
            
        elseif trainingMode==2
            if reference(z0,6)>1.0
                segmentyTrain=[segmentyTrain; 1];
                sickSegFlag=1;
            elseif reference(z0,7)>1
                segmentyTrain=[segmentyTrain; 1];
                sickSegFlag=1;
            else
                segmentyTrain=[segmentyTrain; 0];
            end
        elseif trainingMode==3
            if reference(z0,6)>1.0
                segmentyTrain=[segmentyTrain; 1];
                sickSegFlag=1;
            elseif reference(z0,7)>1
                segmentyTrain=[segmentyTrain; 2];
                sickSegFlag=1;
            else
                segmentyTrain=[segmentyTrain; 0];
            end
        elseif trainingMode==4
            if reference(z0,7)>1.0
                segmentyTrain=[segmentyTrain; 1];
                sickSegFlag=1;
            else
                segmentyTrain=[segmentyTrain; 0];
            end
        else
            segmentyTrain=[segmentyTrain; reference(z0,6)];
        end

    end   %rof process vessel
    
    if(sickSegFlag==1)
        ssCounter=ssCounter+1;
        ssList=[MHD(j); ssList];
    else
        hsCounter=hsCounter+1;
        hsList=[MHD(j); hsList];
        %hsList{hsCounter}={MHD(j).name(1:end-4)};
    end
    yTrain=[yTrain; segmentyTrain];
end
toc

% print healthy Segment list
for i=1:numel(hsList)
   hsList(i).name(1:end-4)
end

% print sick Segment list
for i=1:numel(ssList)
   ssList(i).name(1:end-4)
end

[trainingMode numseg ssCounter 14]
%     2   245    63    14
fname='tm0_segments.txt';
fid=fopen(fname,'W');
if fid~=-1
    for i=1:numel(ssList)
       fprintf(fid,'%s\n',ssList(i).name(1:end-4));
    end
end
fclose(fid);

% Random permuation of positive and negative data points
p = randperm(numel(ssList));
n = randperm(numel(hsList));

% 80-20 split for training and test
tstpf = p(1:round(numel(ssList)/10));
tstnf = n(1:round(numel(hsList)/10));
trpf = setdiff(p, tstpf);
trnf = setdiff(n, tstnf);
% 
% c = randperm(numel(ssList));
%     s = randperm(numel(hsList));
%     h = randperm(numel(hsList));
% 
%     % 66-33 split for training and test
%     tstsf = s(1:round(row_soft/5));
%     tstcf = c(1:round(row_calc/5));
%     tsthf = h(1:round(row_heal/5));
%     trsf = setdiff(s, tstsf);
%     trcf = setdiff(c, tstcf);
%     trhf = setdiff(h, tsthf);

save 'tmp_data/TM0_90_oob_list' tstpf tstnf trpf trnf ssList hsList
