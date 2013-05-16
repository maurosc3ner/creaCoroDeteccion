%***************************
%*  3D Model Training Binary class:
%*  I,Gr,Gt,Radon
%   Longitudinal 
%
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
inDir = 'Training_vessels/';

%busqueda de carpetas
DT= dir(fullfile(inDir,'dt*'));

%% Flags
visualDebug=false;
%% Training mode 
% 0-> binary 0=healthy,soft|1=calc,mix
% 1-> binary 0=healthy,soft,mix|1=calc
trainingMode=0;

%% Cylinder mask creation



t = pi/4:pi/4:2*pi;
theta = 0:180;
trainData=[];
yTrain=[];
tic
%% dataset lecture
for j =1:numel(DT),
    
    %% Files lecture
    inDT=strcat(inDir,DT(j).name);
    %busqueda de vasos
    MHD= dir(fullfile(inDT,'*.mhd'));
    
    vesselTrainData=[];
    vesselyTrain=[];
    
    for vessel_i=1:numel(MHD),
       
        refFilename=fullfile(inDT,[MHD(vessel_i).name(1:end-4) '.txt'])
        
        cprFilename=fullfile(inDT,MHD(vessel_i).name);
        reference=load(refFilename);
        info = mha_read_header(cprFilename)
        V = mha_read_volume(info);

        %Process
        [dims]=size(V);
        [fx fy fz]=gradient(V);
        x0 = round(dims(1)/2); y0 = round(dims(2)/2);
        sliceOffset=round(dims(3)*.05);

        for steps=sliceOffset:1:dims(3)-sliceOffset       % loop through CPR

            %% roisetup
            t = 0:pi/20:2*pi;
            radio=4;
            R0 = radio/info.PixelDimensions(1); 
            x0 = round(dims(1)/2); y0 = round(dims(2)/2);
            xi = R0*cos(t)+x0;
            yi = R0*sin(t)+y0;
            roimask = poly2mask(xi,yi, dims(1),dims(2));
            pr_r = find(roimask);
            
            Slice=V(:,:,steps);
            [sliceFeat, xcenters] =hist(Slice(pr_r),255);
            
            
            %w = waitforbuttonpress;

            vesselTrainData=[vesselTrainData;sliceFeat];
            if trainingMode==0
                if reference(steps,6)>1.0
                    vesselyTrain=[vesselyTrain; logical(1)];
                else
                    vesselyTrain=[vesselyTrain; logical(0)];
                end
            elseif trainingMode==1
                if reference(steps,6)==2.0
                    vesselyTrain=[vesselyTrain; logical(1)];
                else
                    vesselyTrain=[vesselyTrain; logical(0)];
                end
            else
                vesselyTrain=[vesselyTrain; reference(steps,6)];
            end

        end   %rof process vessel
        
        
        %key=input('key');
    end      %rof process dataset
    %w = waitforbuttonpress;
    trainData=[trainData;vesselTrainData];
    yTrain=[yTrain; vesselyTrain];
end
toc
%% Writing dataset
% outFile=strcat('tmp_data/','training1')
% save(outFile, 'trainData','yTrain');
%matlab2csv('BT_L5_train.csv',trainData,[],2);

if trainingMode==0
    save 'tmp_data/BC_Circle_allvesselsTrain' trainData yTrain
elseif trainingMode==1
   save 'tmp_data/BC_Circle_allvesselsTrain_OnlyC' trainData yTrain
else
    save 'tmp_data/MC_Circle_allvesselsTrain' trainData yTrain
end






