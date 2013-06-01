%***************************
%*  3D Cylinder Model:
%*  I,Gr,Gt,Radon
%   Longitudinal 
%   Selective segment training
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

%% Lecture mode (decision criters)
% 0-> .mhd slowly
% 1-> .mat compressed (faster)
lectureMode=0;

%% Cylinder mask creation
min_r=1;
max_r=3.5;
rSampl=3;
radiusStep=linspace(min_r,max_r,rSampl);
L=5;

t = pi/4:pi/4:2*pi;
theta = 0:180;
trainData=[];
yTrain=[];
numves=0;

%% Selective Segment lecture

if trainingMode==0
load 'tmp_data/TM0_90_oob_list'
elseif trainingMode==2
load 'tmp_data/TM2_oob_list'
elseif trainingMode==4
load 'tmp_data/TM4_oob_list'
end
MHD=[hsList(trnf);ssList(trpf)];

tic
for j =1:numel(MHD),
    
    %% Files lecture
    if lectureMode==0
        refFilename=fullfile(inDir,[MHD(j).name(1:end-4) '.txt'])
        cprFilename=fullfile(inDir,MHD(j).name);
        reference=load(refFilename);
        info = mha_read_header(cprFilename)
        V = mha_read_volume(info);
    elseif lectureMode==1
       load ([inDir MHD(j).name(1:end-4) '.mat'])
    end
    numves=numves+1
    %Process
    [dims]=size(V);
    [fx fy fz]=gradient(V);
    x0 = round(dims(1)/2); y0 = round(dims(2)/2);
    
    %key=input('key input');
    segmentTrainData=[];
    segmentyTrain=[];
    for z0=round(L/2):1:dims(3)-round(L/2)+1     % loop through CPR
        
        % 3D Sampling pattern
        longitudinalIntensityFeature=zeros(27,L);
        feature=[];
        lIndex=1;
        for k=z0-round(L/2)+1:z0+round(L/2)-1;

            cylFeature=[];
            Slice=V(:,:,k);
            patternSlice=uint8(zeros(dims(1),dims(2)));
            if visualDebug
                imagesc(patternSlice),axis square;
            end
            
            for radio=radiusStep                % loop through radius' scales
                R0 = radio/info.PixelDimensions(1);
                yi = R0*cos(t);
                xi = R0*sin(t);
                
                % 2D Sampling pattern
                pr_r=sub2ind([dims(1) dims(2)],round(x0+xi),round(y0+yi));
                patternSlice(pr_r)=1;
                
                %% Intensity feature
                intensitySlice=Slice(pr_r)';

                %% Gradient feature
                gradxSlice=fx(:,:,k);
                gradySlice=fy(:,:,k);
                grad_xi=[gradxSlice(pr_r)' gradySlice(pr_r)'];
                
                %radial direction points
                ui=[[(x0+xi)-x0]' [(y0+yi)-y0]'];
                ui=ui./norm(ui);
                ui=ui.*2;
                % radial gradients array
                radialGrad=diag(grad_xi(:,:)*ui(:,:)');

                %tangent direction points
                ti=[-ui(:,2) ui(:,1)];
                % tangent gradients array
                tangentGrad=diag(grad_xi(:,:)*ti(:,:)');

                %radon Feature
                [RT,xp] = radon(Slice(round(dims(1)/2-R0):round(dims(1)/2+R0),...
                    round(dims(2)/2-R0):round(dims(2)/2+R0)),theta);
                
                %quin
                radialFeature=[[min(intensitySlice) max(intensitySlice) mean(intensitySlice)]...
                    [min(tangentGrad) max(tangentGrad) mean(tangentGrad)]...
                    [min(radialGrad) max(radialGrad) mean(radialGrad)]...
                    min(RT(:)) max(RT(:)) mean(RT(:))];
                
                %Queuing 12 features to the 36 features
                cylFeature=[cylFeature radialFeature];
                
                if visualDebug
                    hold on
                    imagesc(patternSlice),axis square;
                    quiver(x0+xi',y0+yi',ui(:,1),ui(:,2))
                    quiver(x0+xi',y0+yi',ti(:,1),ti(:,2),'Color','r');
                    colormap gray
                    hold off
                    pause(0.5)
                end
                %key=input('input');
                
            end    %rof radius step
            %radio=radiusStep(end);
            %R0 = radio/info.PixelDimensions(1);
            
            %[RT,xp] = radon(Slice(round(dims(1)/2-R0):round(dims(1)/2+R0),...
            %       round(dims(2)/2-R0):round(dims(2)/2+R0)),theta);
            
            
            longitudinalFeature(lIndex,:)=[cylFeature];
            feature=[feature cylFeature];
            lIndex=lIndex+1;
            %w = waitforbuttonpress;
        end    %rof cylinder height=L
        % 1:9 first circle, 10:18 second one...
        A=[min(longitudinalFeature,[],1) max(longitudinalFeature,[],1) mean(longitudinalFeature,1) ];
        
        feature=[feature A];
        %w = waitforbuttonpress;
        
        segmentTrainData=[segmentTrainData;feature];
        if trainingMode==0
            if reference(z0,6)>1.0
                segmentyTrain=[segmentyTrain; 1];
            else
                segmentyTrain=[segmentyTrain; 0];
            end
        elseif trainingMode==1
            if reference(z0,6)==2.0
                segmentyTrain=[segmentyTrain; 1];
            else
                segmentyTrain=[segmentyTrain; 0];
            end
            
        elseif trainingMode==2
            if reference(z0,6)>1.0
                segmentyTrain=[segmentyTrain; 1];
            elseif reference(z0,7)>1
                segmentyTrain=[segmentyTrain; 1];
            else
                segmentyTrain=[segmentyTrain; 0];
            end
        elseif trainingMode==3
            if reference(z0,6)>1.0
                segmentyTrain=[segmentyTrain; 1];
            elseif reference(z0,7)>1
                segmentyTrain=[segmentyTrain; 2];
            else
                segmentyTrain=[segmentyTrain; 0];
            end
        elseif trainingMode==4
            if reference(z0,7)>1.0
                segmentyTrain=[segmentyTrain; 1];
            else
                segmentyTrain=[segmentyTrain; 0];
            end
        else
            segmentyTrain=[segmentyTrain; reference(z0,6)];
        end

        
    end   %rof process vessel
    
    trainData=[trainData;segmentTrainData];
    yTrain=[yTrain; segmentyTrain];
end
toc

%% Selective Testing Segment lecture

MHD=[hsList(tstnf);ssList(tstpf)];
testData=[];
yTest=[];
tic
for j =1:numel(MHD),
    
    %% Files lecture
    if lectureMode==0
        refFilename=fullfile(inDir,[MHD(j).name(1:end-4) '.txt'])
        cprFilename=fullfile(inDir,MHD(j).name);
        reference=load(refFilename);
        info = mha_read_header(cprFilename)
        V = mha_read_volume(info);
    elseif lectureMode==1
       load ([inDir MHD(j).name(1:end-4) '.mat'])
    end
    numves=numves+1
    %Process
    [dims]=size(V);
    [fx fy fz]=gradient(V);
    x0 = round(dims(1)/2); y0 = round(dims(2)/2);
    
    %key=input('key input');
    segmentTrainData=[];
    segmentyTrain=[];
    for z0=round(L/2):1:dims(3)-round(L/2)+1     % loop through CPR
        
        % 3D Sampling pattern
        longitudinalIntensityFeature=zeros(27,L);
        feature=[];
        lIndex=1;
        for k=z0-round(L/2)+1:z0+round(L/2)-1;

            cylFeature=[];
            Slice=V(:,:,k);
            patternSlice=uint8(zeros(dims(1),dims(2)));
            if visualDebug
                imagesc(patternSlice),axis square;
            end
            
            for radio=radiusStep                % loop through radius' scales
                R0 = radio/info.PixelDimensions(1);
                yi = R0*cos(t);
                xi = R0*sin(t);
                
                % 2D Sampling pattern
                pr_r=sub2ind([dims(1) dims(2)],round(x0+xi),round(y0+yi));
                patternSlice(pr_r)=1;
                
                %% Intensity feature
                intensitySlice=Slice(pr_r)';

                %% Gradient feature
                gradxSlice=fx(:,:,k);
                gradySlice=fy(:,:,k);
                grad_xi=[gradxSlice(pr_r)' gradySlice(pr_r)'];
                
                %radial direction points
                ui=[[(x0+xi)-x0]' [(y0+yi)-y0]'];
                ui=ui./norm(ui);
                ui=ui.*2;
                % radial gradients array
                radialGrad=diag(grad_xi(:,:)*ui(:,:)');

                %tangent direction points
                ti=[-ui(:,2) ui(:,1)];
                % tangent gradients array
                tangentGrad=diag(grad_xi(:,:)*ti(:,:)');

                %radon Feature
                [RT,xp] = radon(Slice(round(dims(1)/2-R0):round(dims(1)/2+R0),...
                    round(dims(2)/2-R0):round(dims(2)/2+R0)),theta);
                
                %quin
                radialFeature=[[min(intensitySlice) max(intensitySlice) mean(intensitySlice)]...
                    [min(tangentGrad) max(tangentGrad) mean(tangentGrad)]...
                    [min(radialGrad) max(radialGrad) mean(radialGrad)]...
                    min(RT(:)) max(RT(:)) mean(RT(:))];
                
                %Queuing 12 features to the 36 features
                cylFeature=[cylFeature radialFeature];
                
                if visualDebug
                    hold on
                    imagesc(patternSlice),axis square;
                    quiver(x0+xi',y0+yi',ui(:,1),ui(:,2))
                    quiver(x0+xi',y0+yi',ti(:,1),ti(:,2),'Color','r');
                    colormap gray
                    hold off
                    pause(0.5)
                end
                %key=input('input');
                
            end    %rof radius step
            %radio=radiusStep(end);
            %R0 = radio/info.PixelDimensions(1);
            
            %[RT,xp] = radon(Slice(round(dims(1)/2-R0):round(dims(1)/2+R0),...
            %       round(dims(2)/2-R0):round(dims(2)/2+R0)),theta);
            
            
            longitudinalFeature(lIndex,:)=[cylFeature];
            feature=[feature cylFeature];
            lIndex=lIndex+1;
            %w = waitforbuttonpress;
        end    %rof cylinder height=L
        % 1:9 first circle, 10:18 second one...
        A=[min(longitudinalFeature,[],1) max(longitudinalFeature,[],1) mean(longitudinalFeature,1) ];
        
        feature=[feature A];
        %w = waitforbuttonpress;
        
        segmentTrainData=[segmentTrainData;feature];
        if trainingMode==0
            if reference(z0,6)>1.0
                segmentyTrain=[segmentyTrain; 1];
            else
                segmentyTrain=[segmentyTrain; 0];
            end
        elseif trainingMode==1
            if reference(z0,6)==2.0
                segmentyTrain=[segmentyTrain; 1];
            else
                segmentyTrain=[segmentyTrain; 0];
            end
            
        elseif trainingMode==2
            if reference(z0,6)>1.0
                segmentyTrain=[segmentyTrain; 1];
            elseif reference(z0,7)>1
                segmentyTrain=[segmentyTrain; 1];
            else
                segmentyTrain=[segmentyTrain; 0];
            end
        elseif trainingMode==3
            if reference(z0,6)>1.0
                segmentyTrain=[segmentyTrain; 1];
            elseif reference(z0,7)>1
                segmentyTrain=[segmentyTrain; 2];
            else
                segmentyTrain=[segmentyTrain; 0];
            end
        elseif trainingMode==4
            if reference(z0,7)>1.0
                segmentyTrain=[segmentyTrain; 1];
            else
                segmentyTrain=[segmentyTrain; 0];
            end
        else
            segmentyTrain=[segmentyTrain; reference(z0,6)];
        end

        
    end   %rof process vessel
    
    testData=[testData;segmentTrainData];
    yTest=[yTest; segmentyTrain];
end
toc


if trainingMode==0
    save 'tmp_data/tm0_90_oob_L5_Segments288' trainData yTrain testData yTest
elseif trainingMode==1
   save 'tmp_data/tm1_oob_L5_Segments288' trainData yTrain testData yTest
elseif trainingMode==2
    save 'tmp_data/tm2_oob_L5_Segments288' trainData yTrain testData yTest
elseif trainingMode==3
    save 'tmp_data/tm3_oob_L5_Segments288' trainData yTrain testData yTest
elseif trainingMode==4
    save 'tmp_data/tm4_oob_L5_Segments288' trainData yTrain testData yTest
else
    save 'tmp_data/MC_L7_allvesselsTrain' trainData yTrain
end