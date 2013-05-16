%***************************
%*  3D Cylinder Model:
%*  I,Gr,Gt,Radon
%   Longitudinal 
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
trainingMode=2;

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
tic
%% Segment lecture
    vesselTrainData=[];
    vesselyTrain=[];
%busqueda de segmentos
MHD= dir(fullfile(inDir,'*.mhd'));
for j =1:numel(MHD),
    
    %% Files lecture

    refFilename=fullfile(inDir,[MHD(j).name(1:end-4) '.txt'])
    
    cprFilename=fullfile(inDir,MHD(j).name);
    reference=load(refFilename);
    info = mha_read_header(cprFilename)
    V = mha_read_volume(info);
    numves=numves+1;
    %Process
    [dims]=size(V);
    [fx fy fz]=gradient(V);
    x0 = round(dims(1)/2); y0 = round(dims(2)/2);
    
    
    %key=input('key input');
    
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
        
        vesselTrainData=[vesselTrainData;feature];
        if trainingMode==0
            if reference(z0,6)>1.0
                vesselyTrain=[vesselyTrain; 1];
            else
                vesselyTrain=[vesselyTrain; 0];
            end
        elseif trainingMode==1
            if reference(z0,6)==2.0
                vesselyTrain=[vesselyTrain; 1];
            else
                vesselyTrain=[vesselyTrain; 0];
            end
            
        elseif trainingMode==2
            if reference(z0,6)>1.0
                vesselyTrain=[vesselyTrain; 1];
            elseif reference(z0,7)>1
                vesselyTrain=[vesselyTrain; 1];
            else
                vesselyTrain=[vesselyTrain; 0];
            end
        elseif trainingMode==3
            if reference(z0,6)>1.0
                vesselyTrain=[vesselyTrain; 1];
            elseif reference(z0,7)>1
                vesselyTrain=[vesselyTrain; 2];
            else
                vesselyTrain=[vesselyTrain; 0];
            end
        elseif trainingMode==4
            if reference(z0,7)>1.0
                vesselyTrain=[vesselyTrain; 1];
            else
                vesselyTrain=[vesselyTrain; 0];
            end
        else
            vesselyTrain=[vesselyTrain; reference(z0,6)];
        end

        
    end   %rof process vessel
    
    trainData=[trainData;vesselTrainData];
    yTrain=[yTrain; vesselyTrain];
end
toc



if trainingMode==0
    save 'tmp_data/tm0_L5_Segments288' trainData yTrain
elseif trainingMode==1
   save 'tmp_data/tm1_L5_Segments288' trainData yTrain
elseif trainingMode==2
    save 'tmp_data/tm2_L5_Segments288' trainData yTrain
elseif trainingMode==3
    save 'tmp_data/tm3_L5_Segments288' trainData yTrain
elseif trainingMode==4
    save 'tmp_data/tm4_L5_Segments288' trainData yTrain
else
    save 'tmp_data/MC_L7_allvesselsTrain' trainData yTrain
end