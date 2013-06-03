%***************************
%*  3D Cylinder Model for Training and Testing
%*  Features: I,Gr,Gt,Radon, Longitudinal, Distance.
%   
%   Validation strategy: Out of bag
%*  2013,4,11
%*  References :
%*  Mittal10
%*  
%***************************
%
%******************************


mkdir tmp_data;

clear all; clc; close all;
addpath src
addpath ../ReadData3D_version1k/mha
inDir = 'Training_SelVes/';

%busqueda de carpetas
DT= dir(fullfile(inDir,'dt*'));

%% Flags
visualDebug=false;
includeDistance=false;
%% Training mode (decision criters)
% 0-> binary 0=healthy,soft|1=calc,mix
% 1-> binary 0=healthy,soft,mix|1=calc
% 2-> binary 0 healthy|1=calc,mix,narr>50
% 3-> multi 0=healthy|1=calc,mix|2=narr>50
% 4-> binary 0=healthy|1=grade narrowing>50
trainingMode=0;

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
        % calling Ostium distance feature
        dist=OstDistance(reference);
        % attaching to the reference file
        reference=[reference, dist];
        
        info = mha_read_header(cprFilename)
        V = mha_read_volume(info);
        numves=numves+1;
        %Process
        [dims]=size(V);
        [fx fy fz]=gradient(V);
        x0 = round(dims(1)/2); y0 = round(dims(2)/2);
        sliceOffset=round(dims(3)*.05);

        for steps=sliceOffset:1:dims(3)-L-sliceOffset       % loop through CPR

            z0=round(L/2)+steps;
            % 3D Sampling pattern
            longitudinalIntensityFeature=zeros(27,L);
            feature=[];
            lIndex=1;
            for k=z0-round(L/2)+1:z0+round(L/2)-1;
                %VM(:,:,k)=patternSlice(:,:);
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

                    %Mask initialization
                    %             VM=uint8(zeros(dims(1),dims(2),20));
                    %             patternSlice=uint8(zeros(dims(1),dims(2)));

                    % 2D Sampling pattern
                    pr_r=sub2ind([dims(1) dims(2)],round(x0+xi),round(y0+yi));
                    patternSlice(pr_r)=1;

                    %% Intensity feature
                    intensitySlice=Slice(pr_r)';
                    %[min(intensitySlice) max(intensitySlice) mean(intensitySlice)];

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
                    % radial feature
                    %[min(radialGrad) max(radialGrad) mean(radialGrad)];

                    %tangent direction points
                    ti=[-ui(:,2) ui(:,1)];
                    % tangent gradients array
                    tangentGrad=diag(grad_xi(:,:)*ti(:,:)');
                    % tangent gradient feature
                    %[min(tangentGrad) max(tangentGrad) mean(tangentGrad)];

                    %radon Feature

                    [RT,xp] = radon(Slice(round(dims(1)/2-R0):round(dims(1)/2+R0),...
                        round(dims(2)/2-R0):round(dims(2)/2+R0)),theta);

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
            
            % adding distance feature
            if includeDistance
                feature=[feature reference(z0,8)];
            end
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
           % yTrain=[yTrain; reference(z0,6)];

            %imagesc(V(:,:,z0)),axis square;colormap gray, title(strcat('cylinder pos:',num2str(z0)));
            %pause(0.1)
            
        end   %rof process vessel
        
        
        %key=input('key');
    end      %rof process dataset
    %w = waitforbuttonpress;
    trainData=[trainData;vesselTrainData];
    yTrain=[yTrain; vesselyTrain];
end
toc

numves

%% Selective Testing Segment lecture
inDir = 'Testing_SelVes/';

%busqueda de carpetas
DT= dir(fullfile(inDir,'dt*'));
testData=[];
yTest=[];

tic
%% dataset lecture
for j =1:numel(DT),
    
    %% Files lecture
    inDT=strcat(inDir,DT(j).name);
    %busqueda de vasos
    MHD= dir(fullfile(inDT,'*.mhd'));
    
    vesselTestData=[];
    vesselyTest=[];
    
    for vessel_i=1:numel(MHD),
       
        refFilename=fullfile(inDT,[MHD(vessel_i).name(1:end-4) '.txt'])
        
        cprFilename=fullfile(inDT,MHD(vessel_i).name);
        reference=load(refFilename);
        % calling Ostium distance feature
        dist=OstDistance(reference);
        % attaching to the reference file
        reference=[reference, dist];
        
        info = mha_read_header(cprFilename)
        V = mha_read_volume(info);

        %Process
        [dims]=size(V);
        [fx fy fz]=gradient(V);
        x0 = round(dims(1)/2); y0 = round(dims(2)/2);
        sliceOffset=round(dims(3)*.05);

        for steps=sliceOffset:1:dims(3)-L-sliceOffset       % loop through CPR

            z0=round(L/2)+steps;
            % 3D Sampling pattern
            longitudinalIntensityFeature=zeros(27,L);
            feature=[];
            lIndex=1;
            for k=z0-round(L/2)+1:z0+round(L/2)-1;
                %VM(:,:,k)=patternSlice(:,:);
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
                    %[min(intensitySlice) max(intensitySlice) mean(intensitySlice)];

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
                    % radial feature
                    %[min(radialGrad) max(radialGrad) mean(radialGrad)];

                    %tangent direction points
                    ti=[-ui(:,2) ui(:,1)];
                    % tangent gradients array
                    tangentGrad=diag(grad_xi(:,:)*ti(:,:)');
                    % tangent gradient feature
                    %[min(tangentGrad) max(tangentGrad) mean(tangentGrad)];

                    %radon Feature
                    [RT,xp] = radon(Slice(round(dims(1)/2-R0):round(dims(1)/2+R0),...
                        round(dims(2)/2-R0):round(dims(2)/2+R0)),theta);

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

                longitudinalFeature(lIndex,:)=[cylFeature];
                feature=[feature cylFeature];
                lIndex=lIndex+1;
                %w = waitforbuttonpress;
            end    %rof cylinder height=L
            % 1:9 first circle, 10:18 second one...
            A=[min(longitudinalFeature,[],1) max(longitudinalFeature,[],1) mean(longitudinalFeature,1) ];

            feature=[feature A];
            
            % adding distance feature
            if includeDistance
                feature=[feature reference(z0,8)];
            end
            vesselTestData=[vesselTestData;feature];
            %%
            if trainingMode==0
                if reference(z0,6)>1.0
                    vesselyTest=[vesselyTest;...
                    [reference(z0,1) reference(z0,2) reference(z0,3) 1.0]];
                else
                    vesselyTest=[vesselyTest;...
                    [reference(z0,1) reference(z0,2) reference(z0,3) 0.0]];
                end 
            %% binary 0=healthy,soft,mix|1=calc
            elseif trainingMode==1
                disp('modo 1')
                if reference(z0,6)==2.0
                    vesselyTest=[vesselyTest;...
                    [reference(z0,1) reference(z0,2) reference(z0,3) 1.0]];
                else
                    vesselyTest=[vesselyTest;...
                    [reference(z0,1) reference(z0,2) reference(z0,3) 0.0]];
                end
            %%    
            elseif trainingMode==2
                if reference(z0,6)>1.0
                    vesselyTest=[vesselyTest; ...
                        [reference(z0,1) reference(z0,2) reference(z0,3) 1.0]];
                elseif reference(z0,7)>1
                    vesselyTest=[vesselyTest;...
                        [reference(z0,1) reference(z0,2) reference(z0,3) 1.0]];
                else
                    vesselyTest=[vesselyTest; ...
                        [reference(z0,1) reference(z0,2) reference(z0,3) 0.0]];
                end    
            %% multi 0=healthy|1=calc,mix|2=narr>50
            elseif trainingMode==3
                if reference(z0,6)>1.0
                    vesselyTest=[vesselyTest; ...
                        [reference(z0,1) reference(z0,2) reference(z0,3) 1.0]];
                elseif reference(z0,7)>1
                    vesselyTest=[vesselyTest;...
                        [reference(z0,1) reference(z0,2) reference(z0,3) 2.0]];
                else
                    vesselyTest=[vesselyTest; ...
                        [reference(z0,1) reference(z0,2) reference(z0,3) 0.0]];
                end
            %% binary 0=healthy|1=grade narrowing>50
            elseif trainingMode==4
                disp('modo 4')
                if reference(z0,7)>1.0
                    vesselyTest=[vesselyTest;...
                    [reference(z0,1) reference(z0,2) reference(z0,3) 1.0]];
                else
                    vesselyTest=[vesselyTest;...
                    [reference(z0,1) reference(z0,2) reference(z0,3) 0.0]];
                end    
            %% Others 
            else
                vesselyTest=[vesselyTest;...
                [reference(z0,1) reference(z0,2) reference(z0,3) reference(z0,6)]];
            end
            
        end   %rof process vessel
        
        
        %key=input('key');
    end      %rof process dataset
    %w = waitforbuttonpress;
    testData=[testData;vesselTestData];
    yTest=[yTest; vesselyTest];
end
toc


if trainingMode==0
    save 'tmp_data/tm0_oob_L5_SelVes288' trainData yTrain testData yTest
end