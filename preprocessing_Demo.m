%***************************
%*  3D Cylinder Model for Training and Testing feature matrix
%*  Features: I,Gr,Gt,Radon, Longitudinal, Distance.
%   
%   Validation strategy: Out of bag
%*  2013,4,11
%*  References :
%*  Mittal10
%*  
%***************************
%   Author: Esteban Correa, maurosc3ner@gmail.com
%******************************
mkdir tmp_data;
clear all; clc; close all;
addpath src
addpath ../ReadData3D_version1k/mha

%% Flags
visualDebug=false;
includeDistance=false;
%% Training mode (decision criteria)
% 0-> binary 0=healthy,soft|1=calc,mix
% 1-> binary 0=healthy,soft,mix|1=calc
trainingMode=0;

%% Cylinder mask creation
% Parameters to define the cylinder: minimun radius, maximum radius,
% number of radiuses to evaluate, cylinder height.
min_r=1;
max_r=3.5;
rSampl=3;
radiusStep=linspace(min_r,max_r,rSampl);
L=5;
t = pi/4:pi/4:2*pi;
theta = 0:180;

%% dataset lecture
%dataset indexing
inDir = 'Training_SelVes/';
DT= dir(fullfile(inDir,'dt*'));

tic
trainData=[];
yTrain=[];
for j =1:numel(DT),
    
    %% Vessel lecture
    inDT=strcat(inDir,DT(j).name);
    % Every vessel within the dataset folder is indexed.
    MHD= dir(fullfile(inDT,'*.mhd'));
    
    vesselTrainData=[];
    vesselyTrain=[];
    for vessel_i=1:numel(MHD),
        % A vessel is loaded (reference file points, image.mhd)
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
        % Volume gradient is computed one time.
        [fx fy fz]=gradient(V);
        % Cylinder coordinates are defined, center, limits, among
        % others.
        x0 = round(dims(1)/2); y0 = round(dims(2)/2);
        sliceOffset=round(dims(3)*.05);
        % Surfing through the CPR
        for steps=sliceOffset:1:dims(3)-L-sliceOffset       
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
                 % loop through radius' scales
                for radio=radiusStep               
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

                end    %rof radius step
                
                longitudinalFeature(lIndex,:)=[cylFeature];
                feature=[feature cylFeature];
                lIndex=lIndex+1;
                
            end    %rof cylinder height=L
            % 1:9 first circle, 10:18 second one...
            A=[min(longitudinalFeature,[],1) max(longitudinalFeature,[],1) mean(longitudinalFeature,1) ];

            feature=[feature A];
            
            % adding distance feature
            if includeDistance
                feature=[feature reference(z0,8)];
            end
            vesselTrainData=[vesselTrainData;feature];
            
            %% Adding labels
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
            else
                vesselyTrain=[vesselyTrain; reference(z0,6)];
            end
          
        end   %rof process vessel
        
    end      %rof process dataset
    trainData=[trainData;vesselTrainData];
    yTrain=[yTrain; vesselyTrain];
end
toc


%% Out-of-bag Testing
inDir = 'Testing_SelVes/';

% dataset indexing
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
                end    %rof radius step

                longitudinalFeature(lIndex,:)=[cylFeature];
                feature=[feature cylFeature];
                lIndex=lIndex+1;
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
            %% others    
            else
                vesselyTest=[vesselyTest;...
                [reference(z0,1) reference(z0,2) reference(z0,3) reference(z0,6)]];
            end
            
        end   %rof process vessel
    end      %rof process dataset
    testData=[testData;vesselTestData];
    yTest=[yTest; vesselyTest];
end
toc


if trainingMode==0
    save 'tmp_data/featureMatrix' trainData yTrain testData yTest
end