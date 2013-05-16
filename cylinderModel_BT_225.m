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

min_r=1;
max_r=3.5;
rSampl=3;
radiusStep=linspace(min_r,max_r,rSampl);
L=5;

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

                    radialFeature=[[min(intensitySlice) max(intensitySlice) mean(intensitySlice)]...
                        [min(tangentGrad) max(tangentGrad) mean(tangentGrad)]...
                        [min(radialGrad) max(radialGrad) mean(radialGrad)]];...
                        %min(RT(:)) max(RT(:)) mean(RT(:))];

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
               
                %radon feature
                R0 = radiusStep(end)/info.PixelDimensions(1);
                [RT,xp] = radon(Slice(round(dims(1)/2-R0):round(dims(1)/2+R0),...
                        round(dims(2)/2-R0):round(dims(2)/2+R0)),theta);
                
                longitudinalFeature(lIndex,:)=[cylFeature min(RT(:)) max(RT(:)) mean(RT(:))];
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
                    vesselyTrain=[vesselyTrain; logical(1)];
                else
                    vesselyTrain=[vesselyTrain; logical(0)];
                end
            elseif trainingMode==1
                if reference(z0,6)==2.0
                    vesselyTrain=[vesselyTrain; logical(1)];
                else
                    vesselyTrain=[vesselyTrain; logical(0)];
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
%% Writing dataset
% outFile=strcat('tmp_data/','training1')
% save(outFile, 'trainData','yTrain');
%matlab2csv('BT_L5_train.csv',trainData,[],2);

if trainingMode==0
    save 'tmp_data/BC_L5_Train225' trainData yTrain
elseif trainingMode==1
   save 'tmp_data/BC_L5_Train225_OnlyC' trainData yTrain
else
    save 'tmp_data/MC_L5_allvesselsTrain' trainData yTrain
end






