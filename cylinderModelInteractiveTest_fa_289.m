%***************************
%*  3D Model Test Binary class:
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
mkdir results;

clear all; clc; close all;
addpath src
addpath ../ReadData3D_version1k/mha
inDir = 'Testing_SelVes/';

%busqueda de carpetas
DT= dir(fullfile(inDir,'dt*'));

% GUI to select dataset
prompt={};
prompt = {'Select the dataset to evaluate:'}
%prompt{1}={[prompt{1} sprintf('\n') 'dataset1' sprintf('\n') 'dataset2']}
for i=1:numel(DT),

    prompt{1}=[prompt{1} sprintf('\n') strcat(num2str(i),')',DT(i).name)];
end

dlg_title = 'Datasets available';
num_lines = 1;
def = {'1'};
answer = inputdlg(prompt{1},dlg_title,num_lines,def);

%% Flags
visualDebug=false;
globalVessel=true;

%% Classifier
% 0-> AdaboostM1
% 1-> RUSBoost
% 2-> Random Forest
classifierMode=0;

%% Cylinder mask creation
min_r=1;
max_r=3.5;
rSampl=3;
radiusStep=linspace(min_r,max_r,rSampl);
L=5;

t = pi/4:pi/4:2*pi;
theta = 0:180;
testData=[];
yTarget=[];


tic
%% dataset lecture
for j=str2num(answer{1}),
    
    %% Files lecture
    inDT=strcat(inDir,DT(j).name);
    %busqueda de vasos
    MHD= dir(fullfile(inDT,'*.mhd'));
    
    % GUI to select vessel
    prompt={};
    prompt = {'Select the vessel to evaluate:'}
    for i=1:numel(MHD),

        prompt{1}=[prompt{1} sprintf('\n') strcat(num2str(i),')',MHD(i).name)];
    end

    dlg_title = 'Vessels available';
    num_lines = 1;
    def = {'1'};
    answer2 = inputdlg(prompt{1},dlg_title,num_lines,def);
    
    vesselTestData=[];
    vesselyTest=[];
    
    for vessel_i=str2num(answer2{1}),
       
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
            %w = waitforbuttonpress;
            feature=[feature reference(z0,4)];
            vesselTestData=[vesselTestData;feature];
            
            vesselyTest=[vesselyTest;...
                [reference(z0,1) reference(z0,2) reference(z0,3)]];
        end   %rof process vessel
        
        
        %key=input('key');
    end      %rof process dataset
    %w = waitforbuttonpress;
    testData=[testData;vesselTestData];
    yTarget=[yTarget; vesselyTest];
end
toc

modelFile = 'rbVSrf_250_tm0_oob_SelVes_289'

if classifierMode==0
   ML = 'adb';
else if classifierMode==1
   ML = 'rusTree';
    end
end

vars={ML};
load(strcat('models/',modelFile,'.mat'),vars{:});

%% Prediction
tic
if classifierMode==0
    Yfit = predict(adb,testData);
else if classifierMode==1
    Yfit = predict(rusTree,testData);
    end
end
toc


%offset insertion to 3D points

yTarget=[reference(1:round(L/2)+sliceOffset-1,1)...
    reference(1:round(L/2)+sliceOffset-1,2)...
    reference(1:round(L/2)+sliceOffset-1,3);
    yTarget(:,1:3);...
    [reference(z0+1:steps+L+sliceOffset-1,1)...
    reference(z0+1:steps+L+sliceOffset-1,2)...
    reference(z0+1:steps+L+sliceOffset-1,3)]]; %reste -1 porque es par ARREGLAR!

Yfit=[zeros(size(1:round(L/2)+sliceOffset-1,2),1); Yfit ;zeros(size(z0+1:steps+L+sliceOffset-1,2),1)];


%% finding illness points to interpolate
indInt=find(Yfit~=0);
indGr=[];
lesCounter=1;
lesionSize=[];
for i=1:size(indInt,1)-1
    if(indInt(i)+1~=indInt(i+1))
        lesionSize=[lesionSize lesCounter];
        disp('diferente grupo')
        indGr=[indGr indInt(i)];
        indInt(i);
        indInt(i+1);
        lesCounter=1;
    else
        lesCounter=lesCounter+1;
    end
end
lesionSize=[lesionSize lesCounter]
indGr=[indGr indInt(i+1)]

%%Post-Processing (deleting small lesions)
lesiones=[];
for i=1:numel(indGr)
    if lesionSize(i)>6
        lesiones=[lesiones; round(mean( (indGr(i)-lesionSize(i)+1:indGr(i)) ))];
    else
        Yfit(indGr(i)-lesionSize(i)+1:indGr(i))=0;
    end
end

yffl=strcat('results/',DT(j).name,MHD(vessel_i).name(1:end-4),ML,'.txt');
dlmwrite(yffl,[yTarget (Yfit)*40],'Delimiter',' ',...
                'Newline','unix','precision',6);
                 
            
% Saving medium points of lesions
yffl=strcat('results/',DT(j).name,MHD(vessel_i).name(1:end-4),ML,'_LPos.txt');
dlmwrite(yffl,[yTarget(lesiones(:)',1:3) Yfit(lesiones(:)')],'Delimiter',' ',...
                'Newline','unix','precision',6);

