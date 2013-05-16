%***************************
%*  3D Model Test Binary class:
%*  225 Features:
%*  I,Gr,Gt,Radon, Longitudinal 
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
inDir = 'test/';

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
                    
                    radialFeature=[[min(intensitySlice) max(intensitySlice) mean(intensitySlice)]...
                        [min(tangentGrad) max(tangentGrad) mean(tangentGrad)]...
                        [min(radialGrad) max(radialGrad) mean(radialGrad)]];%...
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

            vesselTestData=[vesselTestData;feature];

            if trainingMode==0
                if reference(z0,6)>1.0
                    vesselyTest=[vesselyTest;...
                    [reference(z0,1) reference(z0,2) reference(z0,3) 1.0]];
                else
                    vesselyTest=[vesselyTest;...
                    [reference(z0,1) reference(z0,2) reference(z0,3) 0.0]];
                end
            elseif trainingMode==1
                disp('modo 1')
                if reference(z0,6)==2.0
                    vesselyTest=[vesselyTest;...
                    [reference(z0,1) reference(z0,2) reference(z0,3) 1.0]];
                else
                    vesselyTest=[vesselyTest;...
                    [reference(z0,1) reference(z0,2) reference(z0,3) 0.0]];
                end
            else
                vesselyTest=[vesselyTest;...
                [reference(z0,1) reference(z0,2) reference(z0,3) reference(z0,6)]];
            end
            
        end   %rof process vessel
        
        
        %key=input('key');
    end      %rof process dataset
    %w = waitforbuttonpress;
    testData=[testData;vesselTestData];
    yTarget=[yTarget; vesselyTest];
end
toc

%yTarget(:,4)=yTarget(:,4)+1;
%modelFile = 'rb500_CMBC_MP_WC51';
%modelFile = 'rusboost500_BT_allvesselsTrain'
modelFile = 'rb500_CMBC_MP80'
vars={'rusTree'};
load(strcat('models/',modelFile,'.mat'),vars{:});

%% Test

% check confusion matrix
tic
Yfit = predict(rusTree,testData);
toc
tab2=ones(2,1);
tab = tabulate(yTarget(:,4));
for i=1:size(tab,1)
    tab2(i,1)=tab(i,2);
end
cm=confusionmat(yTarget(:,4),Yfit)
cm2=bsxfun(@rdivide,cm,tab2(:))*100

%saving dataset results
if globalVessel
yTargetTotal=[];
yFitTotal=[];
    %load tmp_data/gvEval

    yTargetTotal=[yTargetTotal; yTarget];
    yFitTotal=[yFitTotal; Yfit];
    save tmp_data/gvEval yTargetTotal yFitTotal
end

% Measures
TP=cm(1,1)
TN=cm(2,2)
FP=cm(1,2)
FN=cm(2,1)

SEN=TP/(TP+FN)
FPR=FP/(FP+TN);
SPC=TN/(FP+TN)
ACC=(TP+TN)/(TP+FN+FP+TN)
PPV=TP/(TP+FP)
NPV=TN/(FN+TN)

%offset insertion
yTarget=[reference(1:round(L/2)+sliceOffset-1,1)...
    reference(1:round(L/2)+sliceOffset-1,2)...
    reference(1:round(L/2)+sliceOffset-1,3);
    yTarget(:,1:3);...
    [reference(z0+1:steps+L+sliceOffset,1)...
    reference(z0+1:steps+L+sliceOffset,2)...
    reference(z0+1:steps+L+sliceOffset,3)]]; 
Yfit=[reference(1:round(L/2)+sliceOffset-1,6); Yfit ;reference(z0+1:steps+L+sliceOffset,6)];

yffl=strcat('results/',DT(j).name,MHD(vessel_i).name(1:end-4),modelFile,'.txt');
dlmwrite(yffl,[yTarget Yfit],'Delimiter',' ',...
                'Newline','unix','precision',6);

%% finding illness points to interpolate
indInt=find(Yfit~=0);
indGr=[];
for i=1:size(indInt,1)-1
    if(indInt(i)+1~=indInt(i+1))
        disp('diferente grupo')
        indGr=[indGr i];
        indInt(i);
        indInt(i+1);
    end
    
end
indGr=[indGr i+1];

izq=1;
lesiones=[];
for i=1:numel(indGr)
    lesiones=[lesiones; round(mean(indInt(izq:indGr(i))))];
    izq=indGr(i)+1;
end

% Saving medium points of lesions
yffl=strcat('results/',DT(j).name,MHD(vessel_i).name(1:end-4),modelFile,'_LPos.txt');
dlmwrite(yffl,[yTarget(lesiones(:)',1:3) Yfit(lesiones(:)')],'Delimiter',' ',...
                'Newline','unix','precision',6);

