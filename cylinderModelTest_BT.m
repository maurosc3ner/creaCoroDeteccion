mkdir tmp_data;

clear all; clc; close all;
addpath src
addpath ../ReadData3D_version1k/mha
inDir = 'test/';
D= dir(fullfile(inDir,'dt*'));

%% Flags
visualDebug=false;
saveResults=true;

%% Cylinder mask creation
min_r=1;
max_r=3.5;
rSampl=3;
radiusStep=linspace(min_r,max_r,rSampl);
L=5;

t = pi/4:pi/4:2*pi;
theta = 0:180;
testData=[];
yTest=[];
%% dataset lecture
for j=1:numel(D),
    
    %% Files lecture
    datasetDirectory=fullfile(inDir,D(j).name)
    refFilename=strcat(datasetDirectory,'/','_creandpuj_reference_rca_path.txt');
    cprFilename=strcat(datasetDirectory,'/','_creandpuj_3Dimage_rca.mhd');
    reference=load(refFilename);
    info = mha_read_header(cprFilename);
    V = mha_read_volume(info);
    [dims]=size(V);
    [fx fy fz]=gradient(V);
    x0 = round(dims(1)/2); y0 = round(dims(2)/2);
    
    
    datasetData=[];
    datasetTarget=[];
    for steps=0:1:dims(3)-L       % loop through CPR
        
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
            
            longitudinalFeature(lIndex,:)=[cylFeature];
            feature=[feature cylFeature];
            lIndex=lIndex+1;
            %w = waitforbuttonpress;
        end    %rof cylinder height=L
        % 1:9 first circle, 10:18 second one...
        A=[min(longitudinalFeature,[],1) max(longitudinalFeature,[],1) mean(longitudinalFeature,1) ];
        
        feature=[feature A];
        %w = waitforbuttonpress;
        datasetData=[datasetData;feature];
        
        if reference(z0,6)>1.0
            datasetTarget=[datasetTarget;...
            [reference(z0,1) reference(z0,2) reference(z0,3) 1.0]];
        
        else
            datasetTarget=[datasetTarget;...
            [reference(z0,1) reference(z0,2) reference(z0,3) 0.0]];
        
        end
        
        %datasetTarget=[datasetTarget;...
        %    [reference(z0,1) reference(z0,2) reference(z0,3) reference(z0,6)]];
        
        imagesc(V(:,:,z0)),axis square;colormap gray, title(strcat('cylinder pos:',num2str(z0)));
        %pause(0.1)
        
    end   %rof process dataset
    
    if saveResults
            tdfl=strcat('tmp_data/',D(j).name,'TestData.txt');
            dlmwrite(tdfl,datasetData,'Delimiter',' ',...
                'Newline','unix','precision',6);
            ydfl=strcat('tmp_data/',D(j).name,'YData.txt');
            dlmwrite(ydfl,datasetTarget,'Delimiter',' ',...
                'Newline','unix','precision',6);
            
    end
    
    testData=[testData; datasetData];
    yTest=[yTest; datasetTarget];
    %w = waitforbuttonpress;

end

load 'tmp_data/rusboost2000_BT'
%% Test
figure;
tic
plot(loss(rusTree,datasetData,datasetTarget(:,4),'mode','cumulative'));
toc
grid on;
xlabel('Number of trees');
ylabel('Test classification error');

% check confusion matrix
tic
Yfit = predict(rusTree,datasetData);
toc
tab2=ones(2,1);
tab = tabulate(datasetTarget(:,4));
for i=1:size(tab,1)
    tab2(i,1)=tab(i,2);
end
cm=confusionmat(datasetTarget(:,4),Yfit);
cm2=bsxfun(@rdivide,cm,tab2(:))*100;

% Measures
TPT=diag(cm2);
TP=TPT(1)
TN=TPT(2)
FP=cm(1,2)
FN=cm(2,1)

TPR=TP/(TP+FN)
FPR=FP/(FP+TN)
ACC=(TP+TN)/(TP+FN+FP+TN)
SPC=1-FPR
PPV=TP/(TP+FP)

yffl=strcat('tmp_data/',D(j).name,'BTYFit.txt');
dlmwrite(yffl,[datasetTarget(:,1:3) Yfit],'Delimiter',' ',...
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
yffl=strcat('tmp_data/',D(j).name,'BTL5_LesionPositions.txt');
dlmwrite(yffl,[datasetTarget(lesiones(:)',1:3) Yfit(lesiones(:)')],'Delimiter',' ',...
                'Newline','unix','precision',6);

