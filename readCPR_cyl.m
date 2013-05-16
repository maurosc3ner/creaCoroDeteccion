%***************************
%*  3D Feature Extractor
%*  Esteban Correa
%* 
%*  2013,3,12
%*  References :
%*  
%*  
%***************************
%
% Features (mean, std, gradMean, gradStd) are extracted with a 3D region of
% interest.
% ROI is defined by the cylinder equation:
%
% (x-XC)^2 + (y-YC)^2 = r^2, z ? [0 L] 
%******************************

clear all; clc; close all;
patient='dt08'
datasetDirectory=strcat('training/',patient)
addpath src

%% file lecture
refFilename=strcat(datasetDirectory,'/','_creandpuj_reference_rca_path.txt');

cprFilename=strcat(datasetDirectory,'/','_creandpuj_3Dimage_rca.mhd');

reference=load(refFilename);

info = mha_read_header(cprFilename)
V = mha_read_volume(info);
[fx fy fz]=gradient(V);
%% Flags
visualDebug=false;

%% Cylinder mask creation
[dims]=size(V);
min_r=1;
max_r=3.5;
rSampl=3;
radiusStep=linspace(min_r,max_r,rSampl);

t = pi/4:pi/4:2*pi; 
L=5;
x0 = round(dims(1)/2); y0 = round(dims(2)/2);
theta = 0:180;
trainData=[];
for steps=0:1:dims(3)-L        % loop through CPR
    
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
            
            [RT,xp] = radon(Slice(round(dims(1)/2)-R0:round(dims(1)/2)+R0,...
                        round(dims(2)/2)-R0:round(dims(2)/2)+R0),theta);
            
            
            
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
    
    trainData=[trainData;feature];
    
    imagesc(V(:,:,z0)),axis square;colormap gray, title(strcat('cylinder pos:',num2str(z0)));
    pause(0.1)
    %VM=circshift(VM,[0 0 steps]);
end

mkdir tmp_data;
save tmp_data/dt08 trainData

