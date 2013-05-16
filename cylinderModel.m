mkdir tmp_data;

clear all; clc; close all;
addpath src
addpath ../ReadData3D_version1k/mha
inDir = 'training/';
D= dir(fullfile(inDir,'dt*'));

%% Flags
visualDebug=false;

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
    
    for steps=10:1:dims(3)-L-20       % loop through CPR
        
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
        
        trainData=[trainData;feature];
        
%         if reference(z0,6)>0.0
%             yTrain=[yTrain; 1];
%         else
%             yTrain=[yTrain; 0];
%         end
        yTrain=[yTrain; reference(z0,6)];
        
        imagesc(V(:,:,z0)),axis square;colormap gray, title(strcat('cylinder pos:',num2str(z0)));
        %pause(0.1)
        
    end   %rof process dataset
    
    

    
    
    %w = waitforbuttonpress;
    

    
end

%% Writing dataset
% outFile=strcat('tmp_data/','training1')
% save(outFile, 'trainData','yTrain');
% matlab2csv('test.csv',trainData,[],2);


% soft=find(yTrain==1);
% size(soft)
% calc=find(yTrain==2);
% size(calc)
% mix=find(yTrain==3);
% size(mix)
% negat=find(yTrain<1);
% size(negat)
% yhist=[numel(negat) numel(soft) numel(calc) numel(mix)];
% %[meanVHist,xout3]=hist(yTrain)
% 
% %% Rejection Sampling
% %bar([0 1 2 3],yhist,'histc'); % Histogram Normalization
% f=yhist/sum(yhist);
% fLength=4;
% figure;
% bar([0 1 2 3],f,'histc'); % Histogram normalized
% 
% figure; h = plot(f,'r','Linewidth',2);
% hold on;
% l = plot([1 fLength],[max(f) max(f)],'k','Linewidth',2);
%  
% legend([h,l],{'f(x)','q(x)'},'Location','Southwest');
% xlim([0 fLength + 1])
% xlabel('x');
% ylabel('p(x)');
% title('Target (f(x)) and Proposal (q(x)) Distributions');
% 
% % OUR PROPOSAL IS THE DISCRETE UNIFORM ON THE INTERVAL [1 fLength]
% % SO OUR CONSTANT IS
% c = 8*max(f/(1/fLength));
%  
% nSamples = 2154;
% i = 1;
% while i < nSamples
%    proposal = unidrnd(4);
%    q = c*(1/4); % ENVELOPE DISTRIBUTION
%    if rand > f(proposal)/q
%       samps(i) = proposal;
%       i = i + 1;
%    end
% end
%  
% % DISPLAY THE SAMPLES AND COMPARE TO THE TARGET DISTRIBUTION
% bins = 1:fLength;
% counts = histc(samps,bins)
% figure
% b = bar(1:fLength,counts./sum(counts),'FaceColor',[.8 .8 .8])
% hold on;
% h = plot(f,'r','Linewidth',2)
% legend([h,b],{'f(x)','samples'});
% xlabel('x'); ylabel('p(x)');
% xlim([0 fLength + 1]);
% 
% 
% %% Rejection Sampling 2
% %bar([0 1 2 3],yhist,'histc'); % Histogram Normalization
% f=yhist/sum(yhist);
% fLength=4;
% figure;
% bar([0 1 2 3],f,'histc'); % Histogram normalized
% 
% figure; h = plot(f,'r','Linewidth',2);
% hold on;
% l = plot([1 fLength],[max(f) max(f)],'k','Linewidth',2);
%  
% legend([h,l],{'f(x)','q(x)'},'Location','Southwest');
% xlim([0 fLength + 1])
% xlabel('x');
% ylabel('p(x)');
% title('Target (f(x)) and Proposal (q(x)) Distributions');
% 
% % OUR PROPOSAL IS THE DISCRETE UNIFORM ON THE INTERVAL [1 fLength]
% % SO OUR CONSTANT IS
% x=1:4;
% q = inline('normpdf(x,3,1)','x');
%  
% % DETERMINE SCALING CONSTANT
% c = max(f(x)./q(x))
% 
% nSamples = 2154;
% i = 1;
% while i < nSamples
%    proposal = unidrnd(4);
%    %q = c*(1/4); % ENVELOPE DISTRIBUTION
%    if rand > f(proposal)/q(proposal)
%       samps(i) = proposal;
%       i = i + 1;
%    end
% end
%  
% % DISPLAY THE SAMPLES AND COMPARE TO THE TARGET DISTRIBUTION
% bins = 1:fLength;
% counts = histc(samps,bins)
% figure
% b = bar(1:fLength,counts./sum(counts),'FaceColor',[.8 .8 .8])
% hold on;
% h = plot(f,'r','Linewidth',2)
% legend([h,b],{'f(x)','samples'});
% xlabel('x'); ylabel('p(x)');
% xlim([0 fLength + 1]);

