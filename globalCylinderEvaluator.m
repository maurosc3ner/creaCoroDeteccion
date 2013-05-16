%***************************
%*  Global patient evaluator:
%
%*  2013,4,21
%*  References :
%*  
%*  
%***************************
%
%
% 
%******************************



clear all; clc; close all;
load tmp_data/gvEval
addpath src


% check confusion matrix
tab2=ones(2,1);
tab = tabulate(yTargetTotal(:,4));
for i=1:size(tab,1)
    tab2(i,1)=tab(i,2);
end
cm=confusionmat(yTargetTotal(:,4),yFitTotal);
cm2=bsxfun(@rdivide,cm,tab2(:))*100

% Measures
TP=cm(1,1)
TN=cm(2,2)
FP=cm(1,2)
FN=cm(2,1)

TPR=TP/(TP+FN)
FPR=FP/(FP+TN)
ACC=(TP+TN)/(TP+FN+FP+TN)
SPC=1-FPR
PPV=TP/(TP+FP)
NPV=TN/(FN+TN)

%% finding illness points to interpolate
indInt=find(yFitTotal~=0);
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

% Saving medium points of all patient vessels
yffl=strcat('tmp_data/','dt09_','BTL5_all_LesionPositions.txt');
dlmwrite(yffl,[yTargetTotal(lesiones(:)',1:3) yFitTotal(lesiones(:)')],'Delimiter',' ',...
                'Newline','unix','precision',6);

