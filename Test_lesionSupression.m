clear all; clc; close all;
addpath src

load 'tmp_data/Yfit'

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

        lesCounter=1;
    else
        lesCounter=lesCounter+1;
    end
end
lesionSize=[lesionSize lesCounter];
indGr=[indGr indInt(i+1)];

indGr
lesionSize
indIzq=indGr-[lesionSize-1]

%%Post-Processing (deleting small lesions)
Yfit2=Yfit;
for i=1:numel(lesionSize)
    
    if lesionSize(i)<3
        Yfit2(indIzq(i):indGr(i))=0;
    end
end


indInt2=find(Yfit2~=0);

plot(Yfit,'g')
grid on;
hold on
plot(Yfit2,'dr')
hold off


izq=[];
izq=[1 izq];
lesiones=[];
for i=1:numel(indGr)
    
    if lesionSize(i)<3
        lesiones=[lesiones; round(mean(indInt(izq(i):indGr(i))))];
        izq=[indGr(i)+1 izq];
    
    end
end
izq