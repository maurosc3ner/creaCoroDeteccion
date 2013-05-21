%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Vessel evaluator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: E. Correa, March 18, 2013
% This file will help to find consecutive vessel segments
% just modyfing a index order in reference file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all;

inDir = '../';
D= dir(fullfile(inDir,'dataset*'));

%% Flags
visualDebug=false;
saveResults=true;
vesselList=1:14;
%% Files lecture
%% dataset lecture
for dtn=0:numel(D)-1,
    for vsn=1:numel(vesselList),
        vesselFlag=vesselList(vsn);
        path=strcat(inDir,D(dtn+1).name)
        refFilename='reference_CTA.txt';

        reference=load(strcat(path,'/',refFilename));


        %% con estructuras
        field1 = 'st';  value1 = zeros(1,10);
        field2 = 'en';  value2 = zeros(1,10);
        field3 = 'gr';  value3 = {'a', 'b'};

        s = [];



        %% finding illness points to interpolate
        indInt=find(reference(:,4)~=0);
        iter=1;
        start=1;
        for i=1:size(indInt,1)-1

            if(reference(i,4)~=reference(i+1,4))
                %disp('diferente grupo')
                eND=i;
                s=[s struct(field1,start,field2,eND,field3,reference(i,4))];
                iter=iter+1;
                start=eND+1;
            end   
        end
        eND=i+1;

        s=[s struct(field1,start,field2,eND,field3,reference(i,4))];

        %% vessel evaluator 

        if (vesselFlag==1) %RCA corta
            vessel=[1,2,3];
            segString='RCA corta vessels=[1,2,3]';
        elseif (vesselFlag==2) %RCA y R-PDA
            vessel=[1,2,3,4];
            segString='RCA y R-PDA vessels=[1,2,3,4]';
        elseif (vesselFlag==3) %RCA R-PLB 
            vessel=[1,2,3,17];
            segString='RCA y R-PLB vessels=[1,2,3,17]';
        elseif (vesselFlag==4) %LAD D1
            vessel=[5,6,9];
            segString='LAD D1 vessel=[5,6,9]';
        elseif (vesselFlag==5) %LAD corta
            vessel=[5,6,7];
            segString='LAD corta vessel=[5,6,7]';
        elseif (vesselFlag==6) %LAD D2
            vessel=[5,6,7,10];
            segString='LAD D2 vessel=[5,6,7,10]';
        elseif (vesselFlag==7) %LAD completa
            vessel=[5,6,7,8];
            segString='LAD completa vessels=[5,6,7,8]'    
        elseif (vesselFlag==8) %OM1
            vessel=[5,11,12];
            segString='OM1 vessels=[5,11,12]'
        elseif (vesselFlag==9) %OM2
            vessel=[5,11,13,14];
            segString='OM2 vessels=[5,11,13,14]'
        elseif (vesselFlag==10) %LCX muy corta
            vessel=[5,11];
            segString='LCX corta vessels=[5,11]' 
        elseif (vesselFlag==11) %LCX corta
            vessel=[5,11,13];
            segString='LCX corta vessels=[5,11,13]'  
        elseif (vesselFlag==12) %LCX y L-PDA
            vessel=[5,11,13,15];
            segString='LCX y L-PDA vessels=[5,11,13,15]'
        elseif (vesselFlag==13) %RAMUS
            vessel=[5,16];
            segString='RAMUS vessels=[5,16]'
        elseif (vesselFlag==14) %CASO ESPECIAL RCA
            vessel=[1,2];
            segString='RCA muy corta vessels=[1,2]'        
        end
        distacc=0;
        iter=1;

        while iter<numel(vessel)
            seg1=0;
            seg2=0;
            for j=1:numel(s)
                if(s(j).gr == vessel(iter))
                    seg1=j;
                end
                if(s(j).gr == vessel(iter+1))
                    seg2=j;
                end
            end

            if seg1~=0 && seg2~=0
                [s(seg1).gr eucdist(reference(s(seg1).en,1:3),...
                   reference(s(seg2).st,1:3)) s(seg2).gr]
                distacc=distacc+eucdist(reference(s(seg1).en,1:3),...
                    reference(s(seg2).st,1:3));
            elseif seg1==0
                errmsg=strcat('ERROR: Segmento:',num2str(vessel(iter)),' no existente!!!');
                disp(errmsg);
                distacc=10;
                break,
            elseif seg2==0
                errmsg=strcat('ERROR: Segmento ',num2str(vessel(iter+1)),' no existente!!!');
                disp(errmsg);
                distacc=10;
                break,
            end
            iter=iter+1;
        end

        outputPoints=[];
        vn=[];
        if distacc<10
            for i=1:numel(vessel)
                pr_r=find(reference(:,4)==vessel(i));
                outputPoints=[outputPoints; [reference(pr_r,:)]];
                %key=input('key');
                vn=strcat(vn,num2str(vessel(i)));
            end

            yffl=strcat(path,'/vessel_',vn,'.txt');
            dlmwrite(yffl,outputPoints,'Delimiter',' ',...
                        'Newline','unix','precision',6);
            msg=strcat('Los segmentos de la lista:',segString,' son consecutivos!!!');
            disp(msg);
            disp('y han sido extraidos y guardados correctamente');
        else
            errmsg=strcat('ERROR: Los segmentos de la lista:',segString,' no son consecutivos!!!');
            disp(errmsg);
            disp('Se recomienda extraer otra combinacion de segmentos');
        end
        %key=input('key');

    end
end