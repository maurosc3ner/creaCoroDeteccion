% function [corteString softM softStd indSoft...
%     calcM calcStd indCalc mixM mixStd indMix]= ...
%     sliceStatus(softM, softStd, indSoft,...
%     calcM,calcStd,indCalc,mixM,mixStd,indMix,...
%     pos,illness, roimean,roistd)
% 
% corteString='Sano';
%     if illness==1.0
%         corteString='Suave';
%         softM=[softM, roimean];
%         softStd=[softStd, roistd];
%         indSoft=[indSoft, pos];
%     elseif illness==2.0
%         corteString='Calcificada';
%         calcM=[calcM, roimean];
%         calcStd=[calcStd, roistd];
%         indCalc=[indCalc, pos];
%     elseif illness==3.0
%         corteString='Mixta';
%         mixM=[mixM, roimean];
%         mixStd=[mixStd, roistd];
%         indMix=[indMix, pos];
%     end
% end