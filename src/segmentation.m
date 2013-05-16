function I=segmentation(R,BW,pixelLab,L,W,P,colorTransform,se1,se2,offset)
Lab = applycform(R, colorTransform);
Lab= lab2uint8(Lab);
spW=4;
%subplot(spW,spW,1); imshow(R),colorbar;title('input');

l=Lab(:,:,1);
a=Lab(:,:,2);
b=Lab(:,:,3);
% Extract the a,b information

LRoi=zeros(480,640);
ARoi=zeros(480,640);
BRoi=zeros(480,640);
if pixelLab(1,1)<10
    %disp('negro');
    thlow=pixelLab(1,1)*(1-offset);
    thhigh=pixelLab(1,1)*(1+offset);
    LRoi=roicolor(l, thlow, thhigh);
else if pixelLab(1,1)>180
        %disp('blanco');
        thlow=pixelLab(1,1)*(1-offset);
        thhigh=255;
        LRoi=roicolor(l, thlow, thhigh);
    else
        thlow=pixelLab(1,2)*(1-offset);
        thhigh=pixelLab(1,2)*(1+offset);
        ARoi = roicolor(a, thlow, thhigh);
        thlow=pixelLab(1,3)*(1-offset);
        thhigh=pixelLab(1,3)*(1+offset);
        BRoi = roicolor(b, thlow, thhigh);
        %disp('color');
    end
end

BW=double(BW);
LRoi=BW.*LRoi;
LRoi(1:480,1:L(2)-P)=0;
LRoi(1:480,L(2)+W(2)+P:640)=0;
LRoi(1:L(1)-P,1:640)=0;
LRoi(L(1)+W(1)+P:480,1:640)=0;

ARoi(1:480,1:L(2)-P)=0;
ARoi(1:480,L(2)+W(2)+P:640)=0;
ARoi(1:L(1)-P,1:640)=0;
ARoi(L(1)+W(1)+P:480,1:640)=0;

BRoi(1:480,1:L(2)-P)=0;
BRoi(1:480,L(2)+W(2)+P:640)=0;
BRoi(1:L(1)-P,1:640)=0;
BRoi(L(1)+W(1)+P:480,1:640)=0;
BRoi=imdilate(BRoi,se2);
% subplot(spW,spW,2); imshow(LRoi), colormap(gray),colorbar;title('LRoi');
% subplot(spW,spW,3); imshow(ARoi), colormap(gray),colorbar;title('ARoi');
% subplot(spW,spW,4); imshow(BRoi), colormap(gray),colorbar;title('BRoi');

I=LRoi| (ARoi&BRoi);

%subplot(spW,spW,5); imshow(I), colormap(gray),colorbar;title('LRoi| (ARoi&BRoi)');
I=imdilate(I,se2);
%subplot(spW,spW,6); imshow(I),colorbar;title('Dilate');
I=imerode(I,se1);
%subplot(spW,spW,6); imshow(I), colormap(gray),colorbar;title('Erode');
I=imfill(I,'holes');
%subplot(spW,spW,7); imshow(I), colormap(gray),colorbar;title('fill');



end