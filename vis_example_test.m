

diff=anisodiff(Slice,10,20,0.25,2);
se = strel('disk',12);
J1 = imtophat(diff,se);
J2 = imtophat(Slice,se);

subplot(2,2,1)
imagesc(Slice), colormap(gray),axis image
subplot(2,2,2)
imagesc(diff), colormap(gray),axis image
subplot(2,2,3)
imagesc(J1), colormap(gray),axis image
subplot(2,2,4)
imagesc(J2), colormap(gray),axis image


r=2; n=20;
[x,y,z]=sphere(n);
figure(2)
mesh(r*x,r*y,r*z)
title('Sphere')
axis equal

XC=40;
YC=40;
ZC=40;
RX=5;
RY=5;
RZ=5;
[X,Y,Z]=ELLIPSOID(XC,YC,ZC,RX,RY,RZ);
mesh(X,Y,Z,'FaceColor','red','EdgeColor','red'),title('(x-XC)^2/RX^2+(y-YC)^2/RY^2+(z-ZC)^2/RZ^2)=1')
hidden off;alpha(.1);
hold on
[x3,y3,z3]=ELLIPSOID(42,42,42,0.2,0.2,0.2);
mesh(x3,y3,z3,'FaceColor','red','EdgeColor','none')
[x2,y2,z2]=ELLIPSOID(36,36,36,0.2,0.2,0.2);
mesh(x2,y2,z2,'FaceColor','blue','EdgeColor','none');
hold off


pInside=[42 42 42];
pOutside=[42 42 42];
ellipsoidequation=(pOutside(1)-XC)^2/RX^2+...
    (pOutside(2)-YC)^2/RY^2+(pOutside(3)-ZC)^2/RZ^2
if ellipsoidequation>1
    fprintf('Punto por fuera');
elseif ellipsoidequation<1
    fprintf('Punto por dentro');
end

[x2,y2,z2]=ELLIPSOID(37,37,37,1,1,1);
mesh(x2,y2,z2)
hold on, plot3(37,37,37,'linewidth',3)
   