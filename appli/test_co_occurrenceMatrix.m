I = imread('circuit.tif');

%offset all directions 
% pos     Angle
% [0 1;     0
% -1 1;     45
% -1 0;     90
% -1 -1;    135
% 0 -1;     180
% 1 -1;     225
% 1 0;      270
% 1 1]      315
Iroi=I(100:150,100:170);
[GLCM1 SI] = graycomatrix(I,'Offset',[-1 0;0 1],'NumLevels',8);
stats = graycoprops(GLCM1,{'contrast','homogeneity'})

entropy(I)
entropy(Iroi)
entropy(SI(:))

GLCM2 = graycomatrix(Iroi,'Offset',[-1 0;0 1]);
stats = graycoprops(GLCM2,{'contrast','homogeneity'})