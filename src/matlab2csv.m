function st = matlab2csv(filename,input, target, mode, featureNames)
% Convert matlab data to a compatible .csv file for use by weka.
%
% input          - An n-by-d matrix with n, d-featured examples.
%
% target         - A n-by-1 vector that contains the class type of
%                  each row.
%
% mode           - If 0, the target values are written as float.
%                  If 1, the target values are written as A,B,C,..
%                - If 2, the target values are written as integers.
%                - If 3, no target values (only test datasets)
%                  In the last case, use integers 1 to 26.
%
% featureNames   - A cell array of d strings, naming each feature/attribute
%                - If not specified , the names get default values.
%
% Written by Alex E.

[n,d]=size(input);
file = fopen(filename,'w');

if(nargin < 5)
    for i=1:d
        fprintf(file,['a' int2str(i) ',']);
    end
else
    for i=1:d
        fprintf(file,'%s,',featureNames{i});
    end
end
fprintf(file,'class\n');

for i=1:n
    for j=1:d
        fprintf(file,'%f,',input(i,j));
    end
    if (mode==1)
        asc = char(64+target(i):64+target(i));
        fprintf(file,'%c\n',asc);
    elseif (mode==0)
        fprintf(file,'%f\n',target(i));
    elseif (mode==2)
        fprintf(file,'%d\n',target(i));
    end
end

st = fclose(file);

end

