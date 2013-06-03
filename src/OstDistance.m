function distance=OstDistance(ref)

    daccVector=[];
    accu=0.0;
    daccVector=[daccVector; accu];
    for i=2:size(ref,1)
        dis= sqrt(sum( ([ref(i-1,1:3)]-[ref(i,1:3)] ).^2 ));
        accu=accu+dis;
        daccVector=[daccVector;  accu];
    end
    %Normalize
    distance=daccVector/max(daccVector);
end