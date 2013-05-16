function VM=sphereMask(dims,sphxc,sphyc,sphzc,radx,rady,radz)
% create a 3D ellipsoid shape roi  
% Esteban Correa correa@creatis.insa-lyon.fr

VM=uint8(zeros(dims));
    for i=1:dims(1)
        for j=1:dims(2)
            for k=1:dims(3)
                ellipsoidequation=(i-sphxc)^2/radx^2+...
                    (j-sphyc)^2/rady^2+(k-sphzc)^2/radz^2;
                if ellipsoidequation<1

                    VM(i,j,k)=1;
                end
            end
        end
    end

end