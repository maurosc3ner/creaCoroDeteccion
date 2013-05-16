function diff = anisodiff(im, niter, kappa, lambda, option)

if ndims(im)==3
  error('Anisodiff only operates on 2D grey-scale images');
end

im = double(im);
[rows,cols] = size(im);
diff = im;

%{
var = 2;
x = (-4:4);
g = exp(-x.*x/(2*var)); g  = g/sum(g);

blurred = conv2(im,g,'same');
im_b = conv2(blurred,g','same'); 

%}

for i = 1:niter
 % fprintf('\rIteration %d',i);

  % Construct diffl which is the same as diff but
  % has an extra padding of zeros around it.
  diffl = zeros(rows+2, cols+2);
  diffl(2:rows+1, 2:cols+1) = diff;

  % North, South, East and West differences
  deltaN = diffl(1:rows,2:cols+1)   - diff;
  
  deltaS = diffl(3:rows+2,2:cols+1) - diff;
  
  deltaE = diffl(2:rows+1,3:cols+2) - diff;
  
  deltaW = diffl(2:rows+1,1:cols)   - diff;
  %deltaN = diff;deltaW;
  

  % Conduction

  if option == 1
    cN = exp(-(deltaN/kappa).^2);
    
    cS = exp(-(deltaS/kappa).^2);
    cE = exp(-(deltaE/kappa).^2);
    cW = exp(-(deltaW/kappa).^2);
    
  elseif option == 2
    cN = 1./(1 + (deltaN/kappa).^2);
    cS = 1./(1 + (deltaS/kappa).^2);
    cE = 1./(1 + (deltaE/kappa).^2);
    cW = 1./(1 + (deltaW/kappa).^2);
  end

  % APPLYING FOUR-POINT-TEMPLETE FOR numerical solution of DIFFUSION P.D.E.
  
  diff = diff + lambda*(cN.*deltaN + cS.*deltaS + cE.*deltaE + cW.*deltaW);
%figure();
%  Uncomment the following to see a progression of images
 %subplot(ceil(sqrt(niter)),ceil(sqrt(niter)), i)
% figure();
 %pause(1)
%imagesc(diff), colormap(gray), axis image

end
%fprintf('\n');