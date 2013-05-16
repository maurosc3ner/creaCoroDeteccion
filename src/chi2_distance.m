% this implementation is suboptimal in matlab

function [d] = chi2_distance(hists1,hists2)

if size(hists1,2)~=size(hists2,2),
    error('histograms should have same number of bins!');
end

hists1=repmat(hists1,[1,1,size(hists2,1)]);

d=zeros(size(hists1,1),size(hists2,1));

for h=1:size(hists1,1),
    
    hist1=squeeze(hists1(h,:,:))';
    d(h,:) = 0.5*sum((hist1-hists2).*(hist1-hists2)./((hist1+hists2)+((hist1+hists2)==0)),2);
    
end

