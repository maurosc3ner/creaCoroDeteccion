function features_norm = normalize_features(features, nrm)
if nargin<2, nrm = 1; end

% each row is one feature
if nrm==1,
    featnorm = sum(features,2);
else
    featnorm = (sum(features.^nrm,2)).^(1/nrm);
end
featnorm2 = repmat(featnorm,[1,size(features,2)]);
features_norm = features./(featnorm2+(featnorm2==0));

