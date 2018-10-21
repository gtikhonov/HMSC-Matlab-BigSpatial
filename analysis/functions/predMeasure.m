function [ predMeasure ] = predMeasure( Y, YV )
%PREDMEASURE Summary of this function goes here
%   Detailed explanation goes here

accuracy = mean(mean(abs(Y - YV)));
discrimination = nanmean(sum(Y.*YV)./sum(YV) - sum(Y.*(1-YV))./sum((1-YV)));
sharpness = mean(mean(sqrt(Y.*(1-Y))));

[~,indSort] = sort(Y);
YSort = Y;
YVSort = YV;
for i=1:size(Y,2)
    YSort(:,i) = Y(indSort(:,i),i);
    YVSort(:,i) = YV(indSort(:,i),i);
end
calibration = mean(mean(abs(cumsum(YSort)-cumsum(YVSort))));
deviance = mean(mean(-log(Y.*YV+(1-Y).*(1-YV)),1));

predMeasure = [accuracy,discrimination,calibration,deviance,sharpness];
end

