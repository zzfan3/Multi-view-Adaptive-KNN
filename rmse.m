

function r=rmse(data,estimate)

I = ~isnan(data) & ~isnan(estimate); 
data = data(I); estimate = estimate(I);

r=sqrt(sum((data(:)-estimate(:)).^2)/numel(data));