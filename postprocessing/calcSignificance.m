function [p] = calcSignificance(data,model)

n = length(data);
xi = data(:,1);
yi = data(:,2);

logxi = log10(xi);
logyi = log10(yi);

a = model(1);
b = model(2);

logyiHat = a*logxi+b;

xBar = mean(logxi);

SE = sqrt(sum((logyi-logyiHat).^2)/(n-2))/sqrt(sum((logxi-xBar).^2));

t = a/SE;

p = 1-tcdf(t,n-2)

end