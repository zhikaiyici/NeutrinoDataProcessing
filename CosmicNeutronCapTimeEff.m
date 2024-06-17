clear;
close all;
timeWindow = 292;
coinciRate = timeWindow .* 55.44 .* 0.0134 .* 0.69 .* 16e4 ./ 26.61 ./ 1e6;
% coinciRate = 54.56 ./ 2;
coinciEvent = coinciRate .* 3600 .* 277 .* 0.23 .* 0.2426 .* 0.9275;

taun = 57.91;
deff1 = @(tn) exp((tn - 8) ./ taun) - exp((tn - 300) ./ taun);
deff2 = @(tn) exp(-8 ./ taun) - exp((tn - 300) ./ taun);

figure;
tn = -timeWindow:0.1:0;
plot(tn, deff1(tn));

n1 = sum(deff1(tn));
effCal1 = n1 ./ length(tn);

hold on;
tn = 0:0.1:timeWindow;
plot(tn, deff2(tn));
hold off;

n2 = sum(deff2(tn));
effCal2 = n2 ./ length(tn);

effCal = (n1 + n2) ./ 2 ./ length(tn);

syms tn;
deff1 = exp((tn - 8) ./ taun) - exp((tn - 300) ./ taun);
deff2 = exp(-8 ./ taun) - exp((tn - 300) ./ taun);

eff = double(int(deff1, tn, -timeWindow, 0) + int(deff2, tn, 0, timeWindow)) ./ timeWindow ./ 2;
eff1 = double(int(deff1, tn, -inf, 0)) ./ timeWindow;
eff2 = double(int(deff2, tn, 0, timeWindow)) ./ timeWindow;

cosicNeutron = coinciEvent .* eff1;

%%
y = rand(1e6, 1);
x = - taun .* log(1 - y); 
figure;
histogram(x);
xi = -timeWindow + timeWindow .* 2 .* rand(size(y));
figure;
histogram(xi);
x_ = x + xi;
figure;
histogram(x_);
effMC = sum(x_ > 8 & x_ < 300) ./ length(y);

xi1 = -timeWindow + timeWindow .* rand(size(y));
figure;
histogram(xi1);
x1 = x + xi1;
figure;
histogram(x1);
effMC1 = sum(x1 > 8 & x1 < 300) ./ length(y);

xi2 = timeWindow .* rand(size(y));
figure;
histogram(xi2);
x2 = x + xi2;
figure;
histogram(x2);
effMC2 = sum(x2 > 8 & x2 < 300) ./ length(y);

