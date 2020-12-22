close all
clc

for i=1:100
    GAMMA(i) = rand(1);
end 

x = linspace(0,50,length(GAMMA));

MAX = max(GAMMA);
MIN = min(GAMMA);

DELTA = (MAX - MIN)+1;

figure(1)
hold on
for i=1:length(GAMMA)
    k = GAMMA(i)/DELTA;
    plot(x(i),GAMMA(i),'o','Color',[k*1.3,k*0.8,0],'LineWidth',4);
end