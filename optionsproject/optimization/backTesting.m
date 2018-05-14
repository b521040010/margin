function backTesting
feb=Dynamic('D20170117T150000','D20170217T150000',1,100000)
feb=feb.run;
mar=Dynamic('D20170217T150000','D20170317T150000',1,feb.histPort.D20170216T150000.payoff(2343.01))
mar=mar.run;
apr=Dynamic('D20170321T150000','D20170421T150000',1,mar.histPort.D20170316T150000.payoff(2383.71))
apr=apr.run;
may=Dynamic('D20170421T150000','D20170519T150000',1,apr.histPort.D20170420T150000.payoff(2354.74))
may=may.run;
june=Dynamic('D20170519T150000','D20170616T150000',1,may.histPort.D20170518T150000.payoff(2371.37))
june=june.run;
july=Dynamic('D20170621T150000','D20170721T150000',1,june.histPort.D20170614T150000.payoff(2431.24))
july=july.run;
aug=Dynamic('D20170724T150000','D20170818T150000',1,july.histPort.D20170720T150000.payoff(2467.40))
aug=aug.run;
sep=Dynamic('D20170818T150000','D20170915T150000',1,aug.histPort.D20170817T150000.payoff(2427.64))
sep=sep.run;
oct=Dynamic('D20170920T150000','D20171020T150000',1,sep.histPort.D20170914T150000.payoff(2495.67))
oct=oct.run;
nov=Dynamic('D20171020T150000','D20171117T150000',1,oct.histPort.D20171019T150000.payoff(2567.56))
nov=nov.run;
jan2018=Dynamic('D20171219T150000','D20180119T150000',1,nov.histPort.D20171116T150000.payoff(2582.94))
jan2018=jan2018.run;
feb2018=Dynamic('D20180122T150000','D20180216T150000',1,jan2018.histPort.D20180118T150000.payoff(2802.60))
feb2018=feb2018.run;
mar2018=Dynamic('D20180221T150000','D20180316T150000',1,feb2018.histPort.D20180215T150000.payoff(2727.14))
mar2018=mar2018.run;

[febMark febSpot febUtility ]=feb.getMarkToMarket;
[marMark marSpot marUtility ]=mar.getMarkToMarket;
[aprMark aprSpot aprUtility ]=apr.getMarkToMarket;
[mayMark maySpot mayUtility ]=may.getMarkToMarket;
[juneMark juneSpot juneUtility ]=june.getMarkToMarket;
[julyMark julySpot julyUtility ]=july.getMarkToMarket;
[augMark augSpot augUtility ]=aug.getMarkToMarket;
[sepMark sepSpot sepUtility ]=sep.getMarkToMarket;
[octMark octSpot octUtility ]=oct.getMarkToMarket;
[novMark novSpot novUtility ]=nov.getMarkToMarket;

[jan2018Mark jan2018Spot jan2018Utility ]=jan2018.getMarkToMarket;
[feb2018Mark feb2018Spot feb2018Utility ]=feb2018.getMarkToMarket;
[mar2018Mark mar2018Spot mar2018Utility ]=mar2018.getMarkToMarket;


% mark=[febMark marMark aprMark mayMark juneMark julyMark augMark sepMark octMark novMark nov.histPort.D20171116T150000.payoff(2582.94)];
% spot=[febSpot marSpot aprSpot maySpot juneSpot julySpot augSpot sepSpot octSpot novSpot 2582.94];
% partitions=[length(febMark) length(marMark) length(aprMark) length(mayMark) length(juneMark) length(julyMark) length(augMark) length(sepMark) length(octMark) length(novMark)]
% partitions=cumsum(partitions)

mark=[febMark marMark aprMark mayMark juneMark julyMark augMark sepMark octMark novMark jan2018Mark feb2018Mark mar2018Mark mar2018.histPort.D20180315T150000.payoff(2750.57)];
spot=[febSpot marSpot aprSpot maySpot juneSpot julySpot augSpot sepSpot octSpot novSpot jan2018Spot feb2018Spot mar2018Spot 2750.57];
partitions=[length(febMark) length(marMark) length(aprMark) length(mayMark) length(juneMark) length(julyMark) length(augMark) length(sepMark) length(octMark) length(novMark) length(jan2018Mark) length(feb2018Mark) length(mar2018Mark)]
partitions=cumsum(partitions)


ax1=subplot(2,1,1);
plot(mark);
hold on
markLength=min(mark)-0.1*min(mark):1000:max(mark)+0.1*max(mark);
ylength=1:1000:15*10^7;
for i=1:length(partitions)
plot([partitions(i)*ones(1,length(markLength))],markLength);
end
axis([0 inf min(markLength) max(markLength)]);
ax2=subplot(2,1,2);
plot(spot);
hold on
spotLength=min(spot)-0.05*min(spot):10:max(spot)+0.05*max(spot);
ylength=2200:1:2600;
for i=1:length(partitions)
plot([partitions(i)*ones(1,length(spotLength))],spotLength)
end
axis([0 inf min(spotLength) max(spotLength)])

figure
plot(log(mark/mark(1)))
hold on
plot(log(spot/spot(1)))
logMarkLength=min(log(mark/mark(1)))-0.25*min(log(mark/mark(1))):0.01:max(log(mark/mark(1)))+0.25*max(log(mark/mark(1)));
for i=1:length(partitions)
plot([partitions(i)*ones(1,length(logMarkLength))],logMarkLength)
end
axis([0 inf min(logMarkLength) max(logMarkLength)])

save feb feb
save mar mar
save apr apr
save may may
save june june
save july july
save aug aug
save sep sep
save oct oct
save nov nov
save jan2018 jan2018
save feb2018 feb2018
save mar2018 mar2018

load feb feb
load mar mar
load apr apr
load may may
load june june
load july july
load aug aug
load sep sep
load oct oct
load nov nov
load jan2018 jan2018
load feb2018 feb2018
load mar2018 mar2018


nov.plotHistPortfolio('D20171116T150000')
sep.plotHistPortfolio('D20170914T150000')
plot(2000:3000,sep.histPort.D20170914T150000.payoff(2000:3000),'displayName','payoff on D20170914T150000')
hold on
%plot(2000:3000,100000*ones(1,length(2000:3000)))
plot(2000:3000,sep.histPort.D20170913T150000.computeMarkToMarket(0)*ones(1,length(2000:3000)),'displayName','mark-to-market on D20170913T150000')
plot(2500*ones(1,length(0:10000:14*10^5)),0:10000:14*10^5,'displayName','spot on D20170914T150000')
plot(2495.67*ones(1,length(0:10000:14*10^5)),0:10000:14*10^5,'displayName','openning price on D20170915T150000')

end
