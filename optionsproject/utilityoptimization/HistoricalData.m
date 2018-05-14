classdef HistoricalData
    properties
        data
        dates
        dateNum
        prices
        
    end
    methods
        function [histD] = HistoricalData(fileData)
            histD.data=readtable(fileData);
            dateChar=datetime(histD.data.SPXIndex,'InputFormat','dd/MM/yyyy');
            histD.dates=dateChar;
            histD.dateNum=datenum(dateChar);
            histD.prices=str2double(histD.data.Var2);
        end
        
        function [selectedPrices,selectedDate] = selectTheIntervals(histD,startingDate,maturity)
            numberOfObservations=252*1;
            numberOfObservations=numberOfObservations+1;
            startingDate=datetime(startingDate,'InputFormat','dd/MM/yyyy');
            maturity=datetime(maturity,'InputFormat','dd/MM/yyyy');
            startingPoint=find(histD.dateNum==datenum(startingDate));
            endingPoint=find(histD.dateNum==datenum(maturity));
            %numberOfDays=endingPoint-startingPoint
            numberOfDays=1;
            selectedDatesNum=histD.dateNum(startingPoint:-numberOfDays:1);
            selectedDatesNum=selectedDatesNum(1:1:numberOfObservations);
            selectedDate=datestr(selectedDatesNum);
            selectedPrices=histD.prices(startingPoint:-numberOfDays:1);
            selectedPrices=selectedPrices(1:1:numberOfObservations);
        end
        
        function [mu,sigma] = calibrateNormal(histD,selectedPrices)
            logReturns=log(selectedPrices(1:end-1))-log(selectedPrices(2:end));
            pd=fitdist(logReturns,'Normal');
            mu=pd.mu
            sigma=pd.sigma
        end
        function [mu,sigma,nu] = calibrateStudentT(histD,selectedPrices)
            logReturns=log(selectedPrices(1:end-1))-log(selectedPrices(2:end));
            pd=fitdist(logReturns,'tLocationScale');
            mu=pd.mu;
            sigma=pd.sigma;
            nu=pd.nu;
        end
        function [sigma,v] = calibrateVarianceGamma(histD,selectedPrices)
            r=log(selectedPrices(1:end-1))-log(selectedPrices(2:end));
            tt=1/252;
            sigma1=sqrt(mean(r.^2)/tt);
            v1=(mean(r.^4)-3*sigma1^4*tt^2)/(3*sigma1^4*tt);
            %a=mle(r,'pdf',@(x,sigma,v,t)pdfVG(x,sigma,v,1/252),'start',[sigma1,v1]);
            [a b]=mle(r,'pdf',@(x,sigma,v,t)pdfVGx(x,sigma,v,1/252) ,'start',[sigma1,0.01]);
            
            
            sigma=a(1);
            v=a(2);
            assert(sigma>0);
            assert(v>0);
        end        
        
    end
end