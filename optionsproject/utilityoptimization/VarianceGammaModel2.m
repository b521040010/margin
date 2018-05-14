classdef VarianceGammaModel2 < Model1D
    
    properties
        S0;
        m;
        theta;
        sigma;
        v;
        T;
    end
    
    methods
        function model = VarianceGammaModel2()            
            model.S0 = 1000.0;
            model.m=0;
            model.theta = 0;
            model.sigma = 0.1257;
            model.v = 0.0027;
            model.T = 1/12;
        end
       
        function [S, times] = simulatePricePaths( model, nPaths, nSteps )
            error('not implemented yet for nSteps>1');
        end
        
        function res = pdf(m,x)
        % Compute the pdf of the model
            if size(x,2)>1
                error('x must be a vector')
            end
            ret=zeros(length(x),1);
            for i=1:length(x)
%                 fun = @(g) (1./(m.sigma.*sqrt(2.*pi.*g))).*exp(-(log(x(i))-(m.theta.*g+log(m.S0)+(m.m+(1/m.v)*log(1-m.theta*m.v-m.sigma^2*m.v/2))*m.T)).^2./(2.*m.sigma.^2.*g)).*((g.^(m.T./m.v-1).*exp(-g./m.v))./(m.v.^(m.T./m.v).*gamma(m.T./m.v)));
%                 ret(i) = integral(fun,0,inf);
                   ret(i) = pdfVG2(x(i),m.m,m.theta,m.sigma,m.v,m.T,m.S0);
            end    
            res=(1./x).*ret;
        end         
        
        function res = logPdf(m,x)
        % Compute the log pdf of the model
            res=log(m.pdf(x));            
        end          
        
        function res = logPdfEstimate2(bsm,x,nPaths)
            [S] = simulatePricePaths( bsm, nPaths, 1 );
            prices = S(:,end);
            [numbersInBins,points]=hist(prices,x);
            res=log(numbersInBins/nPaths/(x(2)-x(1)));
        end
        
        function res = logPdfImpSamp(bsm,x)
        % Compute the log pdf of the model
%             [mu,sigma] = logNormalParameters(bsm);  
%             nu=5;
%             
%             s = log(x);
%             res2 = log(gamma((nu+1)/2)/(gamma(nu/2)*sqrt(pi*nu)*sigma)) ...
%                 + (-(nu+1)/2).*log( 1 + 1/nu*((s-mu)/sigma).^2) ...
%                 - s;
%             
%             res1 = -((-mu + log(x)).^2/(2*sigma^2)) - log((sqrt(2*pi)*sigma*x));    
% 
%             res=res1-res2;
            error('not yet implemented') 
        end    
        
        function r = mean(bsm)
        % Mean of the model
%             [m,s] = logNormalParameters(bsm);   
%             r = exp(m + s^2*0.5);
            error('not yet implemented') 
        end
        
        function r = sd(bsm)
        % SD of the model
%             [m,s] = logNormalParameters(bsm);
%             r = sqrt((exp(s^2)-1)*exp(2*m + s^2 ));
            error('not yet implemented') 
        end
        
        
        function model = fit( bsm, S0, T, returns, weights)
            % Fit the model to historic return data            
%             fittedDist = fitdist( returns + 1, 'LogNormal', 'frequency', weights);
%             m = fittedDist.mu;
%             s = fittedDist.sigma;
%             bsm.S0 = S0;
%             bsm.T = T;
%             model = bsm.fitLogNormalReturn(m,s);
                error('not yet implemented') 
        end
        
        function wayPoints = getWayPoints(model)
        % Returns some standard way points for accurate numeric integration
             mean = model.S0();
             wayPoints = 1000:1000 + mean;
             wayPoints = wayPoints( wayPoints>0 );         
        end
        
        function [ prices, deltas ] = price(o, r, isPut, strike, spot)
            %BLACKSCHOLESPRICE Computes the black scholes price
            %  for an array of European options
%             assert( islogical( isPut ));
%             vol = o.sigma;
%             if (nargin<5)
%                 spot = o.S0;
%             end
% 
%             d1 = 1./(vol .* sqrt(o.T)) .* (log( spot./strike ) + (r+0.5*vol.^2).*o.T);
%             d2 = 1./(vol .* sqrt(o.T)) .* (log( spot./strike ) + (r-0.5*vol.^2).*o.T);
% 
%             callPrice = spot .* normcdf(d1) - exp( -r .* o.T ) .* strike .* normcdf(d2);
%             putPrice = -spot .* normcdf(-d1) + exp( -r .* o.T ) .* strike .* normcdf(-d2);
%             prices = callPrice.*(1-isPut) + isPut.*putPrice;
%             deltas = normcdf(d1) - 1 * isPut;
            error('not yet implemented') 
        end        
        
        function sigma = impliedVolatility( m, r, strike, isPut, optionPrice )
        % Compute the implied volatility of an option
        %optionPrice=optionPrice/100;
%             function priceDiff = f( sigma )
%                 m.sigma = sigma;
%                 priceDiff = m.price( r, isPut, strike ) - optionPrice;
%             end
%             options = optimset('fsolve');
%             options = optimset(options, 'Display', 'off');
%             [sigma,~, ret] = fsolve( @f, 0.1, options );
%             if ret<=0
%                 fprintf('Cannot find implied volatility\n');
%                 fprintf( 'r=%f, T=%f, S0=%f, strike=%f, isPut=%f, optionPrice%f\n', r, m.T, m.S0, strike, isPut, optionPrice );
%                 error('Failed');
%             end
            error('not yet implemented') 
        end
        
    end
end