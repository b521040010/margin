function ret=pdfVG(x,sigma,v,t)
%the pdf was taken from page9 of Madan1998
%     x=z-m.*t-(t./v).*log(1-theta.*v-sigma.^2.*v/2)
%     %a=(2*exp(theta.*x./sigma.^2))
%     %v.^(t./v).*sqrt(2*pi)*sigma*gamma(t./v)
%     a=(2*exp(theta.*x./sigma.^2))./(v.^(t./v).*sqrt(2*pi)*sigma*gamma(t./v))
%     %./(v.^(t./v).*sqrt(2*pi)*sigma*gamma(t./v)
%     b=((x.^2)/(2*sigma.^2./v+theta.^2)).^(t./(2*v)-1/4)
%     c=bessely(t./v-1/2,(1./sigma.^2)*sqrt(x.^2*(2*sigma.^2./v+theta.^2)))
%     ret=a.*b.*c;

%fun = @(g) (1./(sigma.*sqrt(2.*pi.*g))).*exp(-(x-theta.*g).^2./(2*sigma.^2.*g)).*((g.^(t./v-1).*exp(-g./v))./(v.^(t./v)*gamma(t./v)));
% t=1/365;
% sigma=sqrt(mean(r.^2)/t);
% v=(mean(r.^4)-3*sigma^4*t^2)/(3*sigma^4*t)
% sigma=0.1257;
% v=0.0027;
%Command for moment methods
% fileName = strcat( '../SPXFuturesAndOptions/','SPXdataOpen','.xlsx');
% histD = HistoricalData(fileName);
% s=histD.prices
% s=s(2:end)
% r=log(s(2:end))-log(s(1:end-1))
% t=1/365
% sigma=sqrt(mean(r.^2)/t);
% v=(mean(r.^4)-3*sigma^4*t^2)/(3*sigma^4*t)
% sigma=0.2384
% v=0.0093

%command for maximum likelihood
%mle(r,'pdf',@(x,theta,sigma,v,t)pdfVG(x,theta,sigma,v,1/365),'start',[0,0.1257,0.0027])
% a = [sigma , v]
% 
%     0.2243    0.0035
% 
% 
% b =
% 
%     0.2197    0.0032
%     0.2290    0.0039
theta=0;
for i=1:length(x)
    fun = @(g) (1./(sigma.*sqrt(2.*pi.*g))).*exp(-(x(i)-theta.*g).^2./(2.*sigma.^2.*g)).*((g.^(t./v-1).*exp(-g./v))./(v.^(t./v).*gamma(t./v)));
    ret(i) = integral(fun,0,inf);
end
ret=ret';
end