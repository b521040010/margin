function ret=pdfVGx(z,sigma,v,t)  
m=0;
theta=0;

x=z - m.*t - (t./v).*log(1-theta.*v-sigma.^2.*v./2);
     a=(2.* exp((theta.*x)./sigma.^2)) ./  (v.^(t./v)*sqrt(2*pi).*sigma.*gamma(t./v));
     b=( (x.^2) ./ ((2*sigma.^2)./v+theta.^2) ) .^(t./(2.*v)-1/4);
     c=besselk( t./v-1/2 , (1./(sigma.^2)) *sqrt( (x.^2).* ( (2*sigma^2)/v+theta^2 ) ) );
     ret=a.*b.*c;
end