M = 4;
N = 8;
psi = 21/5;
Aq = sqrt(6/((M^2-1)*psi^2 + N^2-1));
m = 4.8;
ai = 2*(1-1/M);
aq = 2*(1-1/N);
snr = -10:1:60;
gb = 10.^(snr/10);
const = m^m./(gb.^m * gamma(m));
nTerms = length(snr);

a = Aq;
I1exact = zeros(1,nTerms);
for t = 1:nTerms
    func = @(x)(erfc(a*sqrt(x)/sqrt(2)).*x.^(m-1).*exp(-m/gb(t).*x));
    I1exact(t) = 0.5*const(t)*integral(func, 0, Inf);
end

a = psi*Aq;
I1exact1 = zeros(1,nTerms);
for t = 1:nTerms
    func = @(x)(erfc(a*sqrt(x)/sqrt(2)).*x.^(m-1).*exp(-m/gb(t).*x));
    I1exact1(t) = 0.5*const(t)*integral(func, 0, Inf);
end


% using kara
a = Aq;
I1kara = (a^2/2+ m./gb).^(-m+1/2)*gamma(m-1/2) - (a^2/2+m./gb).^(-m).*(sqrt(a^2/2+m./gb)*gamma(m-1/2).*hypergeom(m-1/2,1/2,1.98^2*a^2./(8*(a^2/2+m./gb))) - 1.98*a/sqrt(2)*gamma(m)*hypergeom(m,3/2,1.98^2*a^2./(8*(a^2/2+m./gb))));
I1kara = const/(1.135*sqrt(2*pi)*a).*I1kara;


I1Prony2 = m^m./(2*gb.^m).*(0.416*(m./gb+1.942*a^2/2).^(-m) + 0.294*(m./gb+1.050*a^2/2).^(-m));

I1Prony3 = m^m./(2*gb.^m).*(0.336*(m./gb+a^2*1.752/2).^(-m) + 0.288*(m./gb+1.050*a^2/2).^(-m) + 0.004*(m./gb+1.206*a^2/2).^(-m));

I1Olab1 = m^m./(2*gb.^m).*(0.4803*(m./gb+1.1232*a^2/2).^(-m));

I1Olab2 = m^m./(2*gb.^m).*(0.3017*(m./gb+1.0510*a^2/2).^(-m) + 0.4389*(m./gb+2.102*a^2/2).^(-m));

I1Olab3 = m^m./(2*gb.^m).*(0.3357*(m./gb+a^2*1.0649/2).^(-m) + 0.3361*(m./gb+2.1298*a^2/2).^(-m) + 0.0305*(m./gb+3.1947*a^2/2).^(-m));

I1gh = const*0.6405/sqrt(pi).*(m./gb+a^2/2).^(-1/2-m).*(sqrt(m./gb+a^2/2).*gamma(m).*hypergeom(m,1/2,0.6004^2*a^2./(8*(m./gb+a^2/2)))-0.6004*a/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,0.6004^2*a^2./(8*(m./gb+a^2/2))));

I1gh = I1gh + const*0.2457/sqrt(pi).*(m./gb+a^2/2).^(-1/2-m).*(sqrt(m./gb+a^2/2).*gamma(m).*hypergeom(m,1/2,2.5048^2*a^2./(8*(m./gb+a^2/2)))-2.5048*a/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,2.5048^2*a^2./(8*(m./gb+a^2/2))));

I1Chiani = m^m./(gb.^m).*(1/12*(m./gb+a^2/2).^(-m) + 1/4*(m./gb+2/3*a^2).^(-m));

% using SYA3
I1SYA3 = m^m./(gb.^m) .*(1/12*(m./gb+a^2/2).^(-m) + 1/6*(m./gb+2*a^2).^(-m) + 1/6*(m./gb+2/3*a^2).^(-m));

% using SYA4
I1SYA4 = m^m./(gb.^m) .*(1/16*(m./gb+a^2/2).^(-m) + 1/8*(m./gb+a^2).^(-m) + 1/8*(m./gb+10/3*a^2).^(-m) + + 1/8*(m./gb+10/17*a^2).^(-m));

I1mbfc1 = const*exp(-0.6966).*(m./gb+0.769*a^2/2).^(-1/2-m).*(sqrt(m./gb+0.769*a^2/2).*gamma(m).*hypergeom(m,1/2,1.0797^2*a^2./(8*(m./gb+0.769*a^2/2)))-1.0797*a/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,1.0797^2*a^2./(8*(m./gb+0.769*a^2/2))));

I1mbfc2 = const*exp(-0.9049).*(m./gb+0.9548*a^2/2).^(-1/2-m).*(sqrt(m./gb+0.9548*a^2/2).*gamma(m).*hypergeom(m,1/2,0.6341^2*a^2./(8*(m./gb+0.9548*a^2/2)))-0.6341*a/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,0.6341^2*a^2./(8*(m./gb+0.9548*a^2/2))));

I1ampj = const*0.5.*(m./gb+0.748*a^2/2).^(-1/2-m).*(sqrt(m./gb+0.748*a^2/2).*gamma(m).*hypergeom(m,1/2,1.0988^2*a^2./(8*(m./gb+0.748*a^2/2)))-1.0988*a/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,1.0988^2*a^2./(8*(m./gb+0.748*a^2/2))));



a = psi*Aq;
I1kara1 = (a^2/2+ m./gb).^(-m+1/2)*gamma(m-1/2) - (a^2/2+m./gb).^(-m).*(sqrt(a^2/2+m./gb)*gamma(m-1/2).*hypergeom(m-1/2,1/2,1.98^2*a^2./(8*(a^2/2+m./gb))) - 1.98*a/sqrt(2)*gamma(m)*hypergeom(m,3/2,1.98^2*a^2./(8*(a^2/2+m./gb))));
I1kara1 = const/(1.135*sqrt(2*pi)*a).*I1kara1;

I1Prony21 = m^m./(2*gb.^m).*(0.416*(m./gb+1.942*a^2/2).^(-m) + 0.294*(m./gb+1.050*a^2/2).^(-m));

I1Prony31 = m^m./(2*gb.^m).*(0.336*(m./gb+a^2*1.752/2).^(-m) + 0.288*(m./gb+1.050*a^2/2).^(-m) + 0.004*(m./gb+1.206*a^2/2).^(-m));

I1Olab11 = m^m./(2*gb.^m).*(0.4803*(m./gb+1.1232*a^2/2).^(-m));

I1Olab21 = m^m./(2*gb.^m).*(0.3017*(m./gb+1.0510*a^2/2).^(-m) + 0.4389*(m./gb+2.102*a^2/2).^(-m));

I1Olab31 = m^m./(2*gb.^m).*(0.3357*(m./gb+a^2*1.0649/2).^(-m) + 0.3361*(m./gb+2.1298*a^2/2).^(-m) + 0.0305*(m./gb+3.1947*a^2/2).^(-m));

I1gh1 = const*0.6405/sqrt(pi).*(m./gb+a^2/2).^(-1/2-m).*(sqrt(m./gb+a^2/2).*gamma(m).*hypergeom(m,1/2,0.6004^2*a^2./(8*(m./gb+a^2/2)))-0.6004*a/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,0.6004^2*a^2./(8*(m./gb+a^2/2))));

I1gh1 = I1gh1 + const*0.2457/sqrt(pi).*(m./gb+a^2/2).^(-1/2-m).*(sqrt(m./gb+a^2/2).*gamma(m).*hypergeom(m,1/2,2.5048^2*a^2./(8*(m./gb+a^2/2)))-2.5048*a/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,2.5048^2*a^2./(8*(m./gb+a^2/2))));

I1Chiani1 = m^m./(gb.^m).*(1/12*(m./gb+a^2/2).^(-m) + 1/4*(m./gb+2/3*a^2).^(-m));

% using SYA3
I1SYA31 = m^m./(gb.^m) .*(1/12*(m./gb+a^2/2).^(-m) + 1/6*(m./gb+2*a^2).^(-m) + 1/6*(m./gb+2/3*a^2).^(-m));

% using SYA4
I1SYA41 = m^m./(gb.^m) .*(1/16*(m./gb+a^2/2).^(-m) + 1/8*(m./gb+a^2).^(-m) + 1/8*(m./gb+10/3*a^2).^(-m) + + 1/8*(m./gb+10/17*a^2).^(-m));



I1mbfc11 = const*exp(-0.6966).*(m./gb+0.769*a^2/2).^(-1/2-m).*(sqrt(m./gb+0.769*a^2/2).*gamma(m).*hypergeom(m,1/2,1.0797^2*a^2./(8*(m./gb+0.769*a^2/2)))-1.0797*a/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,1.0797^2*a^2./(8*(m./gb+0.769*a^2/2))));

I1mbfc21 = const*exp(-0.9049).*(m./gb+0.9548*a^2/2).^(-1/2-m).*(sqrt(m./gb+0.9548*a^2/2).*gamma(m).*hypergeom(m,1/2,0.6341^2*a^2./(8*(m./gb+0.9548*a^2/2)))-0.6341*a/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,0.6341^2*a^2./(8*(m./gb+0.9548*a^2/2))));

I1ampj1 = const*0.5.*(m./gb+0.748*a^2/2).^(-1/2-m).*(sqrt(m./gb+0.748*a^2/2).*gamma(m).*hypergeom(m,1/2,1.0988^2*a^2./(8*(m./gb+0.748*a^2/2)))-1.0988*a/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,1.0988^2*a^2./(8*(m./gb+0.748*a^2/2))));



a =Aq;
b = psi*Aq;

I2exact = zeros(1,nTerms);
for t = 1:nTerms
    func = @(x)(erfc(a*sqrt(x)/sqrt(2)).*erfc(b*sqrt(x)/sqrt(2)).*x.^(m-1).*exp(-m/gb(t).*x));
    I2exact(t) = 0.25*const(t)*integral(func, 0, Inf);
end


I2kara1 = (a^2/2+ b^2/2+ m./gb).^(-m+1)*gamma(m-1);
I2kara2 = (a^2/2+ b^2/2+m./gb).^(-m+1/2).*(sqrt(a^2/2+ b^2/2+m./gb)*gamma(m-1).*hypergeom(m-1,1/2,1.98^2*a^2./(8*(a^2/2+b^2/2+m./gb))) - 1.98*a/sqrt(2)*gamma(m-1/2)*hypergeom(m-1/2,3/2,1.98^2*a^2./(8*(a^2/2+ b^2/2+m./gb))));
I2kara3 = (a^2/2+ b^2/2+m./gb).^(-m+1/2).*(sqrt(a^2/2+ b^2/2+m./gb)*gamma(m-1).*hypergeom(m-1,1/2,1.98^2*b^2./(8*(a^2/2+b^2/2+m./gb))) - 1.98*b/sqrt(2)*gamma(m-1/2)*hypergeom(m-1/2,3/2,1.98^2*b^2./(8*(a^2/2+ b^2/2+m./gb))));
I2kara4 = (a^2/2+ b^2/2+m./gb).^(-m+1/2).*(sqrt(a^2/2+ b^2/2+m./gb)*gamma(m-1).*hypergeom(m-1,1/2,1.98^2*(a+b)^2./(8*(a^2/2+b^2/2+m./gb))) - 1.98*(a+b)/sqrt(2)*gamma(m-1/2)*hypergeom(m-1/2,3/2,1.98^2*(a+b)^2./(8*(a^2/2+ b^2/2+m./gb))));
I2kara = const/(1.135^2*2*pi*a*b).*(I2kara1 - I2kara2 - I2kara3 + I2kara4);


I2Chiani = m^m./(gb.^m).*(1/144*(m./gb+a^2/2+ b^2/2).^(-m) + 1/48*(m./gb+2/3*a^2 + b^2/2).^(-m)+ 1/16*(m./gb+2*a^2/3+ 2*b^2/3).^(-m) + 1/48*(m./gb+2/3*b^2 + a^2/2).^(-m));

I2Olab2 = m^m./(4*gb.^m).*(0.3017^2*(m./gb+1.0510*(a^2+b^2)/2).^(-m) + 0.4389^2*(m./gb+2.102*(a^2+b^2)/2).^(-m)+ 0.3017*0.4389*((m./gb+1.0510*a^2/2+2.102*b^2/2).^(-m)+ (m./gb+1.0510*b^2/2+2.102*a^2/2).^(-m)));

I2Olab1 = m^m./(4*gb.^m).*(0.4803^2*(m./gb+1.1232*(a^2+b^2)/2).^(-m));

I2Prony2 = m^m./(4*gb.^m).*(0.416^2*(m./gb+1.942*(a^2+b^2)/2).^(-m) + 0.294^2*(m./gb+1.050*(a^2+b^2)/2).^(-m) + 0.416*0.294*((m./gb+1.942*b^2/2+1.050*a^2/2).^(-m)+ (m./gb+1.050*b^2/2+1.942*a^2/2).^(-m)));

I1Prony31 = m^m./(4*gb.^m).*(0.336^2*(m./gb+(a^2+b^2)*1.752/2).^(-m) + 0.288^2*(m./gb+1.050*(a^2+b^2)/2).^(-m) + 0.004^2*(m./gb+1.206*(a^2+b^2)/2).^(-m));
I1Prony32 = m^m./(4*gb.^m).*(0.336*0.288*(m./gb+1.752*a^2/2+1.050*b^2/2).^(-m) + 0.336*0.004*(m./gb+1.752*a^2/2+1.206*b^2/2).^(-m) + 0.288*0.004*(m./gb+1.050*a^2/2+1.206*b^2/2).^(-m));
I1Prony33 = m^m./(4*gb.^m).*(0.336*0.288*(m./gb+1.752*b^2/2+1.050*a^2/2).^(-m) + 0.336*0.004*(m./gb+1.752*b^2/2+1.206*a^2/2).^(-m)+ 0.288*0.004*(m./gb+1.050*b^2/2+1.206*a^2/2).^(-m));

I2Prony3 = I1Prony31 + I1Prony32 + I1Prony33;

I1Olab31 = m^m./(4*gb.^m).*(0.3357^2*(m./gb+(a^2+b^2)*1.0649/2).^(-m) + 0.3361^2*(m./gb+2.1298*(a^2+b^2)/2).^(-m) + 0.0305^2*(m./gb+3.1974*(a^2+b^2)/2).^(-m));
I1Olab32 = m^m./(4*gb.^m).*(0.3357*0.3361*(m./gb+1.0649*a^2/2+2.1298*b^2/2).^(-m) + 0.3357*0.0305*(m./gb+1.0649*a^2/2+3.1947*b^2/2).^(-m) + 0.3361*0.0305*(m./gb+2.1298*a^2/2+3.1947*b^2/2).^(-m));
I1Olab33 = m^m./(4*gb.^m).*(0.3357*0.3361*(m./gb+1.0649*b^2/2+2.1298*a^2/2).^(-m) + 0.3357*0.0305*(m./gb+1.0649*b^2/2+3.1947*a^2/2).^(-m)+ 0.3361*0.0305*(m./gb+2.1298*b^2/2+3.1947*a^2/2).^(-m));

I2Olab3 = I1Olab31 + I1Olab32 + I1Olab33;

I2mbfc1 = const*exp(-2*0.6966).*(m./gb+0.769*(a^2+b^2)/2).^(-1/2-m).*(sqrt(m./gb+0.769*(a^2+b^2)/2).*gamma(m).*hypergeom(m,1/2,1.0797^2*(a+b)^2./(8*(m./gb+0.769*(a^2+b^2)/2)))-1.0797*(a+b)/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,1.0797^2*(a+b)^2./(8*(m./gb+0.769*(a^2+b^2)/2))));

I2mbfc2 = const*exp(-2*0.9049).*(m./gb+0.9548*(a^2+b^2)/2).^(-1/2-m).*(sqrt(m./gb+0.9548*(a^2+b^2)/2).*gamma(m).*hypergeom(m,1/2,0.6341^2*(a+b)^2./(8*(m./gb+0.9548*(a^2+b^2)/2)))-0.6341*(a+b)/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,0.6341^2*(a+b)^2./(8*(m./gb+0.9548*(a^2+b^2)/2))));

I2ampj = const*0.25.*(m./gb+0.748*(a^2+b^2)/2).^(-1/2-m).*(sqrt(m./gb+0.748*(a^2+b^2)/2).*gamma(m).*hypergeom(m,1/2,1.0988^2*(a+b)^2./(8*(m./gb+0.748*(a^2+b^2)/2)))-1.0988*(a+b)/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,1.0988^2*(a+b)^2./(8*(m./gb+0.748*(a^2+b^2)/2))));

I2gh1 = const*0.6405^2/pi.*(m./gb+(a^2+b^2)/2).^(-1/2-m).*(sqrt(m./gb+(a^2+b^2)/2).*gamma(m).*hypergeom(m,1/2,0.6004^2*(a+b)^2./(8*(m./gb+(a^2+b^2)/2)))-0.6004*(a+b)/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,0.6004^2*(a+b)^2./(8*(m./gb+(a^2+b^2)/2))));

I2gh2 = const*0.2457^2/pi.*(m./gb+(a^2+b^2)/2).^(-1/2-m).*(sqrt(m./gb+(a^2+b^2)/2).*gamma(m).*hypergeom(m,1/2,2.5048^2*(a+b)^2./(8*(m./gb+(a^2+b^2)/2)))-2.5048*(a+b)/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,2.5048^2*(a+b)^2./(8*(m./gb+(a^2+b^2)/2))));

I2gh3 = const*0.2457*0.6405/pi.*(m./gb+(a^2+b^2)/2).^(-1/2-m).*(sqrt(m./gb+(a^2+b^2)/2).*gamma(m).*hypergeom(m,1/2,(2.5048*a+0.6004*b)^2./(8*(m./gb+(a^2+b^2)/2)))-(2.5048*a+0.6004*b)/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,(2.5048*a+0.6004*b)^2./(8*(m./gb+(a^2+b^2)/2))));
I2gh4 = const*0.2457*0.6405/pi.*(m./gb+(a^2+b^2)/2).^(-1/2-m).*(sqrt(m./gb+(a^2+b^2)/2).*gamma(m).*hypergeom(m,1/2,(2.5048*b+0.6004*a)^2./(8*(m./gb+(a^2+b^2)/2)))-(2.5048*b+0.6004*a)/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,(2.5048*b+0.6004*a)^2./(8*(m./gb+(a^2+b^2)/2))));

I2gh = I2gh1 + I2gh2 + I2gh3 + I2gh4;


I2SYA3 = m^m./(gb.^m).*(1/144*(m./gb+a^2/2+ b^2/2).^(-m) + 1/36*(m./gb+2*a^2 + 2*b^2).^(-m)+ 1/36*(m./gb+2*a^2/3+ 2*b^2/3).^(-m));
I2SYA3 = I2SYA3 + m^m./(gb.^m).*(1/72*(m./gb+2*b^2 + a^2/2).^(-m) + 1/72*(m./gb+b^2/2 + 2*a^2).^(-m) + 1/36*(m./gb+2*a^2/3+ 2*b^2).^(-m) + 1/36*(m./gb+2*a^2+ 2*b^2/3).^(-m));
I2SYA3 = I2SYA3 + m^m./(gb.^m).*(1/72*(m./gb+2*b^2/3 + a^2/2).^(-m) + 1/72*(m./gb+b^2/2 + 2*a^2/3).^(-m));

I2SYA4 =  m^m./(gb.^m).*(1/256*(m./gb+a^2/2+ b^2/2).^(-m) + 1/64*(m./gb+a^2 + b^2).^(-m)+ 1/64*(m./gb+10*a^2/3+ 10*b^2/3).^(-m) + 1/64*(m./gb+10*a^2/17+ 10*b^2/17).^(-m));
I2SYA4 = I2SYA4 + m^m./(gb.^m).*(1/128*(m./gb+b^2 + a^2/2).^(-m) + 1/128*(m./gb+10*b^2/3 + a^2/2).^(-m) + 1/128*(m./gb+a^2/2+ 10*b^2/17).^(-m));
I2SYA4 = I2SYA4 + m^m./(gb.^m).*(1/128*(m./gb+b^2/2 + a^2).^(-m) + 1/64*(m./gb+10*b^2/3 + a^2).^(-m) + 1/64*(m./gb+10*b^2/17 + a^2).^(-m));
I2SYA4 = I2SYA4 + m^m./(gb.^m).*(1/128*(m./gb+b^2/2 + 10*a^2/3).^(-m) + 1/64*(m./gb+b^2 + 10*a^2/3).^(-m) + 1/64*(m./gb+10*b^2/17 + 10*a^2/3).^(-m));
I2SYA4 = I2SYA4 + m^m./(gb.^m).*(1/128*(m./gb+b^2/2 + 10*a^2/17).^(-m) + 1/64*(m./gb+10*b^2/3 + 10*a^2/17).^(-m) + 1/64*(m./gb+b^2 + 10*a^2/17).^(-m));

sepexact = aq*I1exact + ai*I1exact1 - ai*aq*I2exact;
sepkara = aq*I1kara + ai*I1kara1 - ai*aq*I2kara;
sepprony2 = aq*I1Prony2 + ai*I1Prony21 - ai*aq*I2Prony2;
sepprony3 = aq*I1Prony3 + ai*I1Prony31 - ai*aq*I2Prony3;
sepolab1= aq*I1Olab1 + ai*I1Olab11 - ai*aq*I2Olab1;
sepolab2 = aq*I1Olab2 + ai*I1Olab21 - ai*aq*I2Olab2;
sepolab3 = aq*I1Olab3 + ai*I1Olab31 - ai*aq*I2Olab3;
sepgh = aq*I1gh + ai*I1gh1 - ai*aq*I2gh;
sepchiani = aq*I1Chiani + ai*I1Chiani1 - ai*aq*I2Chiani;
sepmbfc1 = aq*I1mbfc1 + ai*I1mbfc11 - ai*aq*I2mbfc1;
sepmbfc2  = aq*I1mbfc2 + ai*I1mbfc21 - ai*aq*I2mbfc2;
sepampj = aq*I1ampj + ai*I1ampj1 - ai*aq*I2ampj;
sepsya3 = aq*I1SYA3 + ai*I1SYA31 - ai*aq*I2SYA3;
sepsya4 = aq*I1SYA4 + ai*I1SYA41 - ai*aq*I2SYA4;

sepapprox = [sepkara; sepprony2; sepprony3; sepolab1; sepolab2; sepolab3; sepgh; sepchiani; sepmbfc1; sepmbfc2; sepampj; sepsya3; sepsya4];

relerr = zeros(13,nTerms);

for j = 1:13
    relerr(j,:) = abs(sepapprox(j,:) - sepexact)./sepexact;
end

figure, set(gca,'fontsize',14);
semilogy(snr, relerr(1,:), 'b', 'Linewidth', 2);
hold on
semilogy(snr, relerr(2,:), 'k', 'Linewidth', 2);
semilogy(snr, relerr(3,:), '--k', 'Linewidth', 2);
semilogy(snr, relerr(5,:), 'r', 'Linewidth', 2);
semilogy(snr, relerr(6,:), '--r', 'Linewidth', 2);
semilogy(snr, relerr(7,:), 'c', 'Linewidth', 2);
semilogy(snr, relerr(8,:), 'm', 'Linewidth', 2);
semilogy(snr, relerr(9,:), 'y', 'Linewidth', 2);
semilogy(snr, relerr(10,:), '--y', 'Linewidth', 2);
semilogy(snr, relerr(11,:), '--b', 'Linewidth', 2);
semilogy(snr, relerr(12,:), 'g', 'Linewidth', 2);
semilogy(snr, relerr(13,:), '--g', 'Linewidth', 2);
legend('E_{Kara}', 'E_{Prony}^2', 'E_{Prony}^3', 'E_{Olab}^2','E_{Olab}^3', 'E_{GH}^2','E_{Chiani}','E_{MBFC}^{SSE}', 'E_{MBFC}^{MARE}','E_{AMPJ}', 'E_{SYA3}', 'E_{SYA4}');

xlabel('Average SNR (dB)', 'FontSize', 14);
ylabel('Relative Error','FontSize', 14);

% figure, set(gca,'fontsize',14);
% semilogy(snr, sepexact, '--b', 'Linewidth', 2);
% hold on
% semilogy(snr, sepapprox(1,:), 'b', 'Linewidth', 2);
% semilogy(snr, sepapprox(2,:), 'k', 'Linewidth', 2);
% semilogy(snr, sepapprox(3,:), '--k', 'Linewidth', 2);
% semilogy(snr, sepapprox(5,:), 'r', 'Linewidth', 2);
% semilogy(snr, sepapprox(6,:), '--r', 'Linewidth', 2);
% semilogy(snr, sepapprox(7,:), 'c', 'Linewidth', 2);
% semilogy(snr, sepapprox(8,:), 'm', 'Linewidth', 2);
% semilogy(snr, sepapprox(9,:), 'y', 'Linewidth', 2);
% semilogy(snr, sepapprox(10,:), '--y', 'Linewidth', 2);
% semilogy(snr, sepapprox(11,:), '--b', 'Linewidth', 2);
% semilogy(snr, sepapprox(12,:), 'g', 'Linewidth', 2);
% semilogy(snr, sepapprox(13,:), '--g', 'Linewidth', 2);
% 
% legend('Exact', 'E_{Kara}', 'E_{Prony}^2', 'E_{Prony}^3', 'E_{Olab}^2','E_{Olab}^3', 'E_{GH}^2','E_{Chiani}','E_{MBFC}^{SSE}', 'E_{MBFC}^{MARE}','E_{AMPJ}', 'E_{SYA3}', 'E_{SYA4}');
% 
% xlabel('Average SNR (dB)', 'FontSize', 14);
% ylabel('Relative Error','FontSize', 14);
