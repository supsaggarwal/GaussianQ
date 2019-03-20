% sep for HQAM

% m = 0.9;
% no of points in HQAM

M =64;
alpha = 8/141;
K = 163/32;
Kc = 75/16;

% M = 256;
% alpha = 24/(7*M-4);
% K = 2*(3 - 4*M^(-1/2) + 1/M);
% Kc = 6*(1 - M^(-1/2))^2;

% Q(a)*pdf
snr = 45;
gb = 10.^(snr/10);
nTerms = length(snr);
mt = 2.1:0.3:10;
relterms = length(mt);
avgRel = zeros(13, relterms);


for ft = 1:relterms
    m = mt(ft);
    const = m^m./(gb.^m * gamma(m));

    a = sqrt(alpha);
    I1exact = zeros(1,nTerms);
    for t = 1:nTerms
        func = @(x)(erfc(a*sqrt(x)/sqrt(2)).*x.^(m-1).*exp(-m/gb(t).*x));
        I1exact(t) = 0.5*const(t)*integral(func, 0, Inf);
    end


    I1kara = (a^2/2+ m./gb).^(-m+1/2)*gamma(m-1/2) - (a^2/2+m./gb).^(-m).*(sqrt(a^2/2+m./gb)*gamma(m-1/2).*hypergeom(m-1/2,1/2,1.98^2*a^2./(8*(a^2/2+m./gb))) - 1.98*a/sqrt(2)*gamma(m)*hypergeom(m,3/2,1.98^2*a^2./(8*(a^2/2+m./gb))));
    I1kara = const/(1.135*sqrt(2*pi)*a).*I1kara;

    I1Chiani = m^m./(gb.^m).*(1/12*(m./gb+a^2/2).^(-m) + 1/4*(m./gb+2/3*a^2).^(-m));

    I1Olab3 = m^m./(2*gb.^m).*(0.3357*(m./gb+a^2*1.0649/2).^(-m) + 0.3361*(m./gb+2.1298*a^2/2).^(-m) + 0.0305*(m./gb+3.1947*a^2/2).^(-m));

    I1Olab2 = m^m./(2*gb.^m).*(0.3017*(m./gb+1.0510*a^2/2).^(-m) + 0.4389*(m./gb+2.102*a^2/2).^(-m));

    I1Olab1 = m^m./(2*gb.^m).*(0.4803*(m./gb+1.1232*a^2/2).^(-m));

    I1Prony2 = m^m./(2*gb.^m).*(0.416*(m./gb+1.942*a^2/2).^(-m) + 0.294*(m./gb+1.050*a^2/2).^(-m));
    I1Prony3 = m^m./(2*gb.^m).*(0.336*(m./gb+a^2*1.752/2).^(-m) + 0.288*(m./gb+1.050*a^2/2).^(-m) + 0.004*(m./gb+1.206*a^2/2).^(-m));

    I1Isuk = zeros(1,nTerms);
    for i=1:8
        I1Isuk = I1Isuk + 1.98^i*(-1)^(i-1)/factorial(i)*gamma(m+i/2-1/2)*(m./gb+a^2/2).^(-m-i/2+1/2)*a^(i-1)*2^(-i/2+1/2);
    end
    I1Isuk = const/(1.135*2*sqrt(pi)).*I1Isuk;

    I1mbfc1 = const*exp(-0.6966).*(m./gb+0.769*a^2/2).^(-1/2-m).*(sqrt(m./gb+0.769*a^2/2).*gamma(m).*hypergeom(m,1/2,1.0797^2*a^2./(8*(m./gb+0.769*a^2/2)))-1.0797*a/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,1.0797^2*a^2./(8*(m./gb+0.769*a^2/2))));

    I1mbfc2 = const*exp(-0.9049).*(m./gb+0.9548*a^2/2).^(-1/2-m).*(sqrt(m./gb+0.9548*a^2/2).*gamma(m).*hypergeom(m,1/2,0.6341^2*a^2./(8*(m./gb+0.9548*a^2/2)))-0.6341*a/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,0.6341^2*a^2./(8*(m./gb+0.9548*a^2/2))));

    I1ampj = const*0.5.*(m./gb+0.748*a^2/2).^(-1/2-m).*(sqrt(m./gb+0.748*a^2/2).*gamma(m).*hypergeom(m,1/2,1.0988^2*a^2./(8*(m./gb+0.748*a^2/2)))-1.0988*a/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,1.0988^2*a^2./(8*(m./gb+0.748*a^2/2))));

    I1gh = const*0.6405/sqrt(pi).*(m./gb+a^2/2).^(-1/2-m).*(sqrt(m./gb+a^2/2).*gamma(m).*hypergeom(m,1/2,0.6004^2*a^2./(8*(m./gb+a^2/2)))-0.6004*a/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,0.6004^2*a^2./(8*(m./gb+a^2/2))));

    I1gh = I1gh + const*0.2457/sqrt(pi).*(m./gb+a^2/2).^(-1/2-m).*(sqrt(m./gb+a^2/2).*gamma(m).*hypergeom(m,1/2,2.5048^2*a^2./(8*(m./gb+a^2/2)))-2.5048*a/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,2.5048^2*a^2./(8*(m./gb+a^2/2))));

    I1SYA3 = m^m./(gb.^m) .*(1/12*(m./gb+a^2/2).^(-m) + 1/6*(m./gb+2*a^2).^(-m) + 1/6*(m./gb+2/3*a^2).^(-m));
    I1SYA4 = m^m./(gb.^m) .*(1/16*(m./gb+a^2/2).^(-m) + 1/8*(m./gb+a^2).^(-m) + 1/8*(m./gb+10/3*a^2).^(-m) + + 1/8*(m./gb+10/17*a^2).^(-m));





    % Q(a)Q(b)PDF

    a = sqrt(alpha);
    b = sqrt(alpha/3);
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

    I2Isuk = zeros(1,nTerms);
    for i=1:8
        for j = 1:8
            I2Isuk = I2Isuk + 1.98^(i+j)*(-1)^(i+j-2)/(factorial(i)*factorial(j))*gamma(m+i/2+j/2-1)*(m./gb+(a^2+b^2)/2).^(-m-i/2-j/2+1)*a^(i-1)*b^(j-1)*2^(-i/2-j/2+1);
        end
    end
    I2Isuk = const/(1.135^2*4*pi).*I2Isuk;

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



    %Q^2(a)PDF
    a = sqrt(2*alpha/3);
    k = 2;
    b = a;
    I3Isuk = zeros(1,nTerms);
    for i=1:8
        for j = 1:8
            I3Isuk = I3Isuk + 1.98^(i+j)*(-1)^(i+j-2)/(factorial(i)*factorial(j))*gamma(m+i/2+j/2-1)*(m./gb+(a^2+b^2)/2).^(-m-i/2-j/2+1)*a^(i-1)*b^(j-1)*2^(-i/2-j/2+1);
        end
    end
    I3Isuk = const/(1.135^2*4*pi).*I3Isuk;

    I3exact = zeros(1,nTerms);
    for t = 1:nTerms
        func = @(x)((erfc(a*sqrt(x)/sqrt(2))).^k.*x.^(m-1).*exp(-m/gb(t).*x));
        I3exact(t) = 1/2^k*const(t)*integral(func, 0, Inf);
    end



    I3kara = zeros(1,nTerms);

    for i = 0:k
        I3kara = I3kara + (-1)^i*nchoosek(k,i)*((k*a^2/2+m./gb).^(-m+k/2-1/2).*(sqrt(k*a^2/2+m./gb)*gamma(m-k/2).*hypergeom(m-k/2,1/2,1.98^2*a^2*i^2./(8*(k*a^2/2+m./gb))) - i*1.98*a/sqrt(2)*gamma(m-k/2+1/2)*hypergeom(m-k/2+1/2,3/2,i^2*1.98^2*a^2./(8*(k*a^2/2+m./gb)))));
    end

    I3kara = const/(1.135*a*sqrt(2*pi))^k.*I3kara;

    % I3kara1 = (3*a^2/2+ m./gb).^(-m+3/2)*gamma(m-3/2);
    % I3kara2 = (3*a^2/2+m./gb).^(-m+1).*(sqrt(3*a^2/2+m./gb)*gamma(m-3/2).*hypergeom(m-3/2,1/2,9*1.98^2*a^2./(8*(3*a^2/2+m./gb))) - 3*1.98*a/sqrt(2)*gamma(m-1)*hypergeom(m-1,3/2,9*1.98^2*a^2./(8*(3*a^2/2+m./gb))));
    % I3kara3 = (3*a^2/2+m./gb).^(-m+1).*(sqrt(3*a^2/2+m./gb)*gamma(m-3/2).*hypergeom(m-3/2,1/2,1.98^2*a^2./(8*(3*a^2/2+m./gb))) - 1.98*a/sqrt(2)*gamma(m-1)*hypergeom(m-1,3/2,1.98^2*a^2./(8*(3*a^2/2+m./gb))));
    % I3kara4 = (3*a^2/2+m./gb).^(-m+1).*(sqrt(3*a^2/2+m./gb)*gamma(m-3/2).*hypergeom(m-3/2,1/2,4*1.98^2*a^2./(8*(3*a^2/2+m./gb))) - 2*1.98*a/sqrt(2)*gamma(m-1)*hypergeom(m-1,3/2,4*1.98^2*a^2./(8*(3*a^2/2+m./gb))));
    % 
    % I3kara = const/(1.135^3*(sqrt(2*pi))^3*a^3).*(I3kara1 - I3kara2 - 3*I3kara3 + 3*I3kara4);

    I3Chiani = zeros(1,nTerms);
    for i =0:k
        I3Chiani = I3Chiani + nchoosek(k,i)*(1/12)^i*(1/4)^(k-i)*(m./gb+a^2*2*k/3-a^2*i/6).^(-m);
    end

    I3Chiani = m^m./(gb.^m).*I3Chiani;


    I3Prony2 = zeros(1,nTerms);
    for i =0:k
        I3Prony2 = I3Prony2 + nchoosek(k,i)*(0.416/2)^i*(0.294/2)^(k-i)*(m./gb + 1.942*a^2*i/2 + 1.050*a^2*(k-i)/2).^(-m);
    end

    I3Prony2 = m^m./(gb.^m).*I3Prony2;


    I3Olab2 = zeros(1,nTerms);
    for i =0:k
        I3Olab2 = I3Olab2 + nchoosek(k,i)*(0.3017/2)^i*(0.4389/2)^(k-i)*(m./gb + 1.0510*a^2*i/2 + 2.102*a^2*(k-i)/2).^(-m);
    end

    I3Olab2 = m^m./(gb.^m).*I3Olab2;

    I3Olab1 = m^m./(2^k*gb.^m).*(0.4803^k*(m./gb+1.1232*a^2*k/2).^(-m));

    I3Prony3 = zeros(1,nTerms);
    for i =0:k
        for j=0:(k-i)
            I3Prony3 = I3Prony3 + nchoosek(k,i)*nchoosek(k-i,j)*(0.336/2)^i*(0.288/2)^j*(0.004/2)^(k-i-j)*(1.752*a^2*i/2 + 1.050*a^2*j/2 + 1.206*a^2*(k-i-j)/2 + m./gb).^(-m);
        end
    end
    I3Prony3 = m^m./(gb.^m).*I3Prony3; 



    I3Olab3 = zeros(1,nTerms);
    for i =0:k
        for j=0:(k-i)
            I3Olab3 = I3Olab3 + nchoosek(k,i)*nchoosek(k-i,j)*(0.3357/2)^i*(0.3361/2)^j*(0.0305/2)^(k-i-j)*(1.6049*a^2*i/2 + 2.1298*a^2*j/2 + 3.1947*a^2*(k-i-j)/2 + m./gb).^(-m);
        end
    end
    I3Olab3 = m^m./(gb.^m).*I3Olab3;

    I3mbfc1 = const*exp(-k*0.6966).*(m./gb+0.769*(a^2*k)/2).^(-1/2-m).*(sqrt(m./gb+0.769*(a^2*k)/2).*gamma(m).*hypergeom(m,1/2,1.0797^2*(a*k)^2./(8*(m./gb+0.769*(a^2*k)/2)))-1.0797*(a*k)/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,1.0797^2*(a*k)^2./(8*(m./gb+0.769*(a^2*k)/2))));

    I3mbfc2 = const*exp(-k*0.9049).*(m./gb+0.9548*(a^2*k)/2).^(-1/2-m).*(sqrt(m./gb+0.9548*(a^2*k)/2).*gamma(m).*hypergeom(m,1/2,0.6341^2*(a*k)^2./(8*(m./gb+0.9548*(a^2*k)/2)))-0.6341*(a*k)/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,0.6341^2*(a*k)^2./(8*(m./gb+0.9548*(a^2*k)/2))));

    I3ampj = const*(0.5)^k.*(m./gb+0.748*(a^2*k)/2).^(-1/2-m).*(sqrt(m./gb+0.748*(a^2*k)/2).*gamma(m).*hypergeom(m,1/2,1.0988^2*(a*k)^2./(8*(m./gb+0.748*(a^2*k)/2)))-1.0988*(a*k)/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,1.0988^2*(a*k)^2./(8*(m./gb+0.748*(a^2*k)/2))));

    I3gh = zeros(1,nTerms);
    for i = 0:k
        I3gh = I3gh + nchoosek(k,i)*0.6405^i*0.2457^(k-i)*((m./gb+(a^2*k)/2).^(-1/2-m).*(sqrt(m./gb+(a^2*k)/2).*gamma(m).*hypergeom(m,1/2,(0.6004*i+2.5048*(k-i))^2*a^2./(8*(m./gb+(a^2*k)/2)))-(0.6004*i+2.5048*(k-i))*a/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,(0.6004*i+2.5048*(k-i))^2*a^2./(8*(m./gb+(a^2*k)/2)))));
    end

    I3gh = const/(sqrt(pi))^k .*I3gh;


    I3SYA3 = zeros(1,nTerms);

    for i =0:k
        for j = 0:(k-i)
            I3SYA3 = I3SYA3 + 1/12^i*1/6^(k-i)*nchoosek(k,i)*nchoosek(k-i,j).*(i*a^2/2 + 2*j*a^2 + 2*(k-i-j)*a^2/3 + m./gb).^(-m);
        end
    end
    I3SYA3 = m^m./(gb.^m) .* I3SYA3;

    I3SYA4 = zeros(1,nTerms);

    for i =0:k
        for j = 0:(k-i)
            for r = 0:(k-i-j)
                I3SYA4 = I3SYA4 + 1/16^i*1/8^(k-i)*nchoosek(k,i)*nchoosek(k-i,j)*nchoosek(k-i-j,r).*(i*a^2/2 + j*a^2 + 10*r*a^2/3 + 10*(k-i-j-r)*a^2/17 + m./gb).^(-m);
            end
        end
    end
    I3SYA4 = m^m./(gb.^m) .* I3SYA4;


    sepexact = K*I1exact + 2/3*Kc*I3exact -2*Kc*I2exact;

    sepkara = K*I1kara + 2/3*Kc*I3kara -2*Kc*I2kara;

    sepIsuk = K*I1Isuk + 2/3*Kc*I3Isuk -2*Kc*I2Isuk;

    sepprony2 = K*I1Prony2 + 2/3*Kc*I3Prony2 -2*Kc*I2Prony2;

    sepprony3 = K*I1Prony3 + 2/3*Kc*I3Prony3 -2*Kc*I2Prony3;

    sepolab1 = K*I1Olab1 + 2/3*Kc*I3Olab1 -2*Kc*I2Olab1;

    sepolab2 = K*I1Olab2 + 2/3*Kc*I3Olab2 -2*Kc*I2Olab2;

    sepolab3 = K*I1Olab3 + 2/3*Kc*I3Olab3 -2*Kc*I2Olab3;

    sepgh = K*I1gh + 2/3*Kc*I3gh -2*Kc*I2gh;

    sepchiani = K*I1Chiani + 2/3*Kc*I3Chiani -2*Kc*I2Chiani;

    sepmbfc1 = K*I1mbfc1 + 2/3*Kc*I3mbfc1 -2*Kc*I2mbfc1;

    sepmbfc2 = K*I1mbfc2 + 2/3*Kc*I3mbfc2 -2*Kc*I2mbfc2;

    sepampj = K*I1ampj + 2/3*Kc*I3ampj -2*Kc*I2ampj;

    sepsya3 = K*I1SYA3 + 2/3*Kc*I3SYA3 -2*Kc*I2SYA3;

    sepsya4 = K*I1SYA4 + 2/3*Kc*I3SYA4 -2*Kc*I2SYA4;

    sepapprox = [sepkara; sepprony2; sepprony3; sepolab1; sepolab2; sepolab3; sepgh; sepchiani; sepmbfc1; sepmbfc2; sepampj; sepsya3; sepsya4];

    relerr = zeros(13,nTerms);

    for j = 1:13
        relerr(j,:) = abs(sepapprox(j,:) - sepexact)./sepexact;
    end

avgRel(:,ft) = sum(relerr,2)/nTerms;
end
figure, set(gca,'fontsize',14);
semilogy(mt, avgRel(1,:), 'b', 'Linewidth', 2);
hold on
semilogy(mt, avgRel(2,:), 'k', 'Linewidth', 2);
semilogy(mt, avgRel(3,:), '--k', 'Linewidth', 2);
semilogy(mt, avgRel(5,:), 'r', 'Linewidth', 2);
semilogy(mt, avgRel(6,:), '--r', 'Linewidth', 2);
semilogy(mt, avgRel(7,:), 'c', 'Linewidth', 2);
semilogy(mt, avgRel(8,:), 'm', 'Linewidth', 2);
semilogy(mt, avgRel(9,:), 'y', 'Linewidth', 2);
semilogy(mt, avgRel(10,:), '--y', 'Linewidth', 2);
semilogy(mt, avgRel(11,:), '--b', 'Linewidth', 2);
semilogy(mt, avgRel(12,:), 'g', 'Linewidth', 2);
semilogy(mt, avgRel(13,:), '--g', 'Linewidth', 2);


legend('E_{Kara}', 'E_{Prony}^2', 'E_{Prony}^3', 'E_{Olab}^2','E_{Olab}^3', 'E_{GH}^2','E_{Chiani}','E_{MBFC}^{SSE}', 'E_{MBFC}^{MARE}','E_{AMPJ}', 'E_{SYA3}', 'E_{SYA4}');

xlabel('Fading Parameter m', 'FontSize', 14);
ylabel('Average Relative Error','FontSize', 14);











