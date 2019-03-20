M = 16;
N = 8;
m = 7;

Nn = 4- 2*(M+N)/(M*N);
a = sqrt(96/(31*M*N - 32));
Kx = 4 - 4*(M+N)/(M*N) + 8/(M*N);

snr = 45;
gb = 10.^(snr/10);
nTerms = length(snr);
mt = 2.1:0.3:10;
relterms = length(mt);
avgRel = zeros(13, relterms);


for ft = 1:relterms
    m = mt(ft);
    const = m^m./(gb.^m * gamma(m));
    I1exact = zeros(1,nTerms);
    for t = 1:nTerms
        func = @(x)(erfc(a*sqrt(x)/sqrt(2)).*x.^(m-1).*exp(-m/gb(t).*x));
        I1exact(t) = 0.5*const(t)*integral(func, 0, Inf);
    end

    I2exact = zeros(1,nTerms);
    k = 2;
    for t = 1:nTerms
        func = @(x)((erfc(a*sqrt(x)/sqrt(2))).^k.*x.^(m-1).*exp(-m/gb(t).*x));
        I2exact(t) = 1/2^k*const(t)*integral(func, 0, Inf);
    end

    % using kara
    I1kara = (a^2/2+ m./gb).^(-m+1/2)*gamma(m-1/2) - (a^2/2+m./gb).^(-m).*(sqrt(a^2/2+m./gb)*gamma(m-1/2).*hypergeom(m-1/2,1/2,1.98^2*a^2./(8*(a^2/2+m./gb))) - 1.98*a/sqrt(2)*gamma(m)*hypergeom(m,3/2,1.98^2*a^2./(8*(a^2/2+m./gb))));
    I1kara = const/(1.135*sqrt(2*pi)*a).*I1kara;


    k = 2;
    I2kara = zeros(1,nTerms);

    for i = 0:k
        I2kara = I2kara + (-1)^i*nchoosek(k,i)*((k*a^2/2+m./gb).^(-m+k/2-1/2).*(sqrt(k*a^2/2+m./gb)*gamma(m-k/2).*hypergeom(m-k/2,1/2,1.98^2*a^2*i^2./(8*(k*a^2/2+m./gb))) - i*1.98*a/sqrt(2)*gamma(m-k/2+1/2)*hypergeom(m-k/2+1/2,3/2,i^2*1.98^2*a^2./(8*(k*a^2/2+m./gb)))));
    end

    I2kara = const/(1.135*a*sqrt(2*pi))^k.*I2kara;

    I1Prony2 = m^m./(2*gb.^m).*(0.416*(m./gb+1.942*a^2/2).^(-m) + 0.294*(m./gb+1.050*a^2/2).^(-m));

    I2Prony2 = zeros(1,nTerms);
    k = 2;
    for i =0:k
        I2Prony2 = I2Prony2 + nchoosek(k,i)*(0.416/2)^i*(0.294/2)^(k-i)*(m./gb + 1.942*a^2*i/2 + 1.050*a^2*(k-i)/2).^(-m);
    end

    I2Prony2 = m^m./(gb.^m).*I2Prony2;

    I1Prony3 = m^m./(2*gb.^m).*(0.336*(m./gb+a^2*1.752/2).^(-m) + 0.288*(m./gb+1.050*a^2/2).^(-m) + 0.004*(m./gb+1.206*a^2/2).^(-m));

    I2Prony3 = zeros(1,nTerms);
    k = 2;
    for i =0:k
        for j=0:(k-i)
            I2Prony3 = I2Prony3 + nchoosek(k,i)*nchoosek(k-i,j)*(0.336/2)^i*(0.288/2)^j*(0.004/2)^(k-i-j)*(1.752*a^2*i/2 + 1.050*a^2*j/2 + 1.206*a^2*(k-i-j)/2 + m./gb).^(-m);
        end
    end
    I2Prony3 = m^m./(gb.^m).*I2Prony3; 


    I1Olab1 = m^m./(2*gb.^m).*(0.4803*(m./gb+1.1232*a^2/2).^(-m));

    k = 2;
    I2Olab1 = m^m./(2^k*gb.^m).*(0.4803^k*(m./gb+1.1232*a^2*k/2).^(-m));

    I1Olab2 = m^m./(2*gb.^m).*(0.3017*(m./gb+1.0510*a^2/2).^(-m) + 0.4389*(m./gb+2.102*a^2/2).^(-m));

    I2Olab2 = zeros(1,nTerms);
    k =2;
    for i =0:k
        I2Olab2 = I2Olab2 + nchoosek(k,i)*(0.3017/2)^i*(0.4389/2)^(k-i)*(m./gb + 1.0510*a^2*i/2 + 2.102*a^2*(k-i)/2).^(-m);
    end

    I2Olab2 = m^m./(gb.^m).*I2Olab2;

    I1Olab3 = m^m./(2*gb.^m).*(0.3357*(m./gb+a^2*1.0649/2).^(-m) + 0.3361*(m./gb+2.1298*a^2/2).^(-m) + 0.0305*(m./gb+3.1947*a^2/2).^(-m));

    I2Olab3 = zeros(1,nTerms);
    k =2;
    for i =0:k
        for j=0:(k-i)
            I2Olab3 = I2Olab3 + nchoosek(k,i)*nchoosek(k-i,j)*(0.3357/2)^i*(0.3361/2)^j*(0.0305/2)^(k-i-j)*(1.6049*a^2*i/2 + 2.1298*a^2*j/2 + 3.1947*a^2*(k-i-j)/2 + m./gb).^(-m);
        end
    end
    I2Olab3 = m^m./(gb.^m).*I2Olab3;


    I1gh = const*0.6405/sqrt(pi).*(m./gb+a^2/2).^(-1/2-m).*(sqrt(m./gb+a^2/2).*gamma(m).*hypergeom(m,1/2,0.6004^2*a^2./(8*(m./gb+a^2/2)))-0.6004*a/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,0.6004^2*a^2./(8*(m./gb+a^2/2))));

    I1gh = I1gh + const*0.2457/sqrt(pi).*(m./gb+a^2/2).^(-1/2-m).*(sqrt(m./gb+a^2/2).*gamma(m).*hypergeom(m,1/2,2.5048^2*a^2./(8*(m./gb+a^2/2)))-2.5048*a/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,2.5048^2*a^2./(8*(m./gb+a^2/2))));

    I2gh = zeros(1,nTerms);
    k=2;
    for i = 0:k
        I2gh = I2gh + nchoosek(k,i)*0.6405^i*0.2457^(k-i)*((m./gb+(a^2*k)/2).^(-1/2-m).*(sqrt(m./gb+(a^2*k)/2).*gamma(m).*hypergeom(m,1/2,(0.6004*i+2.5048*(k-i))^2*a^2./(8*(m./gb+(a^2*k)/2)))-(0.6004*i+2.5048*(k-i))*a/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,(0.6004*i+2.5048*(k-i))^2*a^2./(8*(m./gb+(a^2*k)/2)))));
    end

    I2gh = const/(sqrt(pi))^k .*I2gh;

    I1Chiani = m^m./(gb.^m).*(1/12*(m./gb+a^2/2).^(-m) + 1/4*(m./gb+2/3*a^2).^(-m));

    I2Chiani = zeros(1,nTerms);
    k=2;
    for i =0:k
        I2Chiani = I2Chiani + nchoosek(k,i)*(1/12)^i*(1/4)^(k-i)*(m./gb+a^2*2*k/3-a^2*i/6).^(-m);
    end

    I2Chiani = m^m./(gb.^m).*I2Chiani;

    % using SYA3
    I1SYA3 = m^m./(gb.^m) .*(1/12*(m./gb+a^2/2).^(-m) + 1/6*(m./gb+2*a^2).^(-m) + 1/6*(m./gb+2/3*a^2).^(-m));

    I2SYA3 = zeros(1,nTerms);

    k=2;
    for i =0:k
        for j = 0:(k-i)
            I2SYA3 = I2SYA3 + 1/12^i*1/6^(k-i)*nchoosek(k,i)*nchoosek(k-i,j).*(i*a^2/2 + 2*j*a^2 + 2*(k-i-j)*a^2/3 + m./gb).^(-m);
        end
    end
    I2SYA3 = m^m./(gb.^m) .* I2SYA3;


    % using SYA4
    I1SYA4 = m^m./(gb.^m) .*(1/16*(m./gb+a^2/2).^(-m) + 1/8*(m./gb+a^2).^(-m) + 1/8*(m./gb+10/3*a^2).^(-m) + + 1/8*(m./gb+10/17*a^2).^(-m));

    I2SYA4 = zeros(1,nTerms);

    k=2;
    for i =0:k
        for j = 0:(k-i)
            for r = 0:(k-i-j)
                I2SYA4 = I2SYA4 + 1/16^i*1/8^(k-i)*nchoosek(k,i)*nchoosek(k-i,j)*nchoosek(k-i-j,r).*(i*a^2/2 + j*a^2 + 10*r*a^2/3 + 10*(k-i-j-r)*a^2/17 + m./gb).^(-m);
            end
        end
    end
    I2SYA4 = m^m./(gb.^m) .* I2SYA4;


    I1mbfc1 = const*exp(-0.6966).*(m./gb+0.769*a^2/2).^(-1/2-m).*(sqrt(m./gb+0.769*a^2/2).*gamma(m).*hypergeom(m,1/2,1.0797^2*a^2./(8*(m./gb+0.769*a^2/2)))-1.0797*a/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,1.0797^2*a^2./(8*(m./gb+0.769*a^2/2))));

    k=2;
    I2mbfc1 = const*exp(-k*0.6966).*(m./gb+0.769*(a^2*k)/2).^(-1/2-m).*(sqrt(m./gb+0.769*(a^2*k)/2).*gamma(m).*hypergeom(m,1/2,1.0797^2*(a*k)^2./(8*(m./gb+0.769*(a^2*k)/2)))-1.0797*(a*k)/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,1.0797^2*(a*k)^2./(8*(m./gb+0.769*(a^2*k)/2))));

    I1mbfc2 = const*exp(-0.9049).*(m./gb+0.9548*a^2/2).^(-1/2-m).*(sqrt(m./gb+0.9548*a^2/2).*gamma(m).*hypergeom(m,1/2,0.6341^2*a^2./(8*(m./gb+0.9548*a^2/2)))-0.6341*a/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,0.6341^2*a^2./(8*(m./gb+0.9548*a^2/2))));

    k=2;
    I2mbfc2 = const*exp(-k*0.9049).*(m./gb+0.9548*(a^2*k)/2).^(-1/2-m).*(sqrt(m./gb+0.9548*(a^2*k)/2).*gamma(m).*hypergeom(m,1/2,0.6341^2*(a*k)^2./(8*(m./gb+0.9548*(a^2*k)/2)))-0.6341*(a*k)/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,0.6341^2*(a*k)^2./(8*(m./gb+0.9548*(a^2*k)/2))));


    I1ampj = const*0.5.*(m./gb+0.748*a^2/2).^(-1/2-m).*(sqrt(m./gb+0.748*a^2/2).*gamma(m).*hypergeom(m,1/2,1.0988^2*a^2./(8*(m./gb+0.748*a^2/2)))-1.0988*a/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,1.0988^2*a^2./(8*(m./gb+0.748*a^2/2))));

    k=2;
    I2ampj = const*(0.5)^k.*(m./gb+0.748*(a^2*k)/2).^(-1/2-m).*(sqrt(m./gb+0.748*(a^2*k)/2).*gamma(m).*hypergeom(m,1/2,1.0988^2*(a*k)^2./(8*(m./gb+0.748*(a^2*k)/2)))-1.0988*(a*k)/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,1.0988^2*(a*k)^2./(8*(m./gb+0.748*(a^2*k)/2))));

    sepexact = Nn*I1exact - Kx*I2exact;
    sepkara = Nn*I1kara - Kx*I2kara;
    sepprony2 = Nn*I1Prony2 - Kx*I2Prony2;
    sepprony3 = Nn*I1Prony3 - Kx*I2Prony3;
    sepolab1= Nn*I1Olab1 - Kx*I2Olab1;
    sepolab2 = Nn*I1Olab2 - Kx*I2Olab2;
    sepolab3 = Nn*I1Olab3 - Kx*I2Olab3;
    sepgh = Nn*I1gh - Kx*I2gh;
    sepchiani = Nn*I1Chiani - Kx*I2Chiani;
    sepmbfc1 = Nn*I1mbfc1 - Kx*I2mbfc1;
    sepmbfc2  = Nn*I1mbfc2 - Kx*I2mbfc2;
    sepampj = Nn*I1ampj - Kx*I2ampj;
    sepsya3 = Nn*I1SYA3 - Kx*I2SYA3;
    sepsya4 = Nn*I1SYA4 - Kx*I2SYA4;

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


