
M = 32;
kmam = (M-1)/M;
a = sqrt(6/(M^2-1));
snr = -10:1:60;
gb = 10.^(snr/10);
nTerms = length(snr);

mt = 0.5:0.3:2;
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
    sepexact = 2*kmam*I1exact;

    % using kara
    I1Kara  = (a^2/2+ m./gb).^(-m+1/2)*gamma(m-1/2) - (a^2/2+m./gb).^(-m).*(sqrt(a^2/2+m./gb)*gamma(m-1/2).*hypergeom(m-1/2,1/2,1.98^2*a^2./(8*(a^2/2+m./gb))) - 1.98*a/sqrt(2)*gamma(m)*hypergeom(m,3/2,1.98^2*a^2./(8*(a^2/2+m./gb))));
    I1Kara = const/(1.135*sqrt(a*2*pi)) .* I1Kara;
    sepkara = 2*kmam*I1Kara;

    % using Chiani
    I1Chiani = m^m./(gb.^m).*(1/12*(m./gb+a^2/2).^(-m) + 1/4*(m./gb+2/3*a^2).^(-m));
    sepchiani = 2*kmam*I1Chiani;

    % using Prony2
    I1Prony2 = m^m./(2*gb.^m).*(0.416*(m./gb+1.942*a^2/2).^(-m) + 0.294*(m./gb+1.050*a^2/2).^(-m));
    sepprony2 = 2*kmam*I1Prony2;

    % using Prony3
    I1Prony3 = m^m./(2*gb.^m).*(0.336*(m./gb+a^2*1.752/2).^(-m) + 0.288*(m./gb+1.050*a^2/2).^(-m) + 0.004*(m./gb+1.206*a^2/2).^(-m));
    sepprony3 = 2*kmam*I1Prony3;

    % using Olab
    I1Olab1 = m^m./(2*gb.^m).*(0.4803*(m./gb+1.1232*a^2/2).^(-m));
    I1Olab2 = m^m./(2*gb.^m).*(0.3017*(m./gb+1.0510*a^2/2).^(-m) + 0.4389*(m./gb+2.102*a^2/2).^(-m));
    I1Olab3 = m^m./(2*gb.^m).*(0.3357*(m./gb+a^2*1.0649/2).^(-m) + 0.3361*(m./gb+2.1298*a^2/2).^(-m) + 0.0305*(m./gb+3.1947*a^2/2).^(-m));

    sepolab1 = 2*kmam*I1Olab1;
    sepolab2 = 2*kmam*I1Olab2;
    sepolab3 = 2*kmam*I1Olab3;

    I1gh = const*0.6405/sqrt(pi).*(m./gb+a^2/2).^(-1/2-m).*(sqrt(m./gb+a^2/2).*gamma(m).*hypergeom(m,1/2,0.6004^2*a^2./(8*(m./gb+a^2/2)))-0.6004*a/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,0.6004^2*a^2./(8*(m./gb+a^2/2))));

    I1gh = I1gh + const*0.2457/sqrt(pi).*(m./gb+a^2/2).^(-1/2-m).*(sqrt(m./gb+a^2/2).*gamma(m).*hypergeom(m,1/2,2.5048^2*a^2./(8*(m./gb+a^2/2)))-2.5048*a/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,2.5048^2*a^2./(8*(m./gb+a^2/2))));
    sepgh = 2*kmam*I1gh;

    I1SYA3 = m^m./(gb.^m) .*(1/12*(m./gb+a^2/2).^(-m) + 1/6*(m./gb+2*a^2).^(-m) + 1/6*(m./gb+2/3*a^2).^(-m));
    sepsya3 = 2*kmam*I1SYA3;

    I1SYA4 = m^m./(gb.^m) .*(1/16*(m./gb+a^2/2).^(-m) + 1/8*(m./gb+a^2).^(-m) + 1/8*(m./gb+10/3*a^2).^(-m) + + 1/8*(m./gb+10/17*a^2).^(-m));
    sepsya4 = 2*kmam*I1SYA4;

    I1mbfc1 = const*exp(-0.6966).*(m./gb+0.769*a^2/2).^(-1/2-m).*(sqrt(m./gb+0.769*a^2/2).*gamma(m).*hypergeom(m,1/2,1.0797^2*a^2./(8*(m./gb+0.769*a^2/2)))-1.0797*a/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,1.0797^2*a^2./(8*(m./gb+0.769*a^2/2))));
    sepmbfc1 = 2*kmam*I1mbfc1;

    I1mbfc2 = const*exp(-0.9049).*(m./gb+0.9548*a^2/2).^(-1/2-m).*(sqrt(m./gb+0.9548*a^2/2).*gamma(m).*hypergeom(m,1/2,0.6341^2*a^2./(8*(m./gb+0.9548*a^2/2)))-0.6341*a/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,0.6341^2*a^2./(8*(m./gb+0.9548*a^2/2))));
    sepmbfc2 = 2*kmam*I1mbfc2;

    I1ampj = const*0.5.*(m./gb+0.748*a^2/2).^(-1/2-m).*(sqrt(m./gb+0.748*a^2/2).*gamma(m).*hypergeom(m,1/2,1.0988^2*a^2./(8*(m./gb+0.748*a^2/2)))-1.0988*a/sqrt(2)*gamma(m+1/2).*hypergeom(m+1/2,3/2,1.0988^2*a^2./(8*(m./gb+0.748*a^2/2))));
    sepampj = 2*kmam*I1ampj;

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
ylabel('Relative Error','FontSize', 14);
