% all existing approximations calculation of erfc(x)

x = 0:0.01:7;

nTerms = length(x);

exact = 0.5*erfc(x/sqrt(2));
% exact =erfc(x);
eBorj = 1/sqrt(pi) .* exp(-x.^2)./(.661*x + .339.*(x.^2+2.755));

eKara = 1/(sqrt(pi)*1.135).*(1- exp(-1.98*x)).*exp(-x.^2)./x;
eKara(1) = 0.984;

eIsuk = zeros(1,nTerms);
for i = 1:8
    eIsuk = eIsuk + (-1)^(i-1)*1.98^i/factorial(i)* x.^(i-1);
end
    
eIsuk = eIsuk.*exp(-x.^2)./(1.135*sqrt(pi));

% eProny2 = 0.416*exp(-1.942*x.^2) + 0.294*exp(-1.050*x.^2);
% eProny3 = 0.336*exp(-1.752*x.^2) + 0.288*exp(-1.050*x.^2) + 0.004*exp(-1.206*x.^2);
% 
% eOlab1 = 0.4803*exp(-1.1232*x.^2);
% eOlab2 = 0.3017*exp(-1.0510*x.^2) + 0.4389*exp(-2.102*x.^2);
% eOlab3 = 0.3357*exp(-1.0649*x.^2) + 0.3361*exp(-2.1298*x.^2) + 0.0305*exp(-3.1947*x.^2);

eProny2 = (0.416*exp(-1.942*x.^2/2) + 0.294*exp(-1.050*x.^2/2))/2;
eProny3 = (0.336*exp(-1.752*x.^2/2) + 0.288*exp(-1.050*x.^2/2) + 0.004*exp(-1.206*x.^2/2))/2;

eOlab1 = (0.4803*exp(-1.1232*x.^2/2))/2;
eOlab2 = (0.3017*exp(-1.0510*x.^2/2) + 0.4389*exp(-2.102*x.^2/2))/2;
eOlab3 = (0.3357*exp(-1.0649*x.^2/2) + 0.3361*exp(-2.1298*x.^2/2) + 0.0305*exp(-3.1947*x.^2/2))/2;

egh2 = 2/sqrt(pi) * exp(-x.^2) .*(0.6405*exp(-0.6004*x) + 0.2457*exp(-2.5048*x));

%eChiani = 1/6*exp(-x.^2) + 1/2*exp(-4/3*x.^2);
eChiani = 1/12*exp(-x.^2/2) + 1/4*exp(-2/3*x.^2);

% esya3 = 1/6*exp(-x.^2) + 1/3*exp(-4*x.^2) + 1/3*exp(-4/3*x.^2);
% esya4 = 1/8*exp(-x.^2)+ 1/4*exp(-2*x.^2)+ 1/4*exp(-20/3*x.^2)+ 1/4*exp(-20/17*x.^2);

esya3 = 1/12*exp(-x.^2/2) + 1/6*exp(-2*x.^2) + 1/6*exp(-2/3*x.^2);
esya4 = 1/16*exp(-x.^2/2)+ 1/8*exp(-x.^2)+ 1/8*exp(-10/3*x.^2)+ 1/8*exp(-10/17*x.^2);


% eMbfc1 = 2* exp(-2*0.3845*x.^2 - sqrt(2)*0.7635*x - 0.6966);
eMbfc1 = exp(-0.3845*x.^2 - 0.7635*x - 0.6966);
%eMbfc1 = 2* exp(-2*0.3842*x.^2 - sqrt(2)*0.7640*x - 0.6964);
% eMbfc2 = 2* exp(-2*0.4774*x.^2 - sqrt(2)*0.4484*x - 0.9049);
eMbfc2 = exp(-0.4774*x.^2 - 0.4484*x - 0.9049);
%eMbfc2 = 2* exp(-2*0.4920*x.^2 - sqrt(2)*0.2887*x - 1.1893);
% etag3p = 2* exp(-0.000015*x.^6 - 0.00062*x.^5 + 0.00895*x.^4 -0.067133*x.^3-0.6893*x.^2-1.0948*x - 0.69874);
etag3p = exp(-0.0000019*x.^6 - 0.000109*x.^5 + 0.002238*x.^4 -0.023735*x.^3-0.344644*x.^2-0.774128*x - 0.69874);

% eAmpj = exp(-0.748*x.^2 - 1.0988*x);
eAmpj = 0.5*exp(-0.374*x.^2 - 0.777*x);

eDE = 0.9702/sqrt(pi) .* exp(-x.^2/1.0028)./(.6013*x + 0.3987.*(x.^2+1.8339)) - (0.0013 * exp(-4.7380*x.^2)+0.0054 * exp(-2.0872*x.^2));

eBeau = zeros(1,nTerms);
for n=1:2:1001
    eBeau = eBeau + (exp(-n^2*(2*pi/28)^2/2)*sin(sqrt(2)*n*2*pi/28.*x))/n;
end

eBeau = 1 - 4/pi*eBeau;

eSora = 2 - 2*2.^(-22.^(1 - 41.^(x*sqrt(2)/10)));

ePoly = zeros(1,nTerms);

for m=0:8
    for p =0:8
        ePoly = ePoly + (-1)^(m+p)*(factorial(8)/(factorial(p)*factorial(8-p)))/(factorial(m)*factorial(8-m))*(4/3)^(p/2)*(4-m)^(8-p)*x.^p.*heaviside(x - sqrt(3/4)*(4-m));
    end
end
ePoly = 2- 2*ePoly;
 
eapproxs = zeros(20,nTerms);
eapproxs(1,:) = eBorj;
eapproxs(2,:) = eKara;
eapproxs(3,:) = eIsuk;
eapproxs(4,:) = eProny2;
eapproxs(5,:) = eProny3;
eapproxs(6,:) = eOlab1;
eapproxs(7,:) = eOlab2;
eapproxs(8,:) = eOlab3;
eapproxs(9,:) = egh2;
eapproxs(10,:) = eChiani;
eapproxs(11,:) = eMbfc1;
eapproxs(12,:) = eMbfc2;
eapproxs(13,:) = etag3p;
eapproxs(14,:) = eAmpj;
eapproxs(15,:) = eDE;
eapproxs(16,:) = eBeau;
eapproxs(17,:) = eSora;
eapproxs(18,:) = ePoly;
eapproxs(19,:) = esya3;
eapproxs(20,:) = esya4;

% maximum absolute error
abserr = zeros(20,nTerms);
relerr = zeros(20,nTerms);
msse = zeros(1,20);
for j=1:20
    abserr(j,:) = (abs(exact - eapproxs(j,:)));
    relerr(j,:) = (abs(exact - eapproxs(j,:))./exact);
    msse(j) = sum((exact - eapproxs(j,:)).^2);
end

msse = msse/nTerms;

% bounds based approximation error
% set(gca,'fontsize',14);
% grid on
figure
semilogy(x, relerr(1,:),'k', 'Linewidth', 2);
hold on
semilogy(x, relerr(2,:),'b', 'Linewidth', 2);
hold on
semilogy(x, relerr(3,:),'--b', 'Linewidth', 2);
xlabel('x', 'FontSize', 14);
ylabel('Relative Error','FontSize', 14);
figure
% Exponential Function approximation error
set(gca,'fontsize',14);
%grid on;
semilogy(x, relerr(4,:),'k', 'Linewidth', 2);
hold on
semilogy(x, relerr(5,:),'--k', 'Linewidth', 2);
semilogy(x, relerr(6,:),'r', 'Linewidth', 2);
semilogy(x, relerr(7,:),'m', 'Linewidth', 2);
semilogy(x, relerr(8,:),'--r', 'Linewidth', 2);
semilogy(x, relerr(9,:),'c', 'Linewidth', 2);
semilogy(x, relerr(10,:),'y', 'Linewidth', 2);
semilogy(x, relerr(19,:),'b', 'Linewidth', 2);
semilogy(x, relerr(20,:),'--b', 'Linewidth', 2);
xlabel('x', 'FontSize', 14);
ylabel('Relative Error','FontSize', 14);
legend('E_{Prony}^2', 'E_{Prony}^3', 'E_{Olab}^1', 'E_{Olab}^2', 'E_{Olab}^3', 'E_{GH}^2', 'E_{Chiani}','E_{SYA}^3', 'E_{SYA}^4');
% %exponent polynomial based
figure
set(gca,'fontsize',14);
grid on
semilogy(x, relerr(11,:),'--k', 'Linewidth', 2);
hold on
semilogy(x, relerr(12,:),'--b', 'Linewidth', 2);
semilogy(x, relerr(13,:),'--r', 'Linewidth', 2);
semilogy(x, relerr(14,:),'--g', 'Linewidth', 2);
xlabel('x', 'FontSize', 14);
ylabel('Relative Error','FontSize', 14);
legend('E_{MBFC}^{SSE}(45)', 'E_{MBFC}^{MARE}(46)', 'E_{TAG3P}(47)','E_{AMPJ}(48)');

% legend('E_{Kara}', 'E_{Prony}^2', 'E_{Prony}^3', 'E_{Olab}^1', 'E_{Olab}^2', 'E_{Olab}^3', 'E_{GH}^2', 'E_{Chiani}','E_{SYA}^3', 'E_{SYA}^4', 'E_{MBFC}^{MARE}','E_1','E_2');
% exponent polynomial based
figure
set(gca,'fontsize',14);
grid on
semilogy(x, relerr(15,:),'k', 'Linewidth', 2);
hold on
semilogy(x, relerr(16,:),'-b', 'Linewidth', 0.5);
semilogy(x, relerr(17,:),'*r', 'Linewidth', 0.5);
semilogy(x, relerr(18,:),'g', 'Linewidth', 2);
xlabel('x', 'FontSize', 14);
ylabel('Relative Error','FontSize', 14);

legend('E_{DE}(49)', 'E_{Beau}(50)', 'E_{Sora}(51)','E_{Poly}(52)');
legend('boxoff')

figure
%exponent polynomial based
% set(gca,'fontsize',14);
% grid on
% semilogy(x, exact,'ok', 'Linewidth', 2);
% hold on
% semilogy(x, eMbfc1,'b', 'Linewidth', 2);
% hold on
% semilogy(x, eMbfc2,'r', 'Linewidth', 2);
% semilogy(x, etag3p,'c', 'Linewidth', 2);
% semilogy(x, eAmpj,'m', 'Linewidth', 2);
% xlabel('x', 'FontSize', 14);
% ylabel('erfc(x)','FontSize', 14);

figure
%exponent polynomial based
set(gca,'fontsize',14);
grid on
semilogy(x, exact,'ok', 'Linewidth', 2);
hold on
semilogy(x, eMbfc1,'b', 'Linewidth', 2);
hold on
semilogy(x, eMbfc2,'r', 'Linewidth', 2);
semilogy(x, etag3p,'c', 'Linewidth', 2);
semilogy(x, eAmpj,'m', 'Linewidth', 2);
xlabel('x', 'FontSize', 14);
ylabel('erfc(x)','FontSize', 14);


legend('Exact', 'E_{MBFC}^{SSE}', 'E_{MBFC}^{MARE}()', 'E_{TAG3P}','E_{AMPJ}');

figure
set(gca,'fontsize',18);
semilogy(x, relerr(1,:),'-k', 'Linewidth', 2);
hold on
semilogy(x, relerr(2,:),'-b', 'Linewidth', 2);
semilogy(x, relerr(10,:),'-r', 'Linewidth', 2);
semilogy(x, relerr(18,:),'-g', 'Linewidth', 2);
xlabel('x', 'FontSize', 18);
ylabel('Absolute Relative Error','FontSize', 18);

legend('E_{Broj}()', 'E_{Kara}()', 'E_{Chiani}()','E_{Poly}()');
legend('boxoff')

% set(gca,'fontsize',14);
% grid on;
% semilogy(x, relerr(4,:),'k', 'Linewidth', 2);
% hold on
% semilogy(x, relerr(5,:),'--k', 'Linewidth', 2);
% semilogy(x, relerr(6,:),'r', 'Linewidth', 2);
% semilogy(x, relerr(7,:),'--r', 'Linewidth', 2);
% semilogy(x, relerr(8,:),'c', 'Linewidth', 2);
% semilogy(x, relerr(9,:),'y', 'Linewidth', 2);
% semilogy(x, relerr(10,:),'m', 'Linewidth', 2);
% semilogy(x, relerr(2,:),'b', 'Linewidth', 2);
% xlabel('x', 'FontSize', 14);
% ylabel('Relative Error','FontSize', 14);
% legend('E_{Prony}^2', 'E_{Prony}^3', 'E_{Olab}^1', 'E_{Olab}^2', 'E_{Olab}^3', 'E_{GH}^2', 'E_{Chiani}','E_{Kara}');

% hold on
% semilogy(x, relerr(11,:),'r', 'Linewidth', 2);
%hold on
%semilogy(x, relerr(12,:),'g', 'Linewidth', 2);