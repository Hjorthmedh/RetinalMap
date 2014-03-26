clear all, close all


WT = load('Stafford2009-CI-WT.txt');
KO = load('Stafford2009-CI-B2KO.txt');

WTd = WT(:,1)*1e-6;
WTci = WT(:,2);

KOd = KO(:,1)*1e-6;
KOci = KO(:,2);

p = plot(WTd*1e6,WTci,'ko-', ...
         KOd*1e6,KOci,'ro-');

hold on
% The fit from curvefit
WTa = 22.64;
WTb = 1/4302;
WTciFit = WTa*exp(-WTd/WTb) + 1;

p(3) = plot(WTd*1e6,WTciFit,'k*--');

KOa = 9.275;
KOb = 1/1260;
KOciFit = KOa*exp(-KOd/KOb) + 1;

p(4) = plot(KOd*1e6, KOciFit,'r*--');

xlabel('Distance (\mum)')
ylabel('Correlation index')

legend(p,'WT data','B2KO data', 'WT fit', 'B2KO fit')


% For the curve fitting, subtract 1
WTciMinus1 = WTci-1;
KOciMinus1 = KOci-1;


% Looking at Dans data, Sterrat mail dec 2, 2011
%
% P2: 41*2+180 = 262 degrees, radius of eyeball 2.146/2 mm
% P6: 37*2+180 = 254 degress, radius of eyeball 2.450/2 mm
%
% P2 dist side to side = 4.9mm, P6 dist = 5.4 mm
%
% Lets use 5 mm

% Comparison with Koulakov
WTdN = WTd / 5e-3;
KOdN = KOd / 5e-3;

% Here 0 means no correlation, as opposed to 1 for CI
WTciN = (WTciFit-1) / max(WTa-1);
KOciN = (KOciFit-1) / max(WTa-1);

KoulakovD = linspace(0,1);
KoulakovCI = exp(-KoulakovD/0.11);

figure
p2 = plot(WTdN,WTciN,'ko-', ...
          KOdN,KOciN,'ro-', ...
          KoulakovD, KoulakovCI, 'b--');
xlabel('Distance (normalised)')
ylabel('Correlation')

legend(p2,'WT','\beta2KO','Koulakov')
title('Data from Stafford et al 2009')

saveas(gcf,'StaffordFit.pdf','pdf')

d = linspace(0,1,100);

% Just double checking that the curves are right

figure

      WTa = 1; % 22.64/22.64
      WTb = 1/(4302*5e-3); % Assums retina is 5mm (flattened out)
      CWT = WTa*exp(-d/WTb);

      KOa = 9.275/22.64;
      KOb = 1/(1260*5e-3);
      CKO = KOa*exp(-d/KOb);
      hold on
      p3 = plot(d,CWT,'k-',d,CKO,'r-', ...
                KoulakovD,KoulakovCI, 'b--', ...
                'linewidth',2);
      
      legend(p3,'WT','\beta2KO','Triplett (2011)')
      xlabel('Intercell distance (normalised)','fontsize',24)
      ylabel('Correlation','fontsize',24)
      set(gca,'fontsize',20)
      
      saveas(gcf,'StaffordFit-rescalingCI.pdf','pdf')