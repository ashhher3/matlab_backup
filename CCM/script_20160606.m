%% look for difference in convergence between hi & low coh for EC108
%% 

load('/newhome/kderosier/Documents/MATLAB/CCM/data/EC100_EC108_HiVsLoAHBCoh_E10_tau10_maxL400_strictThresh.mat')

%% fit each trace & find where crosses convergence threshold

thresh = 0.98;

hi.EC108.conv_AxmapH = threshCross(hi.EC108.AxmapH,hi.EC108.Lvals,thresh);
hi.EC108.conv_HxmapA = threshCross(hi.EC108.HxmapA,hi.EC108.Lvals,thresh);
lo.EC108.conv_AxmapH = threshCross(lo.EC108.AxmapH,lo.EC108.Lvals,thresh);
lo.EC108.conv_HxmapA = threshCross(lo.EC108.HxmapA,lo.EC108.Lvals,thresh);

%% set the NaN values to the max value

hi.EC108.conv_AxmapH(isnan(hi.EC108.conv_AxmapH)) = maxL;

hi.EC108.conv_HxmapA(isnan(hi.EC108.conv_HxmapA)) = maxL;

lo.EC108.conv_AxmapH(isnan(lo.EC108.conv_AxmapH)) = maxL;

lo.EC108.conv_HxmapA(isnan(lo.EC108.conv_HxmapA)) = maxL;

%% plot AxmapH vs HxmapA

[f1, h_AH, h_HA] = plot2Hists(hi.EC108.conv_AxmapH,hi.EC108.conv_HxmapA,...
    20, 'AxmapH','HxmapA','Length of data used for CCM',...
    'Probability of convergence', ...
    {'Convergence rates for AxmapH vs. HxmapA','(EC108, high Am-Hp beta coherence)'},...
    'EC108_convergenceRate_AxmapHvsHxmapA_hiCoh' );

hi.EC108.hist_AxmapH = h_AH.Values;
hi.EC108.hist_HxmapA = h_HA.Values;

[f2, h2_AH, h2_HA] = plot2Hists(lo.EC108.conv_AxmapH,lo.EC108.conv_HxmapA,...
    20, 'AxmapH','HxmapA','Length of data used for CCM',...
    'Probability of convergence', ...
    {'Convergence rates for AxmapH vs. HxmapA','(EC108, low Am-Hp beta coherence)'},...
    'EC108_convergenceRate_AxmapHvsHxmapA_loCoh' );

lo.EC108.hist_AxmapH = h2_AH.Values;
lo.EC108.hist_HxmapA = h2_HA.Values;

%% plot hi vs lo

[f3, h_hi, h_lo] = plot2Hists(hi.EC108.conv_AxmapH,lo.EC108.conv_AxmapH,...
    20, 'Hi Bcoh','Lo Bcoh','Length of data used for CCM',...
    'Probability of convergence', ...
    {'Convergence rates for high vs. low beta coherence','(EC108, AxmapH (Hp driving Am))'},...
    'EC108_convergenceRate_HiVsLo_AxmapH' );

[f4, h4_hi, h4_lo] = plot2Hists(hi.EC108.conv_HxmapA,lo.EC108.conv_HxmapA,...
    20, 'Hi Bcoh','Lo Bcoh','Length of data used for CCM',...
    'Probability of convergence', ...
    {'Convergence rates for high vs. low beta coherence','(EC108, HxmapA (Am driving Hp))'},...
    'EC108_convergenceRate_HiVsLo_HxmapA' );

%% compute KL div (both directions for each pair than average)

% one of them is too short, fix it
lo.EC108.hist_HxmapA(9:20)=0;

% AxmapH vs HxmapA when b coh is hi
d1 = KLDiv(hi.EC108.hist_AxmapH,hi.EC108.hist_HxmapA);
d2 = KLDiv(hi.EC108.hist_HxmapA,hi.EC108.hist_AxmapH);

dist.AvsH_hi = d1+d2/2;

% AxmapH vs HxmapA when b coh is lo
d1 = KLDiv(lo.EC108.hist_AxmapH,lo.EC108.hist_HxmapA);
d2 = KLDiv(lo.EC108.hist_HxmapA,lo.EC108.hist_AxmapH);

dist.AvsH_lo = d1+d2/2;

% Hi vs Lo for AxmapH
d1 = KLDiv(hi.EC108.hist_AxmapH,lo.EC108.hist_AxmapH);
d2 = KLDiv(lo.EC108.hist_AxmapH,hi.EC108.hist_AxmapH);

dist.HvsL_AxmapH = d1+d2/2;

% Hi vs Lo for HxmapA
d1 = KLDiv(hi.EC108.hist_HxmapA,lo.EC108.hist_HxmapA);
d2 = KLDiv(lo.EC108.hist_HxmapA,hi.EC108.hist_HxmapA);

dist.HvsL_HxmapA = d1+d2/2;
