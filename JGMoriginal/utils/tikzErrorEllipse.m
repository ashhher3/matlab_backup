function tikzErrorEllipse(eStats,conf,clrNames,tikzfilename)
% tikzErrorEllipse
% 
% USAGE:
%   tikzErrorEllipse(eStats,conf,clrNames,tikzfilename)
%
% It's much more efficient to have pgf write the ellipses than to have it
% plot a bunch of points that Matlab provided that "happen" to lie on an
% ellipse.

%-------------------------------------------------------------------------%
% Created: 01/21/15
%   by JGM
%-------------------------------------------------------------------------%


%%% TO DO:
% (1) Make the axis labels into arguments!!
% (2) Make the placement of the coordinate frame more rational/flexible
% (3) Currently, you just scale the ellipses up to 100cm, which are then
% scaled back down by the tiny covariance matrix.  It would be much better
% to pick a reasonable size for the picture (10 cm??), then work backwards
% to see what ellipse size puts you there, given the size (e.g., det) of
% the covariance matrix.  This may be problematic, b/c tikz doesn't like
% ellipses that are even nominally 1000 cm (e.g.); so in fact it may be
% best in the end just to scale up the covariance matrix ("cm=[...")
% themselves, just remembering to do it to all of them......
%%%%%%%%%%%%%%%%%%%%%%



% directory
[blank, name] = system('hostname');
switch strtrim(name)
    case 'kobayashi-maru'
        yrtikzdir = 'C:\Documents and Settings\makin\My Documents\#texs\tikzpics\';
    case {'CUPCAKE','Themistocles'}
        yrtikzdir = 'C:\Users\makin\Documents\#texs\tikzpics\';
    case {'mushroom','keck-phaser1','pepperoni','zamfir'}
        yrtikzdir = 'C:\Users\makin\Documents\';
    case 'domestica'
        yrtikzdir = '~/tikzpics/';
  otherwise
        error('unknown host -- jgm');
end


% you will need the carriagereturn 
load ../toys/filez/carriage.mat
outfile = [yrtikzdir,tikzfilename,'.tex'];



% text to be written
outtxt = ['\begin{tikzpicture}',carriagereturn];

% create the coordinate frame
outtxt = [outtxt,'\newcommand{\coordinateFrame}[2]{',carriagereturn];
outtxt = [outtxt,sprintf('\t'),'\draw [very thick,black] '];
outtxt = [outtxt,'(#1,#2) -- (#1,#2+2cm);',carriagereturn];
outtxt = [outtxt,sprintf('\t'),'\draw [very thick,black] '];
outtxt = [outtxt,'(#1,#2) -- ',carriagereturn];
outtxt = [outtxt,sprintf('\t\t'),'node[label={[label distance=-40, '];
outtxt = [outtxt,'color=black]{-4}:{0.02 rad}}]{}',carriagereturn];
outtxt = [outtxt,sprintf('\t\t'),'(#1+2cm,#2);',carriagereturn];
outtxt = [outtxt,sprintf('\t'),'}',carriagereturn];
%%% 2 = 0.02 b/c you set the ellipses at (100cm and 100cm)

% now do for each ellipse...
for iStat = 1:size(eStats.Xpct,2)
    
    % get transforming matrix via covariance (see error_ellipse.m)
    confScale = qchisq(conf,2); % r is the number of dimensions (d.o.f.)
    SigmaOneHalf = chol(confScale*eStats.Cvrn(:,:,iStat));
    
    % make the opacity set-able from the outside w/a command
    outtxt = [outtxt,...
        '\providecommand{\',clrNames{iStat},'opacity}{1}',carriagereturn];
    
    % write into tikz's matrix-scaling thing
    outtxt = [outtxt,'\begin{scope}[cm={',...
        num2str(SigmaOneHalf(1,1)),',',num2str(SigmaOneHalf(1,2)),',',...
        num2str(SigmaOneHalf(2,1)),',',num2str(SigmaOneHalf(2,2)),',(',...
        num2str(eStats.Xpct(1,iStat)),',',num2str(eStats.Xpct(2,iStat)),')}]',...
        carriagereturn];
    % outtxt = [outtxt,'\draw [',eStats.clr,',line width=2.0,opacity=\',...
    outtxt = [outtxt,sprintf('\t'),'\draw [',clrNames{iStat},...
        ',line width=2.0,opacity=\',clrNames{iStat},...
        'opacity] (0,0) {} ellipse (100 and 100);',carriagereturn];
    outtxt = [outtxt,'\end{scope}',carriagereturn];
end

% now place the coordinate frame
outtxt = [outtxt,'\coordinateFrame{-5.5cm}{-6.8cm}',carriagereturn];
outtxt = [outtxt,'\node[text=black,anchor=west,rotate=90] (yaxislabel) '];
outtxt = [outtxt,'at (-5.5,-4.7) {$\prop_2$ (elbow)};',carriagereturn];
outtxt = [outtxt,'\node[text=black,anchor=west] (yaxislabel) '];
outtxt = [outtxt,'at (-3.5,-6.8) {$\prop_1$ (shoulder)};',carriagereturn];
%%% you should really compute all the numbers in here!!

% close off the picture
outtxt = [outtxt,'\end{tikzpicture}%'];

% write to file
fid = fopen(outfile,'wt+');
fwrite(fid,outtxt);
fclose(fid);

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function x = qchisq(P,n)
% QCHISQ(P,N) - quantile of the chi-square distribution.

%-------------------------------------------------------------------------%
% Cribbed: 01/21/15
%   -from error_ellipse.m
%   -by JGM
%-------------------------------------------------------------------------%

if nargin<2
    n=1;
end

s0 = P==0;
s1 = P==1;
s = P>0 & P<1;
x = 0.5*ones(size(P));
x(s0) = -inf;
x(s1) = inf;
x(~(s0|s1|s))=nan;

for ii=1:14
    dx = -(pchisq(x(s),n)-P(s))./dchisq(x(s),n);
    x(s) = x(s)+dx;
    if all(abs(dx) < 1e-6)
        break;
    end
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function F=pchisq(x,n)
% PCHISQ(X,N) - Probability function of the chi-square distribution.

%-------------------------------------------------------------------------%
% Cribbed: 01/21/15
%   -from error_ellipse.m
%   -by JGM
%-------------------------------------------------------------------------%

if nargin<2
    n=1;
end
F=zeros(size(x));

if rem(n,2) == 0
    s = x>0;
    k = 0;
    for jj = 0:n/2-1;
        k = k + (x(s)/2).^jj/factorial(jj);
    end
    F(s) = 1-exp(-x(s)/2).*k;
else
    for ii=1:numel(x)
        if x(ii) > 0
            F(ii) = quadl(@dchisq,0,x(ii),1e-6,0,n);
        else
            F(ii) = 0;
        end
    end
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function f=dchisq(x,n)
% DCHISQ(X,N) - Density function of the chi-square distribution.
%-------------------------------------------------------------------------%
% Cribbed: 01/21/15
%   -from error_ellipse.m
%   -by JGM
%-------------------------------------------------------------------------%

if nargin<2
    n=1;
end
f=zeros(size(x));
s = x>=0;
f(s) = x(s).^(n/2-1).*exp(-x(s)/2)./(2^(n/2)*gamma(n/2));

end