function h=abline(a,b,varargin)
% FORMAT h = abline(a,b,...)
% Plots y=a+b*x in dotted line
% FORMAT h = abline('h',y,...)
% Plots a horizontal line at y; y can be a vector, & then multiple lines plotted
% FORMAT h = abline('v',x,...)
% Plots a vertical line at x; x can be a vector, & then multiple lines plotted
%
% ...  Other graphics options, e.g. "'LineStyle','-'" or "'LineWidth',2" or
%      "'color',[1 0 0]",  etc
%
% Like Splus' abline.  Line is plotted and then moved behind all other
% points on the graph.
%______________________________________________________________________________
% Author: T. Nichols
% Version: http://github.com/nicholst/matlab/tree/$Format:%h$
%          $Format:%ci$

if (nargin==2) & isstr(a)
  a = lower(a);
else

  if (nargin<1)
    a = 0;
  end
  if (nargin<2)
    b = 0;
  end
  
end

XX=get(gca,'Xlim');
YY=get(gca,'Ylim');
h_exist = get(gca,'children');

g = [];
if isstr(a) & (a=='h')

  for i=1:length(b)
    g=[g;line(XX,[b(i) b(i)],'LineStyle',':','color',0.5*[1 1 1],varargin{:})];
  end

elseif isstr(a) & (a=='v')

  for i=1:length(b)
    g=[g;line([b(i) b(i)],YY,'LineStyle',':','color',0.5*[1 1 1],varargin{:})];
  end

else

  g=line(XX,a+b*XX,'LineStyle',':','color',0.5*[1 1 1],varargin{:});

end

uistack(h_exist,'top');

set(gca,'Xlim',XX);
set(gca,'Ylim',YY);

if (nargout>0)
  h=g;
end
