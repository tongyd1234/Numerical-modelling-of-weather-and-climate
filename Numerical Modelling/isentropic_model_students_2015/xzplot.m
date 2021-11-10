% -------------------------------------------------------------------
% Produces xz-plots of the isentrop model output
%
% Input:
%   simname	    Simulation name
%   time	    output frame number
%
%   Variable to be plotted:
%   'varn', 'u'         u velocity
%                       qv specific humidity
%                       qc cloud water content
%                       qr rain water content
%                       nc cloud droplet number density
%                       nr rain droplet number density
%                       lh latent heating
%
%   Domain:
%   'xlim', [100,400]   set the x-axis window
%   'zlim', [0,9]       set the z-axis window
%
%   Potential temperature contours:
%   'tci', 1            set theta contouring interval to 1K
%   'tlim', [300,350]   restrict the theta contours to [300,350]
%
%   Velocity contours:
%   'vci', 1            set velocity contouring interval to 1 m/s
%   'vlim', [0,40]      restrict the velocity contours to [0,40]
%
%   Specific humidity contours:
%   'qvci', 0.5         set qv contouring interval to 1 g/kg
%   'qvlim', [0,0.5]    restrict the qv contours to [0,0.5]
%
%   Cloud water content contours:
%   'qcci', 0.5         set qc contouring interval to 1 g/kg
%   'qclim', [0,0.5]    restrict the qc contours to [0,0.5]
%
%   Rain water content contours:
%   'qrci', 0.5         set qr contouring interval to 1 g/kg
%   'qrlim', [0,0.5]    restrict the qr contours to [0,0.5]
%
%   Cloud droplet number density contours:
%   'ncci', 10^7        set nc contouring interval to 10^7 1/kg
%   'nclim', [10^7,10^8]    restrict the nc contours to [10^7,10^8]
%
%   Cloud droplet number density contours:
%   'nrci', 10^4        set nr contouring interval to 10^4 1/kg
%   'nrlim', [10^4,10^5]    restrict the nr contours to [10^4,10^5]
%
%   Latent heating contours: 
%   'lhci', 1           set lh contouring interval to 1 K/h
%   'lhlim', [0,2]      add latent heating to [0, 2]
%
% Example:
%   xzplot('long_def',9,'varn','u','vci',4,'vlim',[0,60])
%
% -------------------------------------------------------------------
function xzplot(simname, time, varargin)

%addpath '../isentrop1/dat/'
%addpath '../isentrop2/dat/'
%addpath '../isentrop3/dat/'

% read simulation data
v = readsim(simname);
s = size(v.x);
nx = s(1);
nz = s(2);

v.plotu = 0;
v.plotqv = 0;
v.plotqc = 0;
v.plotqr = 0;
v.plotlh = 0;
v.plotnr = 0;
v.plotnc = 0;

% default values
device = '';
fig_fn = strcat('fig/', simname, '_xz');

% override default values
if nargin > 2
    n = length(varargin);
    for i=1:2:n
	switch varargin{i}
        case 'xlim'
                v.xlim = varargin{i+1};
        case 'zlim'
                v.zlim = varargin{i+1};
	    case 'tci'
                v.tci = varargin{i+1};
        case 'vci'
                v.vci = varargin{i+1}; 
        case 'qvci'
                v.qvci = varargin{i+1};
        case 'qcci'
                v.qcci = varargin{i+1};
        case 'qrci'
                v.qrci = varargin{i+1};
        case 'lhci'
                v.lhci = varargin{i+1};
        case 'ncci'
                v.ncci = varargin{i+1};
        case 'nrci'
                v.nrci = varargin{i+1};
        case 'vlim'
                v.vlim = varargin{i+1};
        case 'tlim'
                v.tlim = varargin{i+1};
        case 'qvlim'
                v.qvlim = varargin{i+1};
        case 'qclim'
                v.qclim = varargin{i+1};
        case 'qrlim'
                v.qrlim = varargin{i+1};
        case 'lhlim'
                v.lhlim = varargin{i+1};
        case 'nclim'
                v.nclim = varargin{i+1};
        case 'nrlim'
                v.nrlim = varargin{i+1};
        case 'v0'
                v.v0 = varargin{i+1};
        case 'dev'
               device = varargin{i+1};
        case 'fign'
                fig_fn = varargin{i+1};
                fig_fn = strcat('fig/', fig_fn);
        case 'varn'
                varn = varargin{i+1};
                switch varn
                    case 'u'
                        v.plotu = 1;
                    case 'qv'
                        v.plotqv = 1;
                    case 'qc'
                        v.plotqc = 1;
                    case 'qr'
                        v.plotqr = 1;
                    case 'nc'
                        v.plotnc = 1;
                    case 'nr'
                        v.plotnr = 1;
                    case 'lh'
                        v.plotlh = 1;
                end %switch
	end %switch
    end %for
end %if

% initialize figure
figure('Name', strcat(simname, ' T=', num2str(time)));
ax1=axes('XLim',v.xlim,'YLim',v.zlim);
title(['TIME = ' num2str(v.t(time)) ' seconds']);
hold on;

if v.plotqr==1
    clf
    ax1=axes('Position',[0.1 0.45 0.8 0.5],'XLim',v.xlim,'YLim',v.zlim);
    ax2=axes('Position',[0.1 0.1 0.8 0.3],'XLim',v.xlim,'YLim',v.preclim);
    % rain mixing ratio qr
    set(ax1,'NextPlot','add') 
    hold on;
end

% produce contour plot
[c,h] =contour(ax1,v.x(:,:)', v.z(:,:,time)', v.theta(:,:)', v.tlim(1):v.tci:v.tlim(2),'LineColor',[0.5 0.5 0.5]); 
set(h, 'LineWidth', 1);
plot(ax1,v.x(:,1), v.topo(:,time), 'k', 'LineWidth', 2);
if v.plotu==1
    % horizontal velocity u
    [c,h] = contour(ax1,v.x(:,:)', v.z(:,:,time)', v.u(:,:,time)', v.vlim(1):v.vci:(v.v0-v.vci), 'g');
    set(h, 'LineWidth', 2) 
    clabel(c,h);
    [c,h] = contour(ax1,v.x(:,:)', v.z(:,:,time)', v.u(:,:,time)', (v.v0+v.vci):v.vci:v.vlim(2), 'r');
    set(h, 'LineWidth', 2) 
    clabel(c,h);
    plot(ax1,v.x(:,1), v.topo(:,time), 'k', 'LineWidth', 2);
end 
if v.plotqv==1
    % water vapor qv
    [c,h] = contour(ax1,v.x(:,:)', v.z(:,:,time)', v.qv(:,:,time)', v.qvlim(1):v.qvci:v.qvlim(2), 'LineColor',[186 212 244]/256.);
    set(h, 'LineWidth', 2)
    clabel(c,h);
    plot(ax1,v.x(:,1), v.topo(:,time), 'k', 'LineWidth', 2);
end
if v.plotqc==1
    % cloud water qc
    [c,h] = contour(ax1,v.x(:,:)', v.z(:,:,time)', v.qc(:,:,time)', v.qclim(1):v.qcci:v.qclim(2), 'LineColor',[0 191 191]/256.);
    set(h, 'LineWidth', 2);
    clabel(c,h);
    plot(ax1,v.x(:,1), v.topo(:,time), 'k', 'LineWidth', 2);
end
if v.plotqr==1
    % rain water qr
    [c,h] = contour(ax1,v.x(:,:)', v.z(:,:,time)', v.qr(:,:,time)', v.qrlim(1):v.qrci:v.qrlim(2), 'LineColor',[0 230 150]/256.);
    clabel(c,h);
    set(h, 'LineWidth', 2)
    clabel(c,h);
    plot(ax1,v.x(:,1), v.topo(:,time), 'k', 'LineWidth', 2);
    plot(ax2,v.x(:,1), v.totprec(:,time),'Color',[0 230 150]/256., 'LineWidth', 2); 
    set(ax2,'XLim',v.xlim,'YLim',v.totpreclim)
    set(ax1,'XLim',v.xlim,'YLim',v.zlim)
end
if v.plotnc==1
    % cloud water droplet density
    [c,h] = contour(ax1,v.x(:,:)', v.z(:,:,time)', v.nc(:,:,time)', v.nclim(1):v.ncci:v.nclim(2), 'LineColor',[255 200 0]/256.);
    set(h, 'LineWidth', 2)
    set(h, 'LineStyle', '--')
    clabel(c,h);
    plot(ax1,v.x(:,1), v.topo(:,time), 'k', 'LineWidth', 2);
end
if v.plotnr==1
    % rain water droplet density
    [c,h] = contour(ax1,v.x(:,:)', v.z(:,:,time)', v.nr(:,:,time)', v.nrlim(1):v.nrci:v.nrlim(2), 'LineColor',[255 69 0]/256.);
    set(h, 'LineWidth', 2)
    set(h, 'LineStyle', '--')
    clabel(c,h);
    plot(ax1,v.x(:,1), v.topo(:,time), 'k', 'LineWidth', 2);
end %if
if v.plotlh == 1
    [c2,h2] = contour(ax1,v.x(:,:)', v.z(:,:,time)', v.lheat(:,:,time)'.*3600., v.lhlim(1):v.lhci:-v.lhci, 'b');
    [c3,h3] = contour(ax1,v.x(:,:)', v.z(:,:,time)', v.lheat(:,:,time)'.*3600., v.lhci:v.lhci:v.lhlim(2), 'y');
    set(h2, 'LineWidth', 2);
    set(h3, 'LineWidth', 2);
    clabel(c2,h2);
    clabel(c3,h3);
end %if
hold off;

xlabel('x [km]')
if (v.plotqr==1)
  ylabel('Acum. Rain [mm]')
else
  ylabel('Height [km]')
end

