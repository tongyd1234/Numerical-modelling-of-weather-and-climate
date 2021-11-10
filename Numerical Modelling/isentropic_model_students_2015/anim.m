% -------------------------------------------------------------------
% Animate the output of the isentrop model 
%
% Input:
%   simname         Simulation name
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
%
% Example:
%   anim('long_def','varn','u','vci',4,'vlim',[0,60])
%
% -------------------------------------------------------------------

function anim(simname, varargin)

v = readsim(simname);

% default values
device = '';

v.plotu = 0;
v.plotqv = 0;
v.plotqc = 0;
v.plotqr = 0;
v.plotlh = 0;
v.plotnr = 0;
v.plotnc = 0;

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

visual('Initialize', v);
