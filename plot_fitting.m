function plot_fitting(varargin)
    % varargin{1} : name of the first file
    % varargin{2} : name of variable in the first file
    % varargin{3} : name of the second file
    % varargin{4} : name of variable in the second file
    % varargin{5} : point coordinates to plot
    % varargin{6} : title of the plot
    % varargin{7} : size of the figure
    
    filename1 = varargin{1};
    varname1  = varargin{2};

    filename2 = varargin{3};
    varname2  = varargin{4};
    subplotsize = varargin{5};
    
    if nargin >= 6
        coordinate_points = varargin{6};
    elseif nargin == 5
        drawpoints = 'y'
    end
    if nargin >=7
        Title = varargin{7};
    end
    if nargin >=8
        FigureSize = varargin{8};
    end
    
    data1 = load(filename1);
    data2 = load(filename2);
    indx  = load(coordinate_points);
    
    signal1 = getfield(data1,varname1);
    signal2 = getfield(data2,varname2);
    x       = getfield(indx,'x');
    y       = getfield(indx,'y');
    
    fig= figure;
    for ii=1:length(x)
        xp = x(ii);
        yp = y(ii);
        subplot(subplotsize(1),subplotsize(2),ii)
        plot(signal1(:,xp,yp),'r-','linewidth',1.5), hold on
        plot(signal2(:,xp,yp),'b--','linewidth',1.0), hold on
        title([num2str(ii)])
        set(gca,'FontSize',14)
        if ii ==1
            legend('orig','recon')
        end
    end
    if nargin >=7
        sgtitle(Title)
    end
    if nargin >= 8 
        set(fig,'Units', 'normalized', 'Position', FigureSize)
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    