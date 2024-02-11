%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sweep over values for kappa of BRW turning angle distribution, and
% weight of CRW component.
% CRW kappa is fixed, determined from experiment.
% Step size and arena size are currently fixed
% Output:
%   R_vec_length_mean: array of mean resultant lengths
% Input:
%   BRW_array: vector of BRW kappa values
%   w_array: vector of weight values
%   CRW_array: vector of CRW kappa values
%   step_size_array: vector of step sizes
%   N_iter: number of random walks to simulate

% Requires installation of Circular Statistics toolbox: https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics
function [R_vec_length_mean] = insect_nav(BRW_array,w_array,CRW_array,step_size_array,N_iter)
rng(1);

arena_size = 20; % arena size

% initial conditions
x_init = 0;
y_init = 0;
theta_init = 0;

Tortuosity_mean = zeros(length(BRW_array), length(w_array), length(step_size_array));
Open_angle_mean = Tortuosity_mean;
R_vec_length_mean = Tortuosity_mean;

for BRW_i = 1:length(BRW_array)
    tic;
    for w_i = 1:length(w_array)
        for step_size_i = 1:length(step_size_array)
            
            kappa_CRW = CRW_array(1); % fixed CRW value
            kappa_BRW = BRW_array(BRW_i);
            w = w_array(w_i);
            step_size = step_size_array(step_size_i);
            N_steps = 1000*floor(arena_size/step_size);
            max_radius = arena_size;
            
            [xx, yy] = BCRW(N_steps, N_iter, step_size, w, kappa_BRW, kappa_CRW, ...
                x_init, y_init, theta_init, max_radius);
            
            x_last = zeros(length(xx),1);
            y_last = zeros(length(yy),1);
            num_last_steps = zeros(length(xx),1);
            
            for i = 1:length(xx)
                x_temp = xx{i}(end);
                y_temp = yy{i}(end);
                dist_temp = sqrt(x_temp^2+y_temp^2);
                x_last(i) = max_radius*x_temp/dist_temp;
                y_last(i) = max_radius*y_temp/dist_temp;
                num_last_steps(i) = length(xx{i});
            end
            
            walk_extent = sqrt(x_last.^2 + y_last.^2); % extent of each walk
            Tortuosity_mean(BRW_i, w_i, step_size_i) = mean(walk_extent ./ (num_last_steps*step_size));
            
            % get von mises distribution of final walk angles, then compute fwhm
            walk_angles = atan2(y_last, x_last);
            walk_angles(isnan(walk_angles)) = [];
            kappa = circ_kappa(walk_angles');
            [p, alpha] = circ_vmpdf(0:0.1:(2*pi), pi, kappa);
            fwhmx = fwhm(rad2deg(alpha),p);
            Open_angle_mean(BRW_i, w_i, step_size_i) = fwhmx;
            R_vec_length_mean(BRW_i, w_i, step_size_i) = circ_r(walk_angles);
            
        end
    end
    toc;
end



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xx, yy] = BCRW(N_steps, N_iter, v, w, kappa_BRW, kappa_CRW, x_init, y_init, theta_init, max_rad)

% Computes biased correlated random walk

% Inputs:
%   N_steps:        Number of steps of each walk
%   N_iter:         Number of iterations to simulate
%   v:              Velocity of walker
%   w:              Weight value on CRW component
%   kappa_BRW:      Kappa value of BRW von Mises distribution
%   kappa_CRW:      Kappa value of CRW von Mises distribution
%   x_init:         Initial x position
%   y_init:         Initial y position
%   theta_init:     Initial angle theta, and BRW biasing angle
%   max_rad:        Size of arena

% Outputs
%   xx:             x trajectory for each walk
%   yy:             y trajectory for each walk

% generate von mises PDFs for BRW and CRW turning angle distributions
[p_BRW, alpha_BRW] = circ_vmpdf(-pi:0.001:pi, 0, kappa_BRW);
[p_CRW, alpha_CRW] = circ_vmpdf(-pi:0.001:pi, 0, kappa_CRW);

xx = cell(N_iter,1);
yy = cell(N_iter,1);

x = zeros(N_steps, N_iter);
y = zeros(N_steps, N_iter);
theta = zeros(N_steps, N_iter);

x(1,:) = x_init;
y(1,:) = y_init;
theta(1,:) = theta_init;

for iter_i = 1:N_iter
    
    % randomly draw from BRW and CRW angle distributions to get arrays of
    % turning angles at each time step
    angles_BRW = randpdf(p_BRW, alpha_BRW,[ N_steps, 1]);
    angles_CRW = randpdf(p_CRW, alpha_CRW,[ N_steps, 1]);
    
    for step_i = 2:N_steps
        
        % compute next step
        dx = v*(w*cos(theta_init + angles_BRW(step_i)) + (1-w)*cos(theta(step_i-1,iter_i)+angles_CRW(step_i)));
        dy = v*(w*sin(theta_init + angles_BRW(step_i)) + (1-w)*sin(theta(step_i-1,iter_i)+angles_CRW(step_i)));
        % normalize steps
        step_length = sqrt(dx^2+dy^2);
        dx = dx/step_length;
        dy = dy/step_length;
        x(step_i, iter_i) = x(step_i-1, iter_i) + dx;
        y(step_i, iter_i) = y(step_i-1, iter_i) + dy;
        
        % angle of current step
        theta(step_i, iter_i) = atan2(dy, dx);
        
        
        % angle of current step
%         theta(step_i, iter_i) = atan2(dy, dx);
        
    end
    
    % find when edge was reached
    dist = sqrt(x(:,iter_i).^2 + y(:,iter_i).^2); % distance to center at each timestep
    edge_reached = find(dist >= max_rad, 1, 'first');
    if ~isempty(edge_reached)
        xx{iter_i} = x(1:edge_reached, iter_i);
        yy{iter_i} = y(1:edge_reached, iter_i);
    else
        error('did not reach the edge')
    end
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x=randpdf(p,px,dim)
% RANDPDF
%   Random numbers from a user defined distribution
%
% SYNTAX:
%   x = randpdf(p, px, dim)
%       randpdf(p, px, dim)
%
% INPUT:
%   p   - probability density,
%   px  - values for probability density,
%   dim - dimension for the output matrix.
%
% OUTPUT:
%   x   - random numbers. Run function without output for some plots.
%
% DESCRIPTION:
%   x = randpdf(p, px, dim) returns the matrix of random numbers from
%   probability density distribution defined in p and px. p are the density
%   (the y axis) and px are the value (the x axis) of the pdf. p and px
%   must be of the same length.
%   dim define the output matrix dimensions, for example dim=[100 3] define
%   the 100x3 two dimensional matrix with 300 random numbers.
%
%   REMEMBER: This is not a realy random number generator but only
%   some kind of transformation of uniformly distributed pseudorandom
%   numbers to desired pdf!
%
% EXAMPLE 1:
%   Generation of normal distributed random numbers. This is not typical
%   normal distribution because is limited from the left and right side,
%   i.e. 0 < px < 80 .
%
%   px=0:80;
%   p=1./(10*sqrt(2*pi))*exp((-(px-40).^2)./(2*10^2));
%   randpdf(p,px,[10000,1])
%
%
% EXAMPLE 2:
%   Generation using user defined pdf.
%
%   px=[1 2 3 4 5 6 7 8 9];
%   p= [0 1 3 0 0 4 5 4 0];
%   randpdf(p,px,[50000,1])
% By Adam Nies?ony, Opole University of Technology, Poland
% check the number of input
error(nargchk(3, 3, nargin))
% vectorization and normalization of the input pdf
px=px(:);
p=p(:)./trapz(px,p(:));
% interpolation of the input pdf for better integration
% in my opinion 10000 point is sufficient...
pxi=[linspace(min(px),max(px),500000)]';
pi=interp1(px,p,pxi,'linear');
% computing the cumulative distribution function for input pdf
cdfp = cumtrapz(pxi,pi);
% finding the parts of cdf parallel to the X axis
ind=[true; not(diff(cdfp)==0)];
% and cut out the parts
cdfp=cdfp(ind);
pi=pi(ind);
pxi=pxi(ind);
% generating the uniform distributed random numbers
uniformDistNum=rand(dim);
% and distributing the numbers using cdf from input pdf
userDistNum=interp1(cdfp,pxi,uniformDistNum(:)','linear');
% making graphs if no output exists
if nargout==0
    subplot(3,4,[1 2 5 6])
    [n,xout]=hist(userDistNum,50);
    n=n./sum(n)./(xout(2)-xout(1));
    bar(xout,n)
    hold on
    plot(pxi, pi./trapz(pxi,pi),'r')
    hold off
    legend('pdf from generated numbers','input pdf')
    subplot(3,4,[3 4 7 8])
    plot(pxi, cdfp,'g')
    ylim([0 1])
    legend('cdf from input pdf')
    subplot(3,4,[9:12])
    plot(userDistNum)
    legend('generated numbers')
else
    x=reshape(userDistNum,dim);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p alpha] = circ_vmpdf(alpha, thetahat, kappa)

% [p alpha] = circ_vmpdf(alpha, w, p)
%   Computes the circular von Mises pdf with preferred direction thetahat
%   and concentration kappa at each of the angles in alpha
%
%   The vmpdf is given by f(phi) =
%   (1/(2pi*I0(kappa))*exp(kappa*cos(phi-thetahat)
%
%   Input:
%     alpha     angles to evaluate pdf at, if empty alphas are chosen to
%               100 uniformly spaced points around the circle
%     [thetahat preferred direction, default is 0]
%     [kappa    concentration parameter, default is 1]
%
%   Output:
%     p         von Mises pdf evaluated at alpha
%     alpha     angles at which pdf was evaluated
%
%
%   References:
%     Statistical analysis of circular data, Fisher
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens and Marc J. Velasco, 2009
% velasco@ccs.fau.edu

% if no angles are supplied, 100 evenly spaced points around the circle are
% chosen
if nargin < 1 || isempty(alpha)
    alpha = linspace(0, 2*pi, 101)';
    alpha = alpha(1:end-1);
end
if nargin < 3
    kappa = 1;
end
if nargin < 2
    thetahat = 0;
end

alpha = alpha(:);

% evaluate pdf
if kappa ==0
    p = unifpdf(alpha,-2*pi,2*pi);
elseif kappa<=700
    C = 1/(2*pi*besseli(0,kappa));
    p = C * exp(kappa*cos(alpha-thetahat));
else
    std = sqrt(1./kappa);
    p = normpdf(alpha,thetahat,std);
end



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fwhmx = fwhm(x,data)

% Find the half max value.
halfMax = (min(data) + max(data)) / 2;
% Find where the data first drops below half the max.
index1 = find(data <= halfMax, 1, 'first');
% Find where the data last rises above half the max.
index2 = find(data >= halfMax, 1, 'last');
%fwhm = index2-index1 + 1; % FWHM in indexes.
% OR, if you have an x vector
fwhmx = x(index2) - x(index1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s s0] = circ_std(alpha, w, d, dim)
% s = circ_std(alpha, w, d, dim)
%   Computes circular standard deviation for circular data
%   (equ. 26.20, Zar).
%
%   Input:
%     alpha	sample of angles in radians
%     [w		weightings in case of binned angle data]
%     [d    spacing of bin centers for binned data, if supplied
%           correction factor is used to correct for bias in
%           estimation of r]
%     [dim  compute along this dimension, default is 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_std(alpha, [], [], dim)
%
%   Output:
%     s     angular deviation
%     s0    circular standard deviation
%
% PHB 6/7/2008
%
% References:
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 4
    dim = 1;
end

if nargin < 3 || isempty(d)
    % per default do not apply correct for binned data
    d = 0;
end

if nargin < 2 || isempty(w)
    % if no specific weighting has been specified
    % assume no binning has taken place
    w = ones(size(alpha));
else
    if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1)
        error('Input dimensions do not match');
    end
end

% compute mean resultant vector length
r = circ_r(alpha,w,d,dim);

s = sqrt(2*(1-r));      % 26.20
s0 = sqrt(-2*log(r));    % 26.21


end