% Set-up:
setup

% Parameters:
k = 2*pi;                      % wavenumber
lambda = 2*pi/k;               % wavelength
nx = 100;                      % number of points in x-dir
ny = nx;                       % number of points un y-dir
scl = 5;                       % distance of measurements
epsilon = 1e-2;                % size of random scaterrer
q = 2;                         % distance parameter for controlled source
p = 1;                         % distance parameter for random scatterer
eta = 5e-3;                    % noise level
xmin = -(scl+1)*lambda;        % domain is [xmin, xmax] x [ymin, ymax]
xmax = (scl+1)*lambda;         % domain is [xmin, xmax] x [ymin, ymax]
ymin = xmin;                   % domain is [xmin, xmax] x [ymin, ymax]
ymax = xmax;                   % domain is [xmin, xmax] x [ymin, ymax]
R = 1;                         % number of small scatterers
L = 150;                       % number of realizations
J = 120;                       % number of measurement points
betamax = 0.1;                 % parameter for beta-perturbed grid
sources = 'beta';              % 'beta', 'rand'
obstacle = 'kite';             % 'circ', 'elli', 'kite'
nobstacle = 1;                 % 1 or 2
aperture = 'full';             % 'full', 'limi'
style = 'paper';               % 'slide', 'paper'
levelset = 'n';                % 'y', 'n'

% Obstacle:
s = lambda/2;
if strcmp(obstacle, 'circ') % circle
    x = @(t) cos(t);
    y = @(t) sin(t);
    c = 2*lambda + 2*lambda*1i;
elseif strcmp(obstacle, 'elli') % ellipse
    x = @(t) 1.5*cos(t);
    y = @(t) sin(t);
    c = -2*lambda - 2*lambda*1i;
elseif strcmp(obstacle, 'kite') % kite
    x = @(t) cos(t) + 0.65*cos(2*t) - 0.65;
    y = @(t) 1.5*sin(t);
    if strcmp(aperture, 'full')
        c = 2*lambda + 2*lambda*1i; % full
    else
        c = -2*lambda + 2*lambda*1i; % limi
    end
end
x = @(t) s*x(t) + real(c);
y = @(t) s*y(t) + imag(c);
z = @(t) x(t) + 1i*y(t);
area = sum(chebfun(x, [0, 2*pi]).*diff(chebfun(y, [0, 2*pi])));

% Point source Z:
Rz = epsilon^(-q)*lambda;
tz = pi;
Z = Rz*[cos(tz), sin(tz), 0];
phi = @(X, Y) 1i/4*besselh(0, k*vecnorm(X - Y, 2, 2));

% Sensors X:
Rx = scl*lambda;
if strcmp(aperture, 'full')
    h = 2*pi/J;
    tx = h*linspace(0, J - 1, J)';
    X = Rx*[cos(tx), sin(tx), zeros(size(tx))];
elseif strcmp(aperture, 'limi')
    %dtx = 2*pi/3;
    %dtx = 3*pi/3;
    dtx = 4*pi/3;
    tx = linspace(-dtx/2, dtx/2, J/2)';
    X = Rx*[cos(tx), sin(tx), zeros(size(tx))];
    X = [X; 1.1*Rx*[cos(tx), sin(tx), zeros(size(tx))]];
end

% Compute near-field matrix (scattered fields):
tic
h = 2*pi/L;
N = zeros(J, L);
for l = 1:L
    
    % Small random scatterers Y:
    if strcmp(sources, 'beta') && R == 1
        h = 2*pi/L;
        ty = 2*pi/L*(l - 1 + betamax*rand);
    elseif strcmp(sources, 'beta') && R > 1
        error('Use ''rand'' for R > 1.')
    elseif strcmp(sources, 'rand')
        ty = 2*pi*rand(R, 1);
    end
    Ry = epsilon^(-p)*lambda;
    Y = Ry.*[cos(ty), sin(ty), zeros(size(ty))];

    % Compute near-field:
     N(:, l) = nearFieldRandScatterer(z, k, X,  @(X) phi(X, Z), Y, epsilon, nobstacle);

end
toc

% Add noise to measurements:
noise = 1 + eta*(2*rand(size(N)) - 1) + 1i*eta*(2*rand(size(N)) - 1);
Neta = N.*noise;

% Remove the average:
Neta = Neta - mean(Neta, 2);
delta = max(1e-2, norm(N - Neta));

% Construct (noisy) modified cross-correlation matrix:
Ceta = zeros(J);
psi = @(X, Y) 1i/2*besselj(0, k*vecnorm(X - Y, 2, 2));
sigma = pi^2*abs(besselh(0, 2*pi*epsilon))^2*epsilon^(-q);
for m = 1:J
    Ceta(:, m) = (2i*k)*(2*pi*Ry/(L*R))*Neta(m, :)*Neta';
    Ceta(:, m) = sigma*Ceta(:, m) - psi(X, X(m, :));
end

% SVD:
[U, S, V] = svd(Ceta);
S = diag(S);
S2 = S.^2;
amin = delta*min(S);
amax = delta*max(S);

% Morozov functional:
morozov = @(a, Y) sum((Y.*(a^2 - (delta^2*S2)))./(a + S2).^2);

% Probe the domain:
tic
[xx, yy] = meshgrid(linspace(xmin, xmax, nx), linspace(ymin, ymax, ny));
I = zeros(ny, nx);
for ix = 1:nx
    for iy = 1:ny
        
        % Construct the rhs:
        Uphi = U'*phi(X, [xx(iy, ix), yy(iy, ix), 0]);
        
        % Find the zero of the Morozov functional:
        alpha = fzero(@(a) morozov(a, abs(Uphi).^2), [amin, amax]);
        
        % Compute indicator function:
        Vg = S./(alpha + S2).*Uphi;
        I(iy, ix) = 1/norm(Vg);
        
    end  
end
toc

% Compute area:
if levelset == 'y'
    e = [];
    prmvalues = 0.35:0.01:0.85;
    for prm = prmvalues
        lvl = I > prm * max(abs(I(:)));
        P = [xx(lvl), yy(lvl)];
        [~, av] = convhull(P);
        e = [e; abs(av - area)/area];
    end
    [minimum, pos] = min(e);
    prm = prmvalues(pos);
    lvl = I > prm * max(abs(I(:)));
    P = [xx(lvl), yy(lvl)];
    [k, av] = convhull(P);
    fprintf('Area error %1.2e with level set at %1.2f*max.\n', minimum, prm)
end

% Plot matrix:
surf(imag(Ceta)), view(0, 90), shading interp
axis equal, axis([1, J, 1, J])
xticks([1 40 80 120])
yticks([1 40 80 120])
if strcmp(style, 'slide')
    xticks(''), yticks('')  
end
if strcmp(style, 'paper') && R == 1
    title('$\mathrm{Im}(\widetilde{C}_\delta)$', 'interpreter', 'Latex')
elseif strcmp(style, 'paper') && R > 1
    title('$\mathrm{Im}(\widetilde{C}^R_\delta)$', 'interpreter', 'Latex')
end
set(gcf, 'Position', [100 100 800 800])
shg, print(gcf, 'images/matrix.eps', '-depsc', '-r280')

% Plot indicator function:
if strcmp(aperture, 'full')
    I(sqrt(xx.^2 + yy.^2) > 1.0*scl*lambda) = 0;
end
figure, surf(xx, yy, I), view(2), shading interp
axis equal, axis([xmin, xmax, ymin, ymax])
xticks([-6 -2 2 6]), yticks([-6 -2 2 6])
if strcmp(style, 'slide')
    xticks(''), yticks('')  
end
clim([0, 1])
t = linspace(0, 2*pi, 100);
hold on, plot3(x(t), y(t), 1e1 + max(I(:))*ones(size(y(t))), '-k', LW, 2)
if levelset == 'y'
    plot3(P(k,1), P(k,2), 1e1 + max(I(:))*zeros(size(P(k,1))), '--k', LW, 1.5)
end
if nobstacle == 2
    hold on, plot3(-x(t), -y(t), 1e1 + max(I(:))*ones(size(y(t))), '-k', LW, 2)
end
hold on, plot3(X(1:1:end,1), X(1:1:end,2), 1e1 + ones(size(X(1:1:end,1))), 'xw', MS, 12)
if strcmp(style, 'paper')
    title('$\Vert g_s\Vert_2^{-1}$', 'interpreter', 'Latex')
end
set(gcf, 'Position', [100 100 800 800])
shg, print(gcf, 'images/lsm.eps', '-depsc', '-r280')