% Set-up:
setup

% Parameters:
k = 2*pi;                      % wavenumber
lambda = 2*pi/k;               % wavelength
eta = 0;                       % noise level
nx = 100;                      % number of points in x-dir
ny = nx;                       % number of points un y-dir
scl = 5;                       % distance of sources & measurements
sclxcor = 50;                  % distance of random sources
betamax = 0.1;                 % distribution of ramdom sources
xmin = -floor((scl+1)*lambda); % domain is [xmin, xmax] x [ymin, ymax]
xmax = floor((scl+1)*lambda);  % domain is [xmin, xmax] x [ymin, ymax]
ymin = xmin;                   % domain is [xmin, xmax] x [ymin, ymax]
ymax = xmax;                   % domain is [xmin, xmax] x [ymin, ymax]
L = 80;                        % number of random sources
J = L;                         % number of measurement points X
M = L;                         % number of point sources Y
imaginary = 'n';               % 'y', 'n'
xcor = 'n';                    % 'y', 'n'
xcorsetup = 'A';               % 'A', 'B'
obstacle = 'circ';             % 'circ', 'elli', 'kite'
aperture = 'full';             % 'full', 'limi'

% Obstacle:
s = lambda/2;
if strcmp(obstacle, 'circ') % circle
    x = @(t) cos(t);
    y = @(t) sin(t);
    c = 0;
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

% Point sources:
if strcmp(xcor, 'y')
    Ry = sclxcor*lambda;
else
    Ry = scl*lambda;
end
h = 2*pi/M;
theta = h*linspace(0, M - 1, M)';
if strcmp(xcor, 'y') && strcmp(xcorsetup, 'A')
    beta = betamax*rand(M, 1);
    theta = theta + h*beta;
end
if strcmp(aperture, 'limi') && ~strcmp(xcor, 'y')
    dtheta = 3*pi/4;
    theta = linspace(-dtheta/2, dtheta/2, M)';
end
Y = Ry*[cos(theta), sin(theta), zeros(size(theta))];

% Measurement points:
Rx = scl*lambda;
if strcmp(aperture, 'full')
    h = 2*pi/J;
    theta = h*linspace(0, J - 1, J)';
elseif strcmp(aperture, 'limi')
    dtheta = 3*pi/4;
    theta = linspace(-dtheta/2, dtheta/2, J)';
end
X = Rx*[cos(theta), sin(theta), zeros(size(theta))];

% Compute near-field matrix (scattered fields):
tic
N = zeros(J, M);
phi = @(X, Y) 1i/4*besselh(0, k*vecnorm(X - Y, 2, 2));
for m = 1:M
    ui = @(X) phi(X, Y(m, :));
    N(:, m) = computeNearField(z, k, X, ui);
end
toc

% Compute the total fields:
if strcmp(xcor, 'y')
    for m = 1:M
        N(:, m) = N(:, m) + phi(X, Y(m, :));
    end
    if strcmp(xcorsetup, 'B')
        Ntmp = zeros(J, M);
        for m = 1:M
            sigma = sqrt(pi*Ry/L);
            n = sigma*(randn(1, L) + 1i*randn(1, L));
            Ntmp(:, m) = Ntmp(:, m) + sum(n.*N, 2);
        end
        N = Ntmp;
    end
end

% Add noise:
noise = 1 + eta*(2*rand(size(N)) - 1) + 1i*eta*(2*rand(size(N)) - 1);
Neta = N.*noise;
delta = max(1e-2, norm(N - Neta));

% Construct matrix for LSM:
if strcmp(imaginary, 'y') % setup with imaginary near-field matrix I
    Ieta = Neta - conj(Neta);
    Aeta = Ieta;
elseif strcmp(xcor, 'y') % setup with cross-correlation matrix C
    Ceta = zeros(J, M);
    psi = @(X, Y) 1i/2*besselj(0, k*vecnorm(X - Y, 2, 2));
    for m = 1:M
        if strcmp(xcorsetup, 'A')
            Ceta(:, m) = (2*pi*Ry/L)*(2i*k)*Neta(m, :)*Neta';
        elseif strcmp(xcorsetup, 'B')
            Ceta(:, m) = (1/M)*(2i*k)*Neta(m, :)*Neta';
        end
        Ceta(:, m) = Ceta(:, m) - psi(X, X(m, :));
    end
    Aeta = Ceta;
else % setup wtih near-field matrix N
    Aeta = Neta;
end

% SVD:
[U, S, V] = svd(Aeta);
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

% Plot matrix:
if strcmp(imaginary, 'y') || strcmp(xcor, 'y')
    surf(imag(Aeta)), view(0, -90), shading interp
else
    surf(abs(Aeta)), view(0, -90), shading interp
end
axis equal, axis([1, M, 1, J])
if strcmp(xcorsetup, 'B')
    xticks([1 50 100 150 200]), yticks([1 50 100 150 200])
else
    xticks([1 20 40 60 80]), yticks([1 20 40 60 80])
end
xticks(''), yticks('')
shg, print(gcf, 'images/matrix.eps', '-depsc', '-r280')

% Plot indicator function:
if strcmp(aperture, 'full')
    I(sqrt(xx.^2 + yy.^2) > 1.0*scl*lambda) = 0;
end
figure, surf(xx, yy, I), view(2), shading interp
axis equal, axis([xmin, xmax, ymin, ymax])
t = linspace(0, 2*pi, 100);
hold on, plot3(x(t), y(t), 1e1 + max(I(:))*ones(size(y(t))), '-k', LW, 1)
hold on, plot3(Y(2:2:end,1), Y(2:2:end,2), 1e1 + ones(size(Y(2:2:end,1))), '.w', MS, 30)
if strcmp(xcor, 'y') && strcmp(xcorsetup, 'B')
    hold on
    plot3(X(1:1:end,1), X(1:1:end,2), 1e1 + ones(size(X(1:1:end,1))), '-w', MS, 12)
elseif strcmp(xcor, 'y') && strcmp(xcorsetup, 'A')
    hold on
    plot3(X(1:1:end,1), X(1:1:end,2), 1e1 + ones(size(X(1:1:end,1))), 'xw', MS, 12)
else
    hold on
    plot3(X(1:2:end,1), X(1:2:end,2), 1e1 + ones(size(X(1:2:end,1))), 'xw', MS, 12)
end
xticks(''), yticks('')
shg, print(gcf, 'images/lsm.eps', '-depsc', '-r280')