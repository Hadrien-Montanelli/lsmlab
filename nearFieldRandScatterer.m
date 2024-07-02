function u = nearFieldRandScatterer(z, k, X, ui, Y, epsilon, nobstacle)

% Parameters:
typ = 'P1';     % type of elements
N = 5e1;        % number of elements
tol = 1e-3;     % error tolerance
gss = 3;        % number of Gauss points
Q = size(Y, 1); % number of small scatterers

% Boundary mesh:
lambda = 2*pi/k;
h = (2*pi)/N;
t = (0:h:(2*pi-h))';
vtx = [real(z(t)), imag(z(t)), zeros(N, 1)];
elt = [[(2:N)'; 1], (1:N)'];
if nobstacle == 1
    for q = 1:Q
        vtx = [vtx; lambda*epsilon*[cos(t), sin(t), zeros(N, 1)] + Y(q, :)];
        elt = [elt; [[(2:N)'; 1] (1:N)'] + q*N];
    end
elseif nobstacle == 2
    vtx = [vtx; -real(z(t)), -imag(z(t)), zeros(N, 1)];
    elt = [elt; [[(2:N)'; 1] (1:N)'] + N];
    for q = 1:Q
        vtx = [vtx; lambda*epsilon*[cos(t), sin(t), zeros(N, 1)] + Y(q, :)];
        elt = [elt; [[(2:N)'; 1] (1:N)'] + (q+1)*N];
    end
end
mesh = msh(vtx, elt);

% Frequency adjusted to maximum esge size:
stp = mesh.stp;
kmax = 1/stp(2);
if (k > kmax)
    warning('Wave number is too high for mesh resolution.')
end

% Green kernel:
Gxy = @(X, Y) femGreenKernel(X, Y, '[H0(kr)]', k);

% Domain:
sigma = dom(mesh, gss);

% Finite elements:
u = fem(mesh, typ);
v = fem(mesh, typ);

% Boundary element operator:
S = (1i/4).*integral(sigma, sigma, u, Gxy, v, tol);

% Regularization:
Sr = -1/(2*pi).*regularize(sigma, sigma, u, '[log(r)]', v);

% Final operator:
LHS = S + Sr;

% Finite element incident wave trace:
RHS = -integral(sigma, u, ui);

% LU factorization:
[Lh, Uh] = lu(LHS);

% Solve linear system:
psi = Uh\(Lh\RHS);

% Compute u:
Snear = 1i/4.*integral(X, sigma, Gxy, v, tol);
Sreg = -1/(2*pi).*regularize(X, sigma, '[log(r)]', v);
Snear = Snear + Sreg;
u = Snear*psi;

end