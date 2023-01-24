function u = nearFieldPaper(z, k, X, ui)

% Parameters:
typ = 'P1';   % type of elements
N = 5e1;      % number of elements
tol = 1e-3;   % error tolerance
gss = 3;      % number of Gauss points

% Boundary mesh:
h = (2*pi)/N;
t = (0:h:(2*pi-h))';
vtx = [real(z(t)), imag(z(t)), zeros(N, 1)];
elt = [[(2:N)'; 1], (1:N)'];
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