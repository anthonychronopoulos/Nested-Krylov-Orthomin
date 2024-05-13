

function [solution, iter, ro]= program_grcmgsmr(H,b,K,j,use_ilu,method)

solution = NaN; iter = NaN; ro = NaN; % this is to avoid errors when returning
maxiter = K;
ndim = size(H,1);
ztmp = zeros(ndim,1);

h_scfn = zeros(ndim,1);
scfn = zeros(ndim,1);
ro = zeros(maxiter+1,1);
zrmn = b;

i0 = j;
h_phi = zeros(ndim,i0);
phi = zeros(ndim,i0);
total_phi = 0.0;
tol = 1e-14;
nrm2 = norm(zrmn);
%fprintf("initial zrm nrm2 %.16e\n", nrm2);
nrm2_ini = nrm2;

if use_ilu == 1
    [L, U] = ilu(H);
    ztmp = L\zrmn; %solve Ly = b
    zrmn = U\ztmp; %solve Ux = y
    %LU = L*U;
end
for n = 1: maxiter
    %disp(n);
    h_zrmn = H * zrmn;
    if use_ilu == 1
        ztmp = L\h_zrmn;
        h_zrmn = U\ztmp;
    end
    gamma = (zrmn' * h_zrmn) / (h_zrmn' * h_zrmn);
    kspn = 2.0 * gamma * zrmn - gamma*gamma * h_zrmn;
    
    h_kspn = H * kspn;
    if use_ilu == 1
        ztmp = L\h_kspn;
        h_kspn = U\ztmp;
    end
    
    %MGS LOOP
    for i = 1: i0-1
        ntmpi = n+i-i0;
        if ntmpi<1; continue; end
        dot = h_kspn' * h_phi(:,mod(ntmpi,i0)+1);
        h_kspn = h_kspn - dot * h_phi(:,mod(ntmpi,i0)+1);
        kspn = kspn - dot * phi(:,mod(ntmpi,i0)+1);
    end
    
    nrm_h_kspn = norm(h_kspn);
    if nrm_h_kspn <= 1e-20
        fprintf("nrm_h_kspn equals zero  %.16e\n", nrm_h_kspn);
        return;
    end
    kspn = kspn / nrm_h_kspn;
    h_kspn = h_kspn / nrm_h_kspn;
    dot = zrmn' * h_kspn;
    scfn = dot * kspn;
    h_scfn = dot * h_kspn;
    
    nrm_h_scfn = norm(h_scfn);
    if nrm_h_scfn <= 1e-20
        fprintf("nrm_h_scfn h_phi equals zero  %.16e\n", nrm_h_scfn);
        return;
    end
    phi(:,mod(n,i0)+1) = scfn / nrm_h_scfn;
    h_phi(:,mod(n,i0)+1) = h_scfn / nrm_h_scfn;
    
    %update solution and residual
    total_phi = total_phi + scfn;
    zrmn = zrmn - h_scfn;
    
    %convergence check
    nrm2 = norm(zrmn)/nrm2_ini;
    ro(n) = nrm2;
    if tol >= nrm2
        break;
    end
end
iter = n;
solution = total_phi;
end

