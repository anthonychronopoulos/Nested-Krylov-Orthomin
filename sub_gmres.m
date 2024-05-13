function [solution, flag, relres, iter, ro] = sub_gmres(H,b,K,j,use_ilu)              
    maxiter = K/j;
    tol = 1e-14;

    if use_ilu == 1
        [L, U] = ilu(H); 
        [solution, flag, relres, iter, ro]=gmres(H,b,j,tol,maxiter,L,U);
    else
       n=size(H,1);
       U=speye(n); 
       L=speye(n); 
       [solution, flag, relres, iter, ro]=gmres(H,b,j,tol,maxiter,L,U);
    end
end 

