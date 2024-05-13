%%% MAIN PROGRAM  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%PROGRAM AUTHORS: D. Slavchev(dimitargslavchev@abv.bg) , 
%%T. Abe  (abe.toshihiko24@gmail.com), A. T. Chronopoulos (antony.tc@gmail.com) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% PROBLEM SOLVED
%% sparse linear system A x=b, with A (non)symmetric pos. definite %% 
%% Programs: (1) main_all, 
%% (2) mmread (from https://math.nist.gov/MatrixMarket/)
%% (3)program_grcmgs (4) program_grcmgsmr (5) program_orthomin 
%% (6) sub_gmres    (https://www.mathworks.com/help/matlab/ref/gmres.html)  
%%% METHODS IMPLEMENTED: 
%%%(3), (4): GRC or Nested-Orthomin(Nested Krylov) method with MGS
%% (5) Orthomin method with MGS 
%% Sparse Matrix in (MM)-format %%
%% example: raefsky2.mtx %%
%% rhs b is created to have solution vector of all entries=1  %%
%% initial solution vector =0  %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Nested Krylov (GRC) (programmed in C) are presented in the  paper:
%%Abe, T. and Chronopoulos, A.T., 2023. The generalized residual cutting %%method and its convergence characteristics.Numerical Linear Algebra %%with Applications,30(6), p.e2517. 

matrix_filename = 'raefsky2.mtx';
A = mmread(matrix_filename); 

K = 10000;     % number of iterations.
j = 10;         %j used in GRC-MGS(j), GMRES(j), Othomin(j)

clc

n=size(A,1);
% create rhs b, so that x1 is the solution
x1=ones(n,1); % known solution
b =A*x1;


j = 5;
method = 'orthomin';
use_ilu = 0;
tic;
[solution, iter, ro]=program_orthomin(A,b,K,j,use_ilu,method);
if use_ilu == 1; precon = 'ILU(0)'; else precon = 'None  ';end
et=toc;
fprintf('precondition  method     j    iter    trueError     NormalizedResidualError ElapsedTime\n')
fprintf('%s        %s   %d    %5d   %e  %e             %f\n', precon, method,j, iter, norm(solution-x1), norm(b-A*solution)/norm(b), et)

use_ilu = 1;
tic;
[solution, iter, ro]=program_orthomin(A,b,K,j,use_ilu,method);
if use_ilu == 1; precon = 'ILU(0)'; else precon = 'None  ';end
et=toc;
fprintf('precondition  method     j    iter    trueError     NormalizedResidualError ElapsedTime\n')
fprintf('%s        %s   %d    %5d   %e  %e             %f\n', precon, method,j, iter, norm(solution-x1), norm(b-A*solution)/norm(b), et)
fprintf('\n');

j = 10;
method = 'orthomin';
use_ilu = 0;
tic;
[solution, iter, ro]=program_orthomin(A,b,K,j,use_ilu,method);
if use_ilu == 1; precon = 'ILU(0)'; else precon = 'None  ';end
et=toc;
fprintf('precondition  method     j    iter    trueError     NormalizedResidualError ElapsedTime\n')
fprintf('%s        %s   %d    %5d   %e  %e             %f\n', precon, method,j, iter, norm(solution-x1), norm(b-A*solution)/norm(b), et)

use_ilu = 1;
tic;
[solution, iter, ro]=program_orthomin(A,b,K,j,use_ilu,method);
if use_ilu == 1; precon = 'ILU(0)'; else precon = 'None  ';end
et=toc;
fprintf('precondition  method     j    iter    trueError     NormalizedResidualError ElapsedTime\n')
fprintf('%s        %s   %d    %5d   %e  %e             %f\n', precon, method,j, iter, norm(solution-x1), norm(b-A*solution)/norm(b), et)
fprintf('\n');

j = 3;
method = 'grcmgs';
use_ilu = 0;
tic;
[solution, iter, ro]=program_grcmgs(A,b,K,j,use_ilu,method);
if use_ilu == 1; precon = 'ILU(0)'; else precon = 'None  ';end
et=toc;
fprintf('precondition  method     j    iter    trueError     NormalizedResidualError ElapsedTime\n')
fprintf('%s        %7s   %d    %5d   %e  %e             %f\n', precon, method,j, iter, norm(solution-x1), norm(b-A*solution)/norm(b), et)

use_ilu = 1;
tic;
[solution, iter, ro]=program_grcmgs(A,b,K,j,use_ilu,method);
if use_ilu == 1; precon = 'ILU(0)'; else precon = 'None  ';end
et=toc;
fprintf('precondition  method     j    iter    trueError     NormalizedResidualError ElapsedTime\n')
fprintf('%s        %7s   %d    %5d   %e  %e             %f\n', precon, method,j, iter, norm(solution-x1), norm(b-A*solution)/norm(b), et)
fprintf('\n');

j = 5;
method = 'grcmgs';
use_ilu = 0;
tic;
[solution, iter, ro]=program_grcmgs(A,b,K,j,use_ilu,method);
if use_ilu == 1; precon = 'ILU(0)'; else precon = 'None  ';end
et=toc;
fprintf('precondition  method     j    iter    trueError     NormalizedResidualError ElapsedTime\n')
fprintf('%s        %7s   %d    %5d   %e  %e             %f\n', precon, method,j, iter, norm(solution-x1), norm(b-A*solution)/norm(b), et)

use_ilu = 1;
tic;
[solution, iter, ro]=program_grcmgs(A,b,K,j,use_ilu,method);
if use_ilu == 1; precon = 'ILU(0)'; else precon = 'None  ';end
et=toc;
fprintf('precondition  method     j    iter    trueError     NormalizedResidualError ElapsedTime\n')
fprintf('%s        %7s   %d    %5d   %e  %e             %f\n', precon, method,j, iter, norm(solution-x1), norm(b-A*solution)/norm(b), et)
fprintf('\n');

j = 3;
method = 'grcmgsmr';
use_ilu = 0;
tic;
[solution, iter, ro]=program_grcmgsmr(A,b,K,j,use_ilu,method);
if use_ilu == 1; precon = 'ILU(0)'; else precon = 'None  ';end
et=toc;
fprintf('precondition  method     j    iter    trueError     NormalizedResidualError ElapsedTime\n')
fprintf('%s        %s   %d    %5d   %e  %e             %f\n', precon, method,j, iter, norm(solution-x1), norm(b-A*solution)/norm(b), et)

use_ilu = 1;
tic;
[solution, iter, ro]=program_grcmgsmr(A,b,K,j,use_ilu,method);
if use_ilu == 1; precon = 'ILU(0)'; else precon = 'None  ';end
et=toc;
fprintf('precondition  method     j    iter    trueError     NormalizedResidualError ElapsedTime\n')
fprintf('%s        %s   %d    %5d   %e  %e             %f\n', precon, method,j, iter, norm(solution-x1), norm(b-A*solution)/norm(b), et)
fprintf('\n');

j = 5;
method = 'grcmgsmr';
use_ilu = 0;
tic;
[solution, iter, ro]=program_grcmgsmr(A,b,K,j,use_ilu,method);
if use_ilu == 1; precon = 'ILU(0)'; else precon = 'None  ';end
et=toc;
fprintf('precondition  method     j    iter    trueError     NormalizedResidualError ElapsedTime\n')
fprintf('%s        %s   %d    %5d   %e  %e             %f\n', precon, method,j, iter, norm(solution-x1), norm(b-A*solution)/norm(b), et)

use_ilu = 1;
tic;
[solution, iter, ro]=program_grcmgsmr(A,b,K,j,use_ilu,method);
if use_ilu == 1; precon = 'ILU(0)'; else precon = 'None  ';end
et=toc;
fprintf('precondition  method     j    iter    trueError     NormalizedResidualError ElapsedTime\n')
fprintf('%s        %s   %d    %5d   %e  %e             %f\n', precon, method,j, iter, norm(solution-x1), norm(b-A*solution)/norm(b), et)
fprintf('\n');


j = 10;
method = 'gmres';
use_ilu = 0;
tic;
[solution, flag, relres, iter, ro]=sub_gmres(A,b,K,j,use_ilu);
if use_ilu == 1; precon = 'ILU(0)'; else precon = 'None  ';end
et=toc;
iter = (iter(1)-1)*j+iter(2);
fprintf('precondition  method     j    iter    trueError     NormalizedResidualError ElapsedTime\n')
fprintf('%s        %7s   %d    %5d   %e  %e             %f\n', precon, method,j, iter, norm(solution-x1), norm(b-A*solution)/norm(b), et)

use_ilu = 1;
tic;
[solution, flag, relres, iter, ro]=sub_gmres(A,b,K,j,use_ilu);
if use_ilu == 1; precon = 'ILU(0)'; else precon = 'None  ';end
et=toc;
iter = (iter(1)-1)*j+iter(2);
fprintf('precondition  method     j    iter    trueError     NormalizedResidualError ElapsedTime\n')
fprintf('%s        %7s   %d    %5d   %e  %e             %f\n', precon, method,j, iter, norm(solution-x1), norm(b-A*solution)/norm(b), et)
fprintf('\n');

j = 40;
method = 'gmres';
use_ilu = 0;
tic;
[solution, flag, relres, iter, ro]=sub_gmres(A,b,K,j,use_ilu);
if use_ilu == 1; precon = 'ILU(0)'; else precon = 'None  ';end
et=toc;
iter = (iter(1)-1)*j+iter(2);
fprintf('precondition  method     j    iter    trueError     NormalizedResidualError ElapsedTime\n')
fprintf('%s        %7s   %d    %5d   %e  %e             %f\n', precon, method,j, iter, norm(solution-x1), norm(b-A*solution)/norm(b), et)

use_ilu = 1;
tic;
[solution, flag, relres, iter, ro]=sub_gmres(A,b,K,j,use_ilu);
if use_ilu == 1; precon = 'ILU(0)'; else precon = 'None  ';end
et=toc;
iter = (iter(1)-1)*j+iter(2);
fprintf('precondition  method     j    iter    trueError     NormalizedResidualError ElapsedTime\n')
fprintf('%s        %7s   %d    %5d   %e  %e             %f\n', precon, method,j, iter, norm(solution-x1), norm(b-A*solution)/norm(b), et)

