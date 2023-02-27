
% Image version Conjugate Gradient method to solve 
% A(x)=b
% here, A is an linear operator,
% x0 is the initial guess
% b is known vector in the size of an image
% f is unknown, also in the size of an image.
function [ x ] = CG(Op, b, x0, tol, nloop)

if nargin < 4
    tol = 1e-5;
    nloop = 50;
end
     % Initialization
     r0 = Op(x0) - b; 
     p = -r0; 
     x=x0; 
     err0 = realmax; 
     if norm(r0(:), 2)<tol
           return;
     end
     for iter=1:nloop
         Ap = Op(p); 
         alpha = r0(:)'*r0(:)/(p(:)'*Ap(:)+eps); 
         x = x + alpha*p; 
         r1 = r0 + alpha * Ap; 
         beta = r1(:)'*r1(:)/(r0(:)'*r0(:)+eps);
         p = -r1 + beta*p; 
         r0=r1; 
         
         % stopping condition
         err1 = norm(r1(:), 2);
%          fprintf('CG: iteration %d, rhs: %f\n', iter, err1); 
         if err1<tol
           break;
         end
         if err1 > err0
             nloop = 10;
         else
             err0 = err1; 
         end
%          iter
%           disp(err1)
     end
