function x = conjugateGradient(A, b, x_0, tol)
% This function takes as input a matrix A, RHS vector b, and initial guess
% x_0 and performs conjugate gradient steps to solve the linear system to a
% specified tolerance. For now, matrix-vector products are simply handled
% by MATLAB, under the assumption that it already has decent routines for
% this when A is sparse (a quick Google suggests it does)

% Initialise
x = x_0;

% Calculate the residual
r = b - A * x;

% Initialise direction as the residual
d = r;

% Store square residual
rsold = r' * r;

% Now loop over finding optimal step length in this direction, and then
% calculating a new search direction, until tolerance is met
iters = 1;
while sqrt(rsold) > tol && iters < 1000
    
    % Store matrix vector product
    Ad = A * d;
    
    % Calculate optimal step length
    alpha = rsold / ( d' * Ad );
    
    % Update estimated solution
    x = x + alpha * d;
    
    % Update residual
    r = r - alpha * Ad;
    
    % Update residual estimate
    rsnew = r' * r;
    
    % Find new search direction
    d = r + (rsnew / rsold) * d;
    
    % Now new information becomes old for start of new loop
    rsold = rsnew;
    
    % Update iteration counter
    iters = iters + 1;

end

