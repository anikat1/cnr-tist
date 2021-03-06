function [ Z, R, funVal, iteration ] = OSC_fast( X, lambda_1, lambda_2, gamma_1, gamma_2, p, maxIterations, diagconstraint)
%% Problem
%
% min L(Z) = 1/2 ||X - XZ||^2_F - lambda_1 * ||Z||_1 + lambda_2 * ||ZR||_2/1
% 
% where ||B||_2/1 = ||b_1||_2 + ||b_2||_2 + ... + ||b_n||_2
%
%% Solution 
% We solve this problem via the ADMM (Alternating direction method of multipliers) variant
% of Augmented Lagrangian method as follows:
% 
% Let S = Z and U = SR then we have
% 
% T(Z,S,U) = 1/2||X - XS||^2_F + lambda_1 * ||Z||_1 + lambda_2 * ||U||_2/1
% + <G,Z - S> + gamma_1/2 * ||Z - S||^2_F + <F,U - SR>
% + gamma_2/2 * ||U - SR||^2_F
%

if (~exist('diagconstraint','var'))
    diagconstraint = 0;
end

funVal = zeros(maxIterations,1);

[~, xn, ~] = size(X);

S = zeros(xn, xn); % S = Z
R = (triu(ones(xn,xn-1),1) - triu(ones(xn, xn-1))) + (triu(ones(xn, xn-1),-1)-triu(ones(xn, xn-1)));
R = sparse(R);

U = zeros(xn, xn-1);

G = zeros(xn, xn);
F = zeros(xn, xn-1); 

Z = zeros(xn, xn);

for iteration=1:maxIterations
    %fprintf('Iteration Number %d \n',iteration);
    %% Step 1
    V = S - (G/gamma_1);
    [vm, vn, ~] = size(V);
    rolled_v = reshape(V,vm*vn,1);

    rolled_z = shrink_l1(rolled_v, lambda_1/gamma_1);

    Z = reshape(rolled_z, vm, vn);

    % Set Z diag to 0
    if (diagconstraint)
        Z(logical(eye(size(Z)))) = 0;
    end

    %% Step 2
    A = X'*X + gamma_1*speye(xn,xn);
    B = gamma_2*(R*R');
    C = -(X'*X + gamma_2*U*R' + gamma_1*Z + G + F*R');

    S = lyap(A, B, C);

    %% Step 3 
    V = S*R - (1/gamma_2)*F;

    U = mysolve_l1l2(V, lambda_2/gamma_2);

    %% Step 4 update G

    G = G + gamma_1 * (Z - S);

    %% Step 5 update G

    F = F + gamma_2 * (U - S*R);

    %% Step 6 update gamma params

    gamma_1 = p * gamma_1;
    gamma_2 = p * gamma_2;

    %% Calculate function values
    funVal(iteration) = .5 * norm(X - X*Z,'fro')^2 + lambda_1*norm(Z,1) +lambda_2*l2l1norm(Z*R); %U=ZR
	fprintf('Iteration %d func val: %.3f\n',iteration,funVal(iteration));
    if iteration > 1
        if funVal(iteration) < 1*10^-3
			fprintf('fro norm %.3f V_1 norm %.3f VR %.3f\n',norm(X - X*Z,'fro')^2,norm(Z,1),l2l1norm(Z*R));
            break
        end
    end
    %fprintf('\nfunction value %.12f, iteration %d\n',funVal(iteration),iteration);
    if iteration > 100
		if funVal(iteration) < 1*10^-4 || funVal(iteration-1) == funVal(iteration) ...
                || funVal(iteration-1) - funVal(iteration) < 1*10^-4
			fprintf('fro norm %.3f V_1 norm %.3f VR %.3f\n',norm(X - X*Z,'fro')^2,norm(Z,1),l2l1norm(Z*R));
            break
        end
    end
    

end

end