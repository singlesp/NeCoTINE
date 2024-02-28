function [E,E_nodal] = MIN_CONTROL_ENERGY(A, WcI, x0, xf, T,normalize)
% Computes minimum control energy for state transition.
% A: System adjacency matrix: n x n
% x0: Initial state
% xf: Final state
% T: Control horizon
% 
% Outputs
% E: Minimum control energy
if ~exist('normalize','var')
	normalize = true;
end

% Normalize
if normalize
	A = (A / (max(eig(A)))) - eye(length(A));
	disp(['After normalization, max eigenvalue of A is ',num2str(max(eig(A)))])
end
% State transition to achieve
Phi = expm(A*T)*x0 - xf;
% Energy - Addition of nodal output by SPS 03/01/23
E_nodal = (WcI*Phi).*Phi;
E = sum(E_nodal);
% E_nodal = E_nodal'; %transpose so that output is same size as MIN_END_CONT.m
end