function MI1 = MI1_dis(probs, q_probs)

%-------------------------------------------------------------------------------------------------
% PURPOSE: calculate mutual (multi) information between subsystems in discretized data
%
% INPUTS:
%   probs: probability distributions computed from X
%       probs.p: probability distribution of X only used for mutual
%       information (MI)
%   q_probs: mismatched probability distributions q
%       q_probs.q: mismtached probability distribution of X
%
% OUTPUTS:
%   MI1: mutual (multi) information between subsystems 
%-------------------------------------------------------------------------------------------------
%
% Masafumi Oizumi, 2018

p = probs.p;
q = q_probs.q;

N_c = length(q);
H_part = zeros(N_c,1);
for i=1: N_c
    H_part(i) = H_dis(q{i});
end
H_joint = H_dis(p);

MI1 = sum(H_part) - H_joint;

end