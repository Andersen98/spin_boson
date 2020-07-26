function psi_prime = anc_H(H_s,H_b,H_sb, psi, num_bath)
%ANC_H Summary of this function goes here
%   Detailed explanation goes here
psi_prime = zeros(2*num_bath*num_bath,1);
H_func = @(psi) H_s(psi) + H_b(psi) + H_sb(psi);
for k = 0:(num_bath-1)
    idx_start = k*num_bath+1;
    idx_end = (k+1)*(num_bath);
    psi_tmp = [psi(idx_start:idx_end,1);psi(idx_start+num_bath^2:idx_end+num_bath^2,1)];
    output_psi_prime =  H_func(psi_tmp);
    psi_prime(idx_start:idx_end,1) = psi_prime(idx_start:idx_end,1) + output_psi_prime(1:num_bath,1);
    psi_prime(idx_start+num_bath^2:idx_end+num_bath^2,1) = psi_prime(idx_start+num_bath^2:idx_end+num_bath^2,1) + output_psi_prime(num_bath+1:2*num_bath,1);
end
    
end

