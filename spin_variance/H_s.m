function H_spin = H_s(n_max,bath_elements,epsilon)
%H_S outputs the spin only operators of the Hamiltonian
%   H_s = epsilon*sigma_z + sigma_x

%define variables
mode_count = n_max +1;
osc_count = length(bath_elements(1,:));
n_bath = mode_count^osc_count;
b_ones = ones(1,n_bath);
n_total = 2*n_bath;
H_spin = zeros(n_total,n_total);
%top left
H_spin(1:n_bath,1:n_bath) = diag(epsilon.*b_ones);
%top right
H_spin(1:n_bath,(n_bath+1):n_total)= diag(b_ones);
%botom left
H_spin((n_bath+1):n_total,1:n_bath)= diag(b_ones);
%bottom right
H_spin((n_bath+1):n_total,(n_bath+1):n_total) = diag(-1.*epsilon.*b_ones);
%;diag(b_ones),diag(-1.*epsilon.*b_ones)];
end

