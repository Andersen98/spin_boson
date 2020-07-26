function psi_prime = H_s(bath_elements,psi,epsilon)

osc_count = size (bath_elements);
osc_count = osc_count(2);
bath_count = prod(bath_elements(3,:));
n_max = bath_elements(3,1)-1;
number_state = zeros(1,osc_count);

psi_prime = zeros(bath_count*2,1);
psi_prime(1:bath_count,1) = epsilon.* psi(1:bath_count,1);
psi_prime(bath_count + 1: 2*bath_count,1) = -1.*epsilon.*psi(bath_count +1: 2*bath_count,1);

psi_prime(1:bath_count,1) = psi_prime(1:bath_count,1) + psi(bath_count + 1:2* bath_count,1);
psi_prime(bath_count+1:2*bath_count,1) = psi_prime(bath_count + 1:2*bath_count,1) + psi(1:bath_count,1);


end