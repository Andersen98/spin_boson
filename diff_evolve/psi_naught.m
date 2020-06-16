function psi = psi_naught(n_max,bath_elements,beta)
%PSI_NAUGHT given harmonic osc energy and inverse temp, generates initial
%global state
%   makes state that gives rho(o) = |1><1|exp(-BetaHb)/Z_b

%define variables
mode_count = n_max +1;
osc_count = length(bath_elements(1,:));
n_bath = mode_count^osc_count;
n_total = 2*n_bath;
psi = zeros(n_total,1);
my_ket = zeros(1,osc_count);
for ii = 1:n_bath
    dec = ii-1;
    %change the base into radix = mode_count
    for jj = 1:osc_count
        quo = floor(dec./mode_count);
        rem = mod(dec,mode_count);
        my_ket(jj) = rem;
        dec = quo;
    end
    psi_energy_ii =  sum(bath_elements(1,:).*my_ket) + sum(bath_elements(1,:)./2);
    psi_prob_ii = exp(psi_energy_ii*-beta); 
    psi(ii) = sqrt(psi_prob_ii);
    
end
psi = (1/norm(psi)).*psi;

end
