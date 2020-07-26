function psi_prime = H_b(bath_elements,psi)
%H_b applies the bath hamiltonian to the state psi


osc_count = size (bath_elements);
osc_count = osc_count(2);
bath_count = prod(bath_elements(3,:));
n_max = bath_elements(3,1)-1;
number_state = zeros(1,osc_count);
psi_prime = zeros(bath_count*2,1);


for ii = 1:bath_count
   
    %iterate through number states by idx
    
    for jj = 1:osc_count
        if number_state(jj) <= n_max
            break;
        end
        number_state(jj) = 0;
        number_state(jj+1) = number_state(jj+1)+ 1;
    end    
    %energy is just w_ii*(n_ii + 1/2)
    energies = number_state + 1/2;
    energies = prod([energies;bath_elements(1,:)],1);
    
    %c_ii = sqrt(prob(state_ii))
    c_ii = sum(energies);
    psi_prime(ii) = c_ii*psi(ii);
    psi_prime(ii+bath_count) = c_ii*psi(ii+bath_count);
    
    %increment number state
    number_state(1) = number_state(1) + 1;
end

end

