function initial_state = psi_naught(bath_elements,beta)
%Makes initial state in number representation

osc_count = size (bath_elements);
osc_count = osc_count(2);
bath_count = prod(bath_elements(3,:));
n_max = bath_elements(3,1)-1;
full_base_vec = [bath_elements(3,:),bath_count,2];
%initialize state vector
initial_state = zeros(2*bath_count*bath_count,1);
state_norm = 0;
for ii = 1:bath_count
    number_state = idx2num(ii-1, bath_elements(3,:));   
    %energy is just w_ii*(n_ii + 1/2)
    energies = number_state + 1/2;
    
    energies = prod([energies;bath_elements(1,:)],1);
    %c_ii = sqrt(prob(state_ii))
    c_ii = exp (-.5.*beta.*sum(energies));
    state_norm = state_norm + exp(-1*beta.*sum(energies));
    
    %get the full indexing for purified state
    full_state = [number_state,ii-1,0];
    full_idx = num2idx(full_state, full_base_vec);
    initial_state(full_idx+1) = c_ii;
end

%normalize state
state_norm = sqrt(state_norm);
for ii = 1:bath_count
    number_state = idx2num(ii-1, bath_elements(3,:));   
    %get the full indexing for purified state
    full_state = [number_state,ii-1,0];
    full_idx = num2idx(full_state, full_base_vec);
    initial_state(full_idx+1) = initial_state(full_idx+1)./state_norm;
end

end

