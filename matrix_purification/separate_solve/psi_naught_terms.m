function terms = psi_naught_terms(bath_elements,beta)
%Makes initial state in number representation


num_bath = prod(bath_elements(3,:));
terms = zeros(1,num_bath);
state_norm = 0;
for ii = 1:num_bath
    number_state = idx2num(ii-1, bath_elements(3,:));   
    %energy is just w_ii*(n_ii + 1/2)
    energies = number_state + 1/2;
    
    energies = prod([energies;bath_elements(1,:)],1);
    %c_ii = sqrt(prob(state_ii))
    c_ii = exp (-.5.*beta.*sum(energies));
    state_norm = state_norm + exp(-1*beta.*sum(energies));
    
    %assign value
    terms(ii) = c_ii;
end

%normalize terms
state_norm = sqrt(state_norm);
terms = terms./state_norm;
end

