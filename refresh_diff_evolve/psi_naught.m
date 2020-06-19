function psi = psi_naught(bath_terms,spin_state, beta)

mode_counts = bath_terms(3,:);
omega_vec = bath_terms(1,:);
osc_count = length(omega_vec);
n_bath_terms = prod(mode_counts);

%construct the bath states
psi = zeros(1,n_bath_term);
counter_vec = zeros(1,osc_count);
for ii = 1:n_bath_terms
     %increment counter by one and then iter. through terms
     counter_vec(1) = counter_vec(1) + 1;
    for jj = 1:1:osc_count
      if mod(counter_vec(jj),mode_count(jj)) == 0
          counter_vec(jj) = 0;
          counter_vec(jj+1) = counter_vec(jj+1) +1;
          
      end 
    end
    counter_vec
    tmp_state = zeros(mode_counts(jj),1);
    tmp_state(counter_vec(jj)) = 1;
    tmp_energy = bath_terms(1,jj).*(counter_vec+1/2);
    therm_prob = exp(-beta.*tmp_energy);
    tmp_state = sqrt(therm_prob).*tmp_state;
    psi = tmp_state;
    for jj = 2:1:osc_count
        
        tmp_state = zeros(mode_counts(jj),1);
        tmp_state(counter_vec(jj)) = 1;
        tmp_energy = bath_terms(1,jj).*(counter_vec+1/2);
        therm_prob = exp(-beta.*tmp_energy);
        tmp_state = sqrt(therm_prob).*tmp_state;
        
        psi = kron(tmp_state,psi);
    end

end

psi = 1/norm(psi).*psi;
psi = kron(spin_state,psi);

end