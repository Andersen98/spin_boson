function initial_state = psi_naught(bath_elements,beta,spin_state)
%Makes initial state in number representation

osc_count = size (bath_elements);
osc_count = osc_count(2);
bath_count = prod(bath_elements(3,:));
n_max = bath_elements(3,1)-1;
number_state = zeros(1,osc_count);

%initialize state vector
initial_state = zeros(2*bath_count,1);
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
    c_ii = exp (-.5.*beta.*sum(energies));
     
    initial_state(ii) = spin_state(1)*c_ii;
    initial_state(ii+bath_count) = spin_state(2)*c_ii;
    
    %increment number state
    number_state(1) = number_state(1) + 1;
end

initial_state = 1/norm(initial_state).*initial_state;
end

