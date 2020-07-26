function psi_prime = H_sb(bath_elements,psi)
%H_sb computes the spin boson hamiltonian acting on the state psi

osc_count = length(bath_elements(1,:));
bath_count = prod(bath_elements(3,:));
psi_prime = zeros(bath_count*2,1);
mode_count = bath_elements(3,1);
n_max = mode_count - 1;
%take care of truncated terms first
%suppose that the mode is k = 3
kk = 3;

%we can iterate through a manifold where all number values are fixed
%except k with this loop
number_state = zeros(osc_count,1);
for ii = 1:bath_count
   
    for jj = 1:osc_count
       if(number_state(jj) ~= mode_count)
           break;
       end
       number_state(jj+1) = number_state(jj+1) + 1;
       number_state(jj) = 0;
    end
    
    
    %iterate through osc terms
    for kk = 1:osc_count
        offset = mode_count^(kk-1);
        n_k = number_state(kk);
        q_k = 1/sqrt(2*bath_elements(1,kk));
        c_k = bath_elements(2,kk);
        if n_k < n_max
           psi_prime(ii) = psi_prime(ii) + c_k*psi(ii + offset).*sqrt(n_k+1)*q_k;
           psi_prime(ii+bath_count) = psi_prime(ii+bath_count) - c_k*psi(ii + offset + bath_count).*sqrt(n_k+1)*q_k;
        end
           
        if number_state(kk) > 0
            psi_prime(ii) = psi_prime(ii) + c_k*psi(ii - offset).*sqrt(n_k)*q_k;
            psi_prime(ii+bath_count) = psi_prime(ii+bath_count) -c_k*psi(ii - offset+bath_count).*sqrt(n_k)*q_k;
        end       
    end
    
    
    number_state(1,1) = number_state(1,1) + 1;
end
       
end

