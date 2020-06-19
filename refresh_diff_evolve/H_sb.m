function H_couple = H_sb(bath_elements)

x_osc_matrix = @(omega,num_modes) 1./sqrt(2*omega).*(diag(sqrt(1:(num_modes-1)),1)+ diag(sqrt(1:(num_modes-1)),-1));
pauli_z = sparse([1,0;0,-1]);
osc_count = length(bath_elements(3,:));
H_couple = sparse(2*prod(bath_elements(3,:)),2*prod(bath_elements(3,:)));
for ii = 1:osc_count
    
    omega_ii = bath_elements(1,ii);
    tmp = sparse(eye(bath_elements(3,osc_count)));
    if ii == osc_count
        tmp = sparse(x_osc_matrix(omega_ii,bath_elements(3,ii)));
    end
    
    %compute 1x1x...xQ_iix...x1 direct product for bath state space
    for jj = (osc_count-1):-1:1
        if(jj == ii)
            tmp = kron(sparse(x_osc_matrix(omega_ii,bath_elements(3,ii))),tmp);
        else
            tmp = kron(sparse(eye(bath_elements(3,jj))),tmp);
        end
    end
    
    %compute the product with spin space
    H_couple = H_couple + bath_elements(2,ii).*kron(pauli_z,tmp);
    
end
