clear variables;
%Hamiltonian
%H=delta/2 pauli_1 + \sum(H_{SHO}_i) + pauli_3 \sum C_i (x_i)

%define the parameters of system
params.n_osc = 3;
params.n_max = 2;
params.delta = .1;
params.omega = 6;
params.coupling = ones(1,params.n_osc)*05;


%setup empty Hamiltonian and # of possible states
N = ((params.n_max+1).^params.n_osc)*2;
H = zeros(N);

%populate the hamiltonian
for ii = 1:(N/2)
   for jj = 1:(N/2)
       bra_osc_str = dec2base(ii-1,params.n_max +1,params.n_osc);
       ket_osc_str = dec2base(jj-1,params.n_max +1,params.n_osc);
       %top half of states are spin up, bottom half are spin down.
       H(ii,jj) = h_el(1/2,bra_osc_str,1/2,ket_osc_str,params);
       %bra is now spin down since it is bottom half
       H(ii+N/2,jj) = h_el(-1/2,bra_osc_str,1/2,ket_osc_str,params);
       %ket is spin down since column is in bottom half 
       H(ii,jj+N/2) = h_el(1/2,bra_osc_str,-1/2,ket_osc_str,params);
       %both are spin down 
       H(ii+N/2,jj+N/2) = h_el(-1/2,bra_osc_str,-1/2,ket_osc_str,params);
   end
end
% for ii = 1:N
%     for jj = 1:N
% 
%         ket_osc = dec2base(floor((ii-1)./2),params.range+1,params.count);
%         bra_osc = dec2base(floor((jj-1)./2),params.range+1,params.count);
% 
%         if ii <= N/2
%             spin_b = -1/2;
%         else
%             spin_b = 1/2;
%         end
%         
%         if jj <= N/2
%             spin_k = -1/2;
%         else
%             spin_k = 1/2;
%         end
%         H(ii,jj) = h_el(spin_b,bra_osc,spin_k,ket_osc,params);
%     end
% end
%     




%CRITICAL COUPLING
% params.n_osc = 3;
% params.n_max = 2;
% params.delta = 5;
% params.omega = 2;
% params.coupling = ones(1,params.n_osc)*.40;
