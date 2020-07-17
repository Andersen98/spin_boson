%problem reduced to real numbers
clear variables;
%energy is in units of Delta so we are left with the
%physics params - (k_b=h_bar=1, energy in units of delta)
params.epsilon = 0;
params.lambda = 2.5;        
params.T = .2; 

%simulation parameters
    %osc_count is the number of oscillators to be used in the bath
    %n_max is the maximal mode for each oscillator i,e, can range from
    %0,1,2,3...n
params.osc_count = 3;
params.mode_count = 2;

%Descretized Debey spectral density - need us to choose a
    %charactaristic frequency for the bath wc
    %core frequency range: min to max [w0,wf]
    %coupling strength (renormalizatin energy) lambda
params.wc= .25;
params.wf = 2*params.wc;
params.w0 = params.wf./(2.*params.osc_count);
params.lambda = 2.5;

%construct the bath terms (frequency and coupling
bath_elements = bath_terms(params.w0,params.wf,params.osc_count,params.wc,params.lambda,params.mode_count);

N = prod(bath_elements(3,:))*2;
%separate the real and imaginary parts
initial_state_re = psi_naught(bath_elements,1/params.T,[1;0]);
initial_state_im = zeros(N,1);
initial_state = [initial_state_re;initial_state_im];
norm(initial_state)
 %H_sb(bath_elements,initial_state)
%H_s(bath_elements,initial_state,params.epsilon)
%H_b(bath_elements,initial_state)
 %H_sb(bath_elements,initial_state)+ H_s(bath_elements,initial_state,params.epsilon) +H_b(bath_elements,initial_state)

 %reduced hamiltonians
 r_Hs = @(yre,yim) [H_s(bath_elements,yim,params.epsilon) ;zeros(N,1)] ...
    -[zeros(N,1);H_s(bath_elements,yre,params.epsilon)] ;
 
r_Hsb = @(yre,yim) [H_sb(bath_elements,yim) ;zeros(N,1)] ...
    -[zeros(N,1);H_sb(bath_elements,yre)] ;

r_Hb = @(yre,yim) [H_b(bath_elements,yim) ;zeros(N,1)] ...
    -[zeros(N,1);H_b(bath_elements,yre)] ;

 FULL_H = @(t,y) r_Hs(y(1:N),y(N+1:2*N,1)) + r_Hb(y(1:N),y(N+1:2*N,1))...
     +r_Hsb(y(1:N),y(N+1:2*N,1));  
[t,y] = ode23s(FULL_H,[0,51],initial_state );

pt = zeros(1,length(t));
for ii = 1:length(t)
    tmp_psi(1:N,1) = y(ii,1:N)' + 1i.*y(ii,N+1:2*N)';
    rho_spin = TrX(tmp_psi,2,[2,prod(bath_elements(3,:))]);
    pt(ii) = trace(rho_spin*[1,0;0,-1]);
end
%plot the solution
plot(t',pt,'-o');
title('Poulation difference Tr[\rho_s \sigma_z] solved with ODE 45')
xlabel('Time t');
ylabel('Population Difference');
