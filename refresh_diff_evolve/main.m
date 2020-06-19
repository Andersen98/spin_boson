clear variables;
%energy is in units of Delta so we are left with the
%physics params - (k_b=h_bar=1, energy in units of delta)
params.epsilon = 0;
params.lambda = 2.5;
params.T = 2; 

%simulation parameters
    %osc_count is the number of oscillators to be used in the bath
    %n_max is the maximal mode for each oscillator i,e, can range from
    %0,1,2,3...n
params.n_max = 1;
params.osc_count = 5;

%Descretized Debey spectral density - need us to choose a
    %charactaristic frequency for the bath wc
    %core frequency range: min to max [w0,wf]
    %coupling strength (renormalizatin energy) lambda
params.wc= 4;
params.wf = 2*params.wc;
params.w0 = params.wf./(2.*params.osc_count);
params.lambda = 2.5;

%construct the bath terms (frequency and coupling
bath_elements = bath_terms(params.w0,params.wf,params.osc_count,params.wc,params.lambda);

%get a truncated initial bath state, where the bath is in "thermal
%equilibrium"


%make an initial state that would give a initial density matrix \rho(0) in
%section IV of Berkelbach 2011
initial_state = psi_naught(bath_elements,1/params.T,[1,0]);


[t,y] = ode45(@(t,y) H_s(t,y,params.n_max,bath_elements,params.epsilon),[0,12],initial_state);

%helper variables
mode_count = params.n_max +1;
osc_count = length(bath_elements(1,:));
n_bath = mode_count^osc_count;
n_total = 2*n_bath;
pt = zeros(length(t),1);
for ii = 1:length(t)
    rho_spin = TrX(y(ii,:)',2,[2,n_bath]);
    %rho_spin = [sum(y(ii,1:n_bath).^2),0;0,sum(y(ii,n_bath+1:n_total).^2)];
    pt(ii) = trace(rho_spin*[1,0;0,-1]);
end
%plot the solution
plot(t,pt(:,1),'-o');
title('Solution site population (\X-Dimentional) with ODE 45')
xlabel('Time t');
ylabel('Solution y');
legend('y_1','y_2')