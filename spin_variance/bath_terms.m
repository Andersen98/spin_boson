function bath= bath_terms(w0,wf,f,wc,eta)
%BATH_TERMS produces 
    %row 1 - frequencies of each harmonic osc 
    %row 2 - Debye coupling terms
    bath = zeros(2,f);
    bath(1,:) = linspace(w0,wf,f);
    J = @(w) 2.*eta.*wc.*(w./(w.^2+wc.^2)); %eq35 berkelbach 2011
    rho = @(w) (f./(pi.*wc)).*2./(1+(w./wc).^2); %eq 38 berkelback 2011
    bath(2,:) = 2/pi.*bath(1,:).*J(bath(1,:))./rho(bath(1,:));
    bath(2,:) = sqrt(bath(2,:));
    %eta_check = .5.*sum(((bath(2,:).^2)./(bath(1,:).^2)))
end

