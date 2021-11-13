function [eta, eta1, eta2, eta3,eta4,eta5] = StokesEta(k, h, a, theta)

% here the a is the first order first harmonic wave amplitude, not
% \widetilde{a} in the manuscript (not aw in the code). 
% Note: 
% result = StokesEta(k, h, a, theta) % will only give the final results eta
    sigma = tanh(k*h);
    alpha1 = cosh(2*k*h);
    % first and second order
    eta1 = a*cos(theta); 
    eta2 = k*a^2/4*(3-sigma^2)/sigma^3 * cos(2*theta);
    % Dingemans corrected, 3rd order
    B31 = (3+8*sigma^2-9*sigma^4)/16/sigma^4;
    B33 = (27-9*sigma^2+9*sigma^4-3*sigma^6)/64/sigma^6;
    eta3 = k^2*a^3*(B31*cos(theta)+B33*cos(3*theta));
    % our derived 4th order equation
    sigma1 =24*(3*alpha1+2)*(alpha1-1)^4*sinh(2*k*h);
    eta4 = k^3*a^4/sigma1.*((60*alpha1^6 ...
    +232*alpha1^5-118*alpha1^4-989*alpha1^3-607*alpha1^2 ...
    +352*alpha1+260).*cos(2*theta)+(24*alpha1^6 ...
    +116 *alpha1^5+214*alpha1^4+188*alpha1^3+133*alpha1^2 ...
    +101*alpha1 +34).*cos(4*theta));
    % our derived 5th order equation
    eta5 = k^4*a^5/(192*(alpha1-1)^5*(alpha1+1)).*(121*alpha1^6 ...
    +432*alpha1^5+543*alpha1^4-1407*alpha1^3-258*alpha1^2 ...
    +2001*alpha1-1108).*cos(theta) ...
    +k^4*a^5/(128*(alpha1-1)^6*(3*alpha1+2)).*9*(57*alpha1^7 ...
    +204 *alpha1^6-53*alpha1^5-782*alpha1^4-741*alpha1^3-52*alpha1^2 ...
    +371*alpha1 +186).*cos(3*theta) ...
    +k^4*a^5/(384*(alpha1-1)^6*(12*alpha1^2+11*alpha1+2))*5.*(300*alpha1^8+1579*alpha1^7 ...
    +3176 *alpha1^6+2949*alpha1^5+1188*alpha1^4+675*alpha1^3+1326*alpha1^2 ...
    +827*alpha1 +130).*cos(5*theta);
    eta = eta1+eta2+eta3+eta4+eta5;


end
