function [etaFen,etaFenH] =FentonEta(k, h, a, modeNo, theta)

    S = 1/cosh(2*k*h);
    alpha1 = 1/S;
    B22Fen   =   coth(k*h)*(1+2*S)/(2*(1-S));
    B33Fen =  3*(1+3*S+3*S^2+2*S^3)/8/(1-S)^3;
    B44Fen   =   coth(k*h) *(24+92*S+122*S^2+66*S^3+67*S^4+34*S^5)/(24*(3+2*S)*(1-S)^4);
    B55Fen   =   5*(300+1579*S+3176*S.^2+2949*S.^3+1188*S.^4+675*S.^5+1326*S.^6+827*S.^7+130*S.^8)./(384*(3+2*S).*(4+S).*(1-S).^6);
    
    switch modeNo
        case 1    
    B31Fen   =   -3*(1+3*S+3*S^2+2*S^3)/8/(1-S)^3;
    B42Fen   =   coth(k* h)*(6-26*S-182*S^2-204*S^3-25*S^4+26*S^5)/(6*(3+2*S)*(1-S)^4);
    B53Fen   =   9*(132+17*S-2216*S.^2-5897*S.^3-6292*S.^4-2687*S.^5+194*S.^6 +467*S.^7+82*S.^8)./(128*(3+2*S).*(4+S).*(1-S).^6);
    B51Fen = -(B53Fen+B55Fen);
        case 2
            B31Fen   =0;
            B42Fen   =1/(24*(alpha1-1)^4*tanh(k*h)).*(17*alpha1^4+5*alpha1^3-135*alpha1^2 ...
    -56*alpha1+88);
            B51Fen = 0;
            B53Fen = 1/(128*(alpha1-1)^6*(3*alpha1+2)).*9*(51*alpha1^7 ...
    +116 *alpha1^6-211*alpha1^5-760*alpha1^4-597*alpha1^3+106*alpha1^2 ...
    +355*alpha1 +130);
    end
    
    
    etaFen = a*cos(theta) + k*a^2*B22Fen*cos(2*theta) + ...
        k^2*a^3*(B31Fen*cos(theta)+B33Fen*cos(3*theta))+ ...
        k^3*a^4*(B42Fen*cos(2*theta) + B44Fen*cos(4*theta))+...
        k^4*a^5*(B51Fen*cos(theta)+B53Fen*cos(3*theta)+B55Fen*cos(5*theta));
    etaH1 = (a + k^2*a^3* B31Fen+ k^4*a^5*B51Fen)*cos(theta);
    etaH2 = (k*a^2*B22Fen + k^3*a^4*B42Fen) *cos(2*theta);
    etaH3 = (k^2*a^3*B33Fen + k^4*a^5*B53Fen)*cos(3*theta);
    etaH4 = k^3*a^4*B44Fen*cos(4*theta);
    etaH5 = k^4*a^5*B55Fen*cos(5*theta);
% each harmonic contribution collected in one variable
    etaFenH = [etaH1', etaH2', etaH3', etaH4', etaH5'];
    
    
    
