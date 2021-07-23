% This is a solver for nonlinear dispersion relationship of Stokes theory.
% The full nonlinear dispersion relationship has h (mean water depth), either T
% (wave period) or omega0 (leading order frequency), and 
% either a or H is given. Results include k, H, a, omega, omega0. 
% Mode 1, the first order wave amplitude/height is given as a, H.
% Mode 2, the first harmonic term's coefficient is given as aw, Hw. 
% Example 1, 
% Results =StokesDispSolver('h',1,'T', 5, 'a',0.05, 'mode',1)
% Example 2, 
% Results =StokesDispSolver('h',1,'T', 5, 'a',0.05, 'mode',2)
% Example 3, 
% Results =StokesDispSolver('h',1,'T', 5, 'H',0.1, 'mode',1)
% Example 4, 
% Results =StokesDispSolver('h',1,'T', 5, 'H',0.1, 'mode',2)
% Example 5, 
% Results =StokesDispSolver('h',1,'omega0', 1.2566, 'a',0.15, 'mode',1)
% Example 6, 
% Results =StokesDispSolver('h',1,'omega0', 1.2566, 'H',0.3, 'mode',1)

function [Results] = StokesDispSolver(varargin)
names = varargin(1:2:end);
values = varargin(2:2:end);
for kIn = 1:numel(names)
    switch names{kIn}
        case 'a'
            a = values{kIn};
        case 'H'
            H = values{kIn};
        case 'h'
            h = values{kIn};
        case 'T'
            T = values{kIn};
        case 'omega0'
            omega0 = values{kIn};
        case 'mode'
            modeNo = values{kIn};
    end
end
Results.h = h;
syms k
g = 9.81; 
alpha1 = cosh(2*k*h);

if and(exist('a', 'var') ==1, exist('T', 'var') ==1)
    omega0 = sqrt(g*k*tanh(k*h));
Results.T = T;
omega = 2*pi/T;
Results.omega = omega;
% correction to omega, L, k with omega2 and omega4
    switch modeNo
        case 1 % in case first order wave amplitude a is given
            Results.a = a;
            omega2 = k^2*a^2*(2*alpha1^2+7)/4/(alpha1 -1)^2;
            omega4 = (a^4*k^4*(20*alpha1^6 + 132*alpha1^5 + 20*alpha1^4 ...
                - 184*alpha1^3 - 243*alpha1^2 - 11*alpha1 + 428))/(32*(alpha1 - 1)^5*(alpha1 + 1));
            omegaFun = omega0*(1+omega2+omega4) - omega;
        case 2 % in case first harmonic wave amplitude aw is given
            aw = a;
            Results.aw = a;
            omega2 = k^2*aw^2*(2*alpha1^2+7)/4/(alpha1 -1)^2;
            omega4 = (aw^4*k^4*(18*alpha1^6 + 108*alpha1^5 + 29*alpha1^4 ...
                - 244*alpha1^3 - 201*alpha1^2 + 73*alpha1 + 379))/(32*(alpha1 - 1)^5*(alpha1 + 1));
            omegaFun = omega0*(1+omega2+omega4) - omega;        
    end
    kResults = double(vpasolve(omegaFun, k, [0 Inf]));
    Results.omega0 = sqrt(g*kResults*tanh(kResults*h));
    
        B33 = (cosh(2*h*kResults)^3 + 3*cosh(2*h*kResults)^2 + 3*cosh(2*h*kResults) + 2)/(8*(cosh(2*h*kResults) - 1)^3);  
        B55 = (300*cosh(2*h*kResults)^8 + 1579*cosh(2*h*kResults)^7 + 3176*cosh(2*h*kResults)^6 + 2949*cosh(2*h*kResults)^5 ...
            + 1188*cosh(2*h*kResults)^4 + 675*cosh(2*h*kResults)^3 + 1326*cosh(2*h*kResults)^2 + 827*cosh(2*h*kResults) ...
            + 130)*5/((cosh(2*h*kResults) - 1)^6*(12*cosh(2*h*kResults)^2 + 11*cosh(2*h*kResults) + 2)*384);
%     [B31,B33,B51,B53,B55] = findBcoeff(h,kResults);
    switch modeNo
        case 1
            B31 = (cosh(2*h*kResults)^2 + 12*cosh(2*h*kResults) - 7)/(8*(cosh(2*h*kResults) - 1)^2);
            B51 = (121*cosh(2*h*kResults)^6 + 432*cosh(2*h*kResults)^5 + 543*cosh(2*h*kResults)^4 - 1407*cosh(2*h*kResults)^3 ...
            - 258*cosh(2*h*kResults)^2 + 2001*cosh(2*h*kResults) - 1108)/(192*(cosh(2*h*kResults) - 1)^5*(cosh(2*h*kResults) + 1));
            B53 = (57*cosh(2*h*kResults)^7 + 204*cosh(2*h*kResults)^6 - 53*cosh(2*h*kResults)^5 ...
            -782*cosh(2*h*kResults)^4 - 741*cosh(2*h*kResults)^3 - 52*cosh(2*h*kResults)^2 ...
            + 371*cosh(2*h*kResults) + 186)*9/((3*cosh(2*h*kResults) + 2)*(cosh(2*h*kResults) - 1)^6*128);
            Hout = 2*a+2*(B31+B33)*kResults^2*a^3+2*(B51+B53+B55)*kResults^4*a^5;
            Results.H = Hout;
        case 2
            B53 = (51*cosh(2*h*kResults)^7 + 116*cosh(2*h*kResults)^6 - 211*cosh(2*h*kResults)^5 ...
            -760*cosh(2*h*kResults)^4 - 597*cosh(2*h*kResults)^3 +106*cosh(2*h*kResults)^2 ...
            + 355*cosh(2*h*kResults) + 130)*9/((3*cosh(2*h*kResults) + 2)*(cosh(2*h*kResults) - 1)^6*128);
            Hout = 2*aw+2*B33*kResults^2*aw^3+2*(B53+B55)*kResults^4*aw^5;
            Results.Hw = Hout;
    end
    
elseif and(exist('H', 'var') ==1, exist('T', 'var') ==1)
    omega0 = sqrt(g*k*tanh(k*h));
    Results.T = T;
    omega = 2*pi/T;
    Results.omega = omega;
    syms a aw
     
        B33 = (cosh(2*h*k)^3 + 3*cosh(2*h*k)^2 + 3*cosh(2*h*k) + 2)/(8*(cosh(2*h*k) - 1)^3);        
        B55 = (300*cosh(2*h*k)^8 + 1579*cosh(2*h*k)^7 + 3176*cosh(2*h*k)^6 + 2949*cosh(2*h*k)^5 ...
            + 1188*cosh(2*h*k)^4 + 675*cosh(2*h*k)^3 + 1326*cosh(2*h*k)^2 + 827*cosh(2*h*k) ...
            + 130)*5/((cosh(2*h*k) - 1)^6*(12*cosh(2*h*k)^2 + 11*cosh(2*h*k) + 2)*384);
%         [B31,B33,B51,B53,B55] = findBcoeff(h,k);
    switch modeNo
        case 1
            B31 = (cosh(2*h*k)^2 + 12*cosh(2*h*k) - 7)/(8*(cosh(2*h*k) - 1)^2);
            B51 = (121*cosh(2*h*k)^6 + 432*cosh(2*h*k)^5 + 543*cosh(2*h*k)^4 - 1407*cosh(2*h*k)^3 ...
            - 258*cosh(2*h*k)^2 + 2001*cosh(2*h*k) - 1108)/(192*(cosh(2*h*k) - 1)^5*(cosh(2*h*k) + 1));
            B53 = (57*cosh(2*h*k)^7 + 204*cosh(2*h*k)^6 - 53*cosh(2*h*k)^5 ...
            -782*cosh(2*h*k)^4 - 741*cosh(2*h*k)^3 - 52*cosh(2*h*k)^2 ...
            + 371*cosh(2*h*k) + 186)*9/((3*cosh(2*h*k) + 2)*(cosh(2*h*k) - 1)^6*128);
            Results.H = H;
            Hfun1 = H -(2*a+2*(B31+B33)*k^2*a^3+2*(B51+B53+B55)*k^4*a^5) ==0;
            omega2 = k^2*a^2*(2*alpha1^2+7)/4/(alpha1 -1)^2;
            omega4 = (a^4*k^4*(20*alpha1^6 + 132*alpha1^5 + 20*alpha1^4 ...
                - 184*alpha1^3 - 243*alpha1^2 - 11*alpha1 + 428))/(32*(alpha1 - 1)^5*(alpha1 + 1));
            omegaFun1 = omega0*(1+omega2+omega4) - omega ==0;
            S = vpasolve([Hfun1 , omegaFun1], [k, a]);
            kResults = abs(double(S.k));
            aResults = double(S.a);
            Results.a = aResults;
        case 2
            B53 = (51*cosh(2*h*k)^7 + 116*cosh(2*h*k)^6 - 211*cosh(2*h*k)^5 ...
            -760*cosh(2*h*k)^4 - 597*cosh(2*h*k)^3 +106*cosh(2*h*k)^2 ...
            + 355*cosh(2*h*k) + 130)*9/((3*cosh(2*h*k) + 2)*(cosh(2*h*k) - 1)^6*128);
            Hw = H;
            Results.Hw = Hw;
            Hfun2 = Hw - (2*aw+2*B33*k^2*aw^3+2*(B53+B55)*k^4*aw^5)==0;
            omega2 = k^2*aw^2*(2*alpha1^2+7)/4/(alpha1 -1)^2;
            omega4 = (aw^4*k^4*(18*alpha1^6 + 108*alpha1^5 + 29*alpha1^4 ...
                - 244*alpha1^3 - 201*alpha1^2 + 73*alpha1 + 379))/(32*(alpha1 - 1)^5*(alpha1 + 1));
            omegaFun2 = omega0*(1+omega2+omega4) - omega ==0;      
            S = vpasolve([Hfun2 ==0, omegaFun2==0], [k, aw]);
            kResults = abs(double(S.k));
            aResults = double(S.aw);
            Results.aw = aResults;
    end
    Results.omega0 = sqrt(g*kResults*tanh(kResults*h));
elseif  and(exist('a', 'var') ==1, exist('omega0', 'var') ==1)
    Results.omega0 = omega0;
    kResults = double(vpasolve(omega0^2 - g*k*tanh(k*h), k, [0 inf]));
    alpha1 = cosh(2*kResults*h);    
    B33 = (cosh(2*h*kResults)^3 + 3*cosh(2*h*kResults)^2 + 3*cosh(2*h*kResults) + 2)/(8*(cosh(2*h*kResults) - 1)^3);    
    B55 = (300*cosh(2*h*kResults)^8 + 1579*cosh(2*h*kResults)^7 + 3176*cosh(2*h*kResults)^6 + 2949*cosh(2*h*kResults)^5 ...
        + 1188*cosh(2*h*kResults)^4 + 675*cosh(2*h*kResults)^3 + 1326*cosh(2*h*kResults)^2 + 827*cosh(2*h*kResults) ...
        + 130)*5/((cosh(2*h*kResults) - 1)^6*(12*cosh(2*h*kResults)^2 + 11*cosh(2*h*kResults) + 2)*384);
    switch modeNo
        case 1
            B31 = (cosh(2*h*kResults)^2 + 12*cosh(2*h*kResults) - 7)/(8*(cosh(2*h*kResults) - 1)^2);
            B51 = (121*cosh(2*h*kResults)^6 + 432*cosh(2*h*kResults)^5 + 543*cosh(2*h*kResults)^4 - 1407*cosh(2*h*kResults)^3 ...
                - 258*cosh(2*h*kResults)^2 + 2001*cosh(2*h*kResults) - 1108)/(192*(cosh(2*h*kResults) - 1)^5*(cosh(2*h*kResults) + 1));
            B53 = (57*cosh(2*h*kResults)^7 + 204*cosh(2*h*kResults)^6 - 53*cosh(2*h*kResults)^5 ...
                -782*cosh(2*h*kResults)^4 - 741*cosh(2*h*kResults)^3 - 52*cosh(2*h*kResults)^2 ...
                + 371*cosh(2*h*kResults) + 186)*9/((3*cosh(2*h*kResults) + 2)*(cosh(2*h*kResults) - 1)^6*128);
            Results.a = a;
            Hout = 2*a+2*(B31+B33)*kResults^2*a^3+2*(B51+B53+B55)*kResults^4*a^5;
            Results.H = Hout;
            omega2 = kResults^2*a^2*(2*alpha1^2+7)/4/(alpha1 -1)^2;
            omega4 = (a^4*kResults^4*(20*alpha1^6 + 132*alpha1^5 + 20*alpha1^4 ...
                - 184*alpha1^3 - 243*alpha1^2 - 11*alpha1 + 428))/(32*(alpha1 - 1)^5*(alpha1 + 1));
            omega = omega0*(1+omega2+omega4);
            Results.omega = double(omega);
            Results.T = 2*pi/omega;
        case 2
            B53 = (51*cosh(2*h*kResults)^7 + 116*cosh(2*h*kResults)^6 - 211*cosh(2*h*kResults)^5 ...
            -760*cosh(2*h*kResults)^4 - 597*cosh(2*h*kResults)^3 +106*cosh(2*h*kResults)^2 ...
            + 355*cosh(2*h*kResults) + 130)*9/((3*cosh(2*h*kResults) + 2)*(cosh(2*h*kResults) - 1)^6*128);
            Hout = 2*aw+2*B33*kResults^2*aw^3+2*(B53+B55)*kResults^4*aw^5;
            Results.Hw = Hout;
            aw = a;
            Results.aw = a;
            omega2 = kResults^2*aw^2*(2*alpha1^2+7)/4/(alpha1 -1)^2;
            omega4 = (aw^4*kResults^4*(18*alpha1^6 + 108*alpha1^5 + 29*alpha1^4 ...
                - 244*alpha1^3 - 201*alpha1^2 + 73*alpha1 + 379))/(32*(alpha1 - 1)^5*(alpha1 + 1));
            omega = omega0*(1+omega2+omega4) ;
            Results.omega = omega;
            Results.T = 2*pi/omega;
    end
elseif and(exist('H', 'var') ==1, exist('omega0', 'var') ==1)
    Results.omega0 = omega0;
    kResults = double(vpasolve(omega0^2 - g*k*tanh(k*h), k, [0 inf]));
    alpha1 = cosh(2*kResults*h);
    syms a aw    
    B33 = (cosh(2*h*kResults)^3 + 3*cosh(2*h*kResults)^2 + 3*cosh(2*h*kResults) + 2)/(8*(cosh(2*h*kResults) - 1)^3);    
    B55 = (300*cosh(2*h*kResults)^8 + 1579*cosh(2*h*kResults)^7 + 3176*cosh(2*h*kResults)^6 + 2949*cosh(2*h*kResults)^5 ...
        + 1188*cosh(2*h*kResults)^4 + 675*cosh(2*h*kResults)^3 + 1326*cosh(2*h*kResults)^2 + 827*cosh(2*h*kResults) ...
        + 130)*5/((cosh(2*h*kResults) - 1)^6*(12*cosh(2*h*kResults)^2 + 11*cosh(2*h*kResults) + 2)*384);
    
%         [B31,B33,B51,B53,B55] = findBcoeff(h,k);
    switch modeNo
        case 1
            B31 = (cosh(2*h*kResults)^2 + 12*cosh(2*h*kResults) - 7)/(8*(cosh(2*h*kResults) - 1)^2);
            B51 = (121*cosh(2*h*kResults)^6 + 432*cosh(2*h*kResults)^5 + 543*cosh(2*h*kResults)^4 - 1407*cosh(2*h*kResults)^3 ...
                - 258*cosh(2*h*kResults)^2 + 2001*cosh(2*h*kResults) - 1108)/(192*(cosh(2*h*kResults) - 1)^5*(cosh(2*h*kResults) + 1));
            B53 = (57*cosh(2*h*kResults)^7 + 204*cosh(2*h*kResults)^6 - 53*cosh(2*h*kResults)^5 ...
                -782*cosh(2*h*kResults)^4 - 741*cosh(2*h*kResults)^3 - 52*cosh(2*h*kResults)^2 ...
                + 371*cosh(2*h*kResults) + 186)*9/((3*cosh(2*h*kResults) + 2)*(cosh(2*h*kResults) - 1)^6*128);              
            Results.H = H;
            Hfun1 = H -(2*a+2*(B31+B33)*kResults^2*a^3+2*(B51+B53+B55)*kResults^4*a^5) ==0;
            asol = double(vpasolve(Hfun1, [0 H]));
            omega2 = kResults^2*asol^2*(2*alpha1^2+7)/4/(alpha1 -1)^2;
            omega4 = (asol^4*kResults^4*(20*alpha1^6 + 132*alpha1^5 + 20*alpha1^4 ...
                - 184*alpha1^3 - 243*alpha1^2 - 11*alpha1 + 428))/(32*(alpha1 - 1)^5*(alpha1 + 1));
            omega = omega0*(1+omega2+omega4);
            Results.a = asol;
            Results.omega = omega;
            Results.T = 2*pi/omega;
        case 2
            B53 = (51*cosh(2*h*kResults)^7 + 116*cosh(2*h*kResults)^6 - 211*cosh(2*h*kResults)^5 ...
            -760*cosh(2*h*kResults)^4 - 597*cosh(2*h*kResults)^3 +106*cosh(2*h*kResults)^2 ...
            + 355*cosh(2*h*kResults) + 130)*9/((3*cosh(2*h*kResults) + 2)*(cosh(2*h*kResults) - 1)^6*128);
            Hw = H;
            Results.Hw = Hw;
            Hfun2 = Hw - (2*aw+2*B33*kResults^2*aw^3+2*(B53+B55)*kResults^4*aw^5)==0;
            awsol = double( vpasolve(Hfun2, [0 Hw]));
            omega2 = kResults^2*awsol^2*(2*alpha1^2+7)/4/(alpha1 -1)^2;
            omega4 = (awsol^4*kResults^4*(18*alpha1^6 + 108*alpha1^5 + 29*alpha1^4 ...
                - 244*alpha1^3 - 201*alpha1^2 + 73*alpha1 + 379))/(32*(alpha1 - 1)^5*(alpha1 + 1));
            omega = omega0*(1+omega2+omega4);  
            Results.aw = awsol;
            Results.omega = omega;
            Results.T = 2*pi/omega;
    end          
else    
end

Results.k = kResults;
Results.L = 2*pi/kResults;

end


