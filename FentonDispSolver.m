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

function [Results] = FentonDispSolver(varargin)
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
%         case 'omega0'
%             omega0 = values{kIn};
        case 'mode'
            modeNo = values{kIn};
    end
end
Results.h = h;
Results.T = T;
omega = 2*pi/T;
Results.omega = omega;
syms k
g = 9.81; 
S = 1/cosh(2*k*h);
alpha = 1/S;
if and(exist('a', 'var') ==1, exist('T', 'var') ==1)    
    omega2Fen = (2+7*S^2)/(4*(1-S)^2);
    switch modeNo
        case 1
            omega4Fen = (4+32*S - 116*S^2-400*S^3-71*S^4+146*S^5)/(32*(1-S)^5);
        case 2
            omega4Fen = (16*alpha^5+68*alpha^4-38*alpha^3-250*alpha^2 ...
                +55*alpha+230)/32/(alpha-1)^5;
    end
    kResults = double(vpasolve(omega - sqrt(g*k*tanh (k*h))*(1+omega2Fen*(k*a)^2+omega4Fen*(k*a)^4)==0, [0 Inf]));
    
    switch modeNo
        case 1
            Results.a = a;
            Results.H = 2*a;
        case 2
            aw =a;
            S = 1/cosh(2*kResults*h); alpha = 1/S;
            B33 =  3*(1+3*S+3*S^2+2*S^3)/8/(1-S)^3;
            B53 = -(-33*alpha^7 + 4*alpha^6 +553*alpha^5 ...
            +1336*alpha^4 +1239*alpha^3 +362*alpha^2 ...
            -139*alpha -82)*9/((3*alpha + 2)*(alpha - 1)^6*128);
            B55 =5*(300+1579*S+3176*S.^2+2949*S.^3+1188*S.^4+675*S.^5+ ...
                1326*S.^6+827*S.^7+130*S.^8)./(384*(3+2*S).*(4+S).*(1-S).^6);
            Results.aw = a;
            Results.Hw = (2*aw+2*B33*kResults^2*aw^3+2*(B53+B55)*kResults^4*aw^5);
            
    end
%     kResults = double(vpasolve(omegaFun, k, [0 Inf]));
    Results.omega0 = sqrt(g*kResults*tanh(kResults*h));    
    
    
elseif and(exist('H', 'var') ==1, exist('T', 'var') ==1)

    omega0 = sqrt(g*k*tanh(k*h));   
    B33 =  3*(1+3*S+3*S^2+2*S^3)/8/(1-S)^3;
    B55   =   5*(300+1579*S+3176*S.^2+2949*S.^3+1188*S.^4+675*S.^5+1326*S.^6+827*S.^7+130*S.^8)./(384*(3+2*S).*(4+S).*(1-S).^6);
    
    switch modeNo
        case 1
            Results.a = H/2;
            Results.H = H;            
            a = H/2;
            omega2 = (2+7*S^2)/(4*(1-S)^2); %(2*alpha^2+7)/4/(alpha -1)^2;
            omega4Fen = (4+32*S - 116*S^2-400*S^3-71*S^4+146*S^5)/(32*(1-S)^5);
            omegaFun1 = omega0*(1+k^2*a^2*omega2+k^4*a^4*omega4Fen) - omega ==0;            
            kResults = double(vpasolve(omegaFun1, k, [0 Inf]));
%             kResults = abs(double(S));
        case 2    
            syms aw
            B53 = 1/(128*(alpha-1)^6*(3*alpha+2)).*9*(51*alpha^7 ...
            +116 *alpha^6-211*alpha^5-760*alpha^4-597*alpha^3+106*alpha^2 ...
            +355*alpha +130);
%             B53 = -(-33*cosh(2*h*k)^7 + 4*cosh(2*h*k)^6 +553*cosh(2*h*k)^5 ...
%             +1336*cosh(2*h*k)^4 +1239*cosh(2*h*k)^3 +362*cosh(2*h*k)^2 ...
%             -139*cosh(2*h*k) -82)*9/((3*cosh(2*h*k) + 2)*(cosh(2*h*k) - 1)^6*128);        
            Hw = H;            
            Hfun2 = Hw - (2*aw+2*B33*k^2*aw^3+2*(B53+B55)*k^4*aw^5)==0;
            omega2 = (2*alpha^2+7)/4/(alpha -1)^2;
            omega4 = (16*alpha^5+68*alpha^4-38*alpha^3-250*alpha^2 ...
                +55*alpha+230)/32/(alpha-1)^5;
            omegaFun2 = omega0*(1+k^2*aw^2*omega2+k^4*aw^4*omega4) - omega ==0;      
            SHcase2 = vpasolve([Hfun2 ==0, omegaFun2==0], [k, aw]);
            kResults = abs(double(SHcase2.k));
            aResults = double(SHcase2.aw);
            Results.aw = aResults;
            Results.Hw = Hw;
    end
    Results.omega0 = sqrt(g*kResults*tanh(kResults*h));


else    
end

Results.k = kResults;
Results.L = 2*pi/kResults;

end


