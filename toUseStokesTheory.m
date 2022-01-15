% % Stokes theory, 
% % 1. input region
% clc,clear

h0 = 10; 
% For the first mode, a is first harmonic first order amplitude
modeNo = 1; 
a = 1.6; % either a or H is given
% H = 3.2; % either a or H is given
T = 5;
theta0 = 0:0.01:2*pi;
%
Result =StokesDispSolver('h', h0, 'T', T, 'a', a,'mode', modeNo)
if  modeNo ==1
    a = Result.a;
    H = Result.H;
elseif modeNo ==2
    a = Result.aw;
    H = Result.Hw;
end

% Below is to call the free surface function StokesEta, make your own edit to select the preferred output. 
[eta, eta1, eta2, eta3,eta4,eta5] = StokesEta(Result.k, h, Result.a, theta); % free surface
plot(theta0,eta)
c = Result.L/T; % phase speed

% velocity and potential function
for i = 1:length(theta)    
% z = [-10:0.5:-1.5, eta(i)];%linspace(-h, eta(i), 10);
z = [linspace(-h, eta(i), 10), -2.5];% 
[phi(i,:), u, u1(i,:),u2(i,:),u3(i,:),u4(i,:),u5(i,:), ...
    w,w1(i,:),w2(i,:),w3(i,:),w4(i,:),w5(i,:)] = StokesU(Result.k, h, Result.a, theta(i), z);
    
end

Result.u = u;
Result.w = w;

% etaFen = FentonEta(Result.k, h0, a, theta0)

% % For the second mode, a is combined first harmonic amplitude
% modeNo = 2; 
% Result2 =StokesDispSolver('h', h0, 'T', T0, 'a', a,'mode', 2);

% Then we can use StokesEta and StokesU functions



