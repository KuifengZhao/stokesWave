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
% k1 = Result1.k; L1 = Result1.L; % wave number and wave length
eta = StokesEta(Result.k, h0, a, modeNo, theta0); % free surface
plot(theta0,eta)
c = Result.L/T; % phase speed

% velocity and potential function
for i = 1:length(theta0)    
z = linspace(-h0, eta(i), 10);
[phi(:,i), u(:,i), w(:,i)] = StokesU(Result.k, h0, a, modeNo, theta0(i), z);
% plot(u,z,'ro', w,z, 'kx'),hold on,
% axis([-4 4 -1 0.4])
% hwl = refline(0, eta(i)); hwl.Color = 'b'; hwl.LineWidth = 2;
% pause(0.02)
% hold off,
end
ResultCase3 = Result;
[a, H, max(eta), min(eta), max(eta)- min(eta), abs(max(eta)- min(eta) - H)/H]
etaCase3= eta;
uCase3 = u;
wCase3 = w;

% etaFen = FentonEta(Result.k, h0, a, theta0)

% % For the second mode, a is combined first harmonic amplitude
% modeNo = 2; 
% Result2 =StokesDispSolver('h', h0, 'T', T0, 'a', a,'mode', 2);
% k2 = Result2.k; L2 = Result2.L;
% eta = StokesEta(k2, h0, a, 2, theta0);
% for i = 1:length(theta0)    
% z = linspace(-h0, eta(i), 10);
% [phi(:,i), u(:,i), w(:,i)] = StokesU(k2, h0, a, modeNo, theta0(i), z);
% end
% % 2. function region


