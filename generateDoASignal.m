%% My signal generation function

function [R,M,antLoc,lambda,X] = generateDoASignal(doas,nAnt,nSymbols,lenCarrier,noisePower)

dt = 1e-12;                                 % Sample period
fs = 1/dt;                                  % Sample frequency
T = lenCarrier;                             % Length of carrier signal
t = (0:dt:(T-1)*dt)';                       % Time vector
fc = fs/64;                                 % Carrier Frequency
lambda = physconst('lightspeed')/fc;        % Wavelength of carrier frequency

% Antenna Parameters
L = nAnt;                                   % Number of antennas
dist = lambda/2;                            % Distance between antennas
antLoc = (0:dist:(L-1)*dist).';                  % Distances among antennas
theta_m = doas;                             % Direction of arrival
w_m = 2*pi*antLoc/lambda*sin(theta_m*pi/180);    % Phase differences between antennas
M = length(theta_m);                        % Number of sources
A = exp(1j*w_m);                            % Steering Vectors

% figure; scatter(antLoc/lambda,zeros(size(antLoc)),'x'); xlabel("\lambda"); title("Linear Uniform Array")

% Signal Parameters
N = nSymbols;                               % Length of message
s = zeros(M,length(t)*N);
m = zeros(M,N);
for iD = 1:M
    m(iD,:) = exp(1j*pi/2*randi([1,4],[1,N])+1j*pi/4); % Generate QPSK weights
    stemp = m(iD,:).*exp(2j*pi*fc*t);
    s(iD,:) = stemp(:); clear stemp
end

% Declare Noise n and normalize to unit power
n = normrnd(0,1,[L,length(s)]) + 1j*normrnd(0,1,[L,length(s)]);
n = n./sqrt(sum(abs(n).^2,2)/nSymbols/noisePower);

n = wgn(L,length(s),noisePower);

% Declare Xtemp and normalize to unit power
Xtemp = A*s;
Xtemp = Xtemp./sqrt(sum(abs(Xtemp).^2,2)/nSymbols);

X = Xtemp + n;
R = 1/(length(t)*N)*(X*X');

end