%% Generate Signal

doas = [-40,40,60];     % Directions of arrival in degrees
nAnt = 128;             % Number of antenna elements
nSymbols = 50;          % Number of symbols
lenCarrier = 1;         % Length of the carrier sequence
noisePower = 10;        % Noise power

[R,nDoa,antLoc,lambda] = generateDoASignalNonuniform(doas,nAnt,nSymbols,lenCarrier,noisePower);

% MVDR DoA Function

Pmvdr = mvdrDoA(R,antLoc,lambda);

% MVDR-SVM DoA Function

Psvm_mvdr = svm_mvdrDoA(R,antLoc,lambda);

% MUSIC DoA Function

Pmusic = musicDoA(R,nDoa,antLoc,lambda);

% MUSIC-SVM DoA Function

Psvm_music = svm_musicDoA(R,nDoa,antLoc,lambda);

% Compare Classical Algorithms with SVM Algorithms

plotCompare(Pmvdr,Psvm_mvdr,Pmusic,Psvm_music,-20);

%% Plot All

plotAll(Pmvdr,Psvm_mvdr,Pmusic,Psvm_music,-20);
