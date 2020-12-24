n=50;
sigmamin = .1;
sigmamax = .5;

Mu1 = ones(n,1)*[1 -1];
IsoNoise1 = randn(n,2);
RotMat = (sqrt(2)/2)*[1 -1; 1 1];
Noise1 = IsoNoise1*diag([sigmamin sigmamax])*RotMat;
X1 = Mu1 + Noise1;

Mu2 = ones(n,1)*[1 1];
IsoNoise2 = randn(n,2);
Noise2 = IsoNoise2*diag([sigmamin sigmamax])*RotMat;
X2 = Mu2 + Noise2;

Mu3 = ones(n,1)*[-1 1];
IsoNoise3 = randn(n,2);
Noise3 = IsoNoise3*diag([sigmamin sigmamax])*RotMat;
X3 = Mu3 + Noise3;

Mu4 = ones(n,1)*[-1 -1];
IsoNoise4 = randn(n,2);
Noise4 = IsoNoise4*diag([sigmamin sigmamax])*RotMat;
X4 = Mu4 + Noise4;

X=[X1; X2; X3; X4];

scatter(X(:,1),X(:,2))
xlim([-2 2])
ylim([-2 2])

%%
SqCov = diag([sigmamin sigmamax])*RotMat;
InvSqCov = inv(SqCov);

NewX = X*InvSqCov;

figure
scatter(NewX(:,1),NewX(:,2))
xlim([-16 16])
ylim([-16 16])