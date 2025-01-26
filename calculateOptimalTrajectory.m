% 
% DISCLAIMER:
% This project is a recreation of the results in Blackmore, Acikmese, and Schaf's  research paper: “Minimum-Landing-Error Powered-Descent Guidance for  Mars Landing Using Convex Optimization."
% I claim no credit for any of the original ideas presented. I have only written the code as a project to learn more about their work. 

function [UVals, deltaT] = calculateOptimalTrajectory(In)
%% Setting Up the System

%Constants:
grav = [-3.7114, 0, 0]';
mDry = 1505;
rho1 = 4972;
rho2 = 13260;
alpha = 4.53e-4;
gammaAngle = 4*pi/180;

%Extracting the states from the input vector
r0 = In(1:3);
rDot0 = In(4:6);
mWet = In(7);
finalTime = In(8);
z0 = log(mWet);
Y0 = [r0; rDot0; z0];


%Setting up the guidance solution
N = 200;
guidanceSize = 200;
deltaT = finalTime/N;


%Continuous Model A and B matrices
A = [zeros(3) eye(3) zeros(3,1);...
         zeros(3) zeros(3) zeros(3,1);...
         zeros(1,3) zeros(1,3) 0];

Bgrav = [zeros(3); eye(3); zeros(1,3)];

Bu = [zeros(3) zeros(3,1);
          eye(3) zeros(3,1);
          zeros(1,3) -alpha];
B = [Bgrav, Bu];



%% Discretizing the System
syms t;
integralB = double(int(expm(A*(deltaT-t))*B,t,0,deltaT));

gamma = [0 0 0 1 0 0 0;
            0 0 0 0 1 0 0;
            0 0 0 0 0 1 0;
            0 0 0 0 0 0 1];

%Running through each index
for i=1:N
    %PHI Matrix
    tempA = expm(i*A*deltaT);
    Phi_K(i).mat = tempA;

    %Lambda and Psi matrices
    tempB = expm((i-1)*A*deltaT)*integralB;

    LambdaPsi = zeros(length(Y0), 7*N);
    if i==1
        LambdaPsi(:,1:length(tempB)) = tempB;
    else
        LambdaPsi(:,1:length(tempB)) = tempB;
        LambdaPsi(:,length(tempB)+1:end) = LambdaPsi_K(i-1).mat(:, 1:end-length(tempB));
    end
    LambdaPsi_K(i).mat = LambdaPsi;

    %Gamma matrices
    GamKIndeces = ((i-1)*size(gamma,2)+1):(i*size(gamma,2));
    tempGamma = zeros(size(gamma,1), 7*N);
    tempGamma(:, GamKIndeces) =  gamma;
    Gamma_k(i).mat = tempGamma;
end

%% Creating input matrix
clear Eta;
uVars = sdpvar(3,N, 'full');
sigmaVars = sdpvar(1,N, 'full');
U = [grav.*ones(3,N); uVars; sigmaVars];

for i=1:N
    indeces = ((i-1)*7+1):(i*7);
    Eta(indeces,1) = U(:,i);
end


%% Adding constraints

%Constants used for constraints
e1_7 = zeros(7,1); e1_7(1) = 1;
e1_3 = zeros(3,1); e1_3(1) = 1;
e4 = zeros(4,1); e4(4) = 1;
Z0_k = log(mWet-alpha*rho2*deltaT*(1:N));
S = [0 1 0; 0 0 1];
c = e1_3./tan(gammaAngle);

E = [eye(3), zeros(3,4)];
F = [zeros(1,6), 1];
Eu = [eye(3), zeros(3,1)];
Ev = [zeros(3), eye(3), zeros(3,1)];

%Creating the state vector variables
for i=1:N
    %State Vectors for all time. This needs to be donw first as next loop
    %references the "end" variable.
    Y_k(i).vec = Phi_K(i).mat*Y0 + LambdaPsi_K(i).mat*Eta;
end

%Placing the constraints
Constraints = [];
for i=1:N
    tempYk = Y_k(i).vec;

    %Constraint 1: Thrust through  all channels <= Thrust slack variable
    Constraints = [Constraints; norm(Eu*Gamma_k(i).mat*Eta,2)<= e4'*Gamma_k(i).mat*Eta];
    
    %Constraint 2: Thrust Slack Variable between some min and max values
    Constraints = [Constraints; rho1*exp(-Z0_k(i))*(1-(F*tempYk-Z0_k(i)) + ((F*tempYk-Z0_k(i))^2)/2) <= e4'*Gamma_k(i).mat*Eta];
    Constraints = [Constraints; e4'*Gamma_k(i).mat*Eta <= rho2*exp(-Z0_k(i))*(1-(F*tempYk-Z0_k(i)))]; 

    %Constraint 3: Descent within a cone
    Constraints = [Constraints; norm(S*( E*tempYk-E*Y_k(end).vec),2)-c'*(E*tempYk-E*Y_k(end).vec)<=0];

    fprintf("\tAdded %d/%d constraints\n", i, N);
end

%Constraint 4: final fuel usage not to exceed total fuel
Constraints = [Constraints; F*Y_k(end).vec>=log(mDry)];

%Constraint 5: final height = 0
Constraints = [Constraints; Y_k(end).vec'*e1_7 == 0];

%Constraint 6: Final Velocity = 0
Constraints = [Constraints; Ev*Y_k(end).vec == 0];

%Variable to optimize
optVar = norm((E*Y_k(end).vec),2);

%% Optimizing
ops = sdpsettings('solver','mosek','verbose',2,'debug',2);
optimize(Constraints,optVar, ops);


%% Extracting out the control variable
uMatrix = zeros(N+1,3);
massVec =  zeros(1, N+1);
massVec(1) = mWet;
ThrustMatrix = zeros(4, N+1);
for i = 1:N
    massVec(i+1) = exp(value(F*value(Y_k(i).vec)));
    
    tempVec = value(Gamma_k(i).mat*Eta);
    ThrustMatrix(:,i) = massVec(i).*tempVec(1:4);
end


timeOut = 0:deltaT:(guidanceSize-1)*deltaT;
ThrustX = zeros(1,guidanceSize);
ThrustY = zeros(1,guidanceSize);
ThrustZ = zeros(1,guidanceSize);
ThrustSlackVariable = zeros(1,guidanceSize);

if N > guidanceSize
    ThrustX = ThrustMatrix(1, 1:guidanceSize);
    ThrustY = ThrustMatrix(2, 1:guidanceSize);
    ThrustZ = ThrustMatrix(3, 1:guidanceSize);
    ThrustSlackVariable(1:N) = ThrustMatrix(4, 1:guidanceSize);
else
    ThrustX(1:N) = ThrustMatrix(1, 1:N);
    ThrustY(1:N) = ThrustMatrix(2, 1:N);
    ThrustZ(1:N) = ThrustMatrix(3, 1:N);
    ThrustSlackVariable(1:N) = ThrustMatrix(4, 1:N);
    
end

Out = [timeOut;
    ThrustX;
    ThrustY; 
    ThrustZ;
    ThrustSlackVariable
    ];


UVals = value(U);

end