% 
% DISCLAIMER:
% This project is a recreation of the results in Blackmore, Acikmese, and Schaf's  research paper: “Minimum-Landing-Error Powered-Descent Guidance for  Mars Landing Using Convex Optimization."
% I claim no credit for any of the original ideas presented. I have only written the code as a project to learn more about their work. 



%% Extracting data into vectors
fprintf("\nExtracting the Data\n");
PosMatrix = zeros(N,3);
VelMatrix = zeros(N,3);
FuelMatrix = zeros(N,1);

uMatrix = zeros(N+1,3);
minThrottle = zeros(N+1,1);
ThrustMatrix = zeros(N+1,1);
maxThrottle = zeros(N+1,1);

PosMatrix(1,:) = Y0(1:3);
VelMatrix(1,:) = Y0(4:6);
FuelMatrix(1,:) = Y0(7);


massVec = [mWet];
thrustMag = [];
capGamma = [];


for i = 1:N
    tempVec = value(Y_k(i).vec);
    PosMatrix(i+1,:) = tempVec(1:3);
    VelMatrix(i+1,:) = tempVec(4:6);
    FuelMatrix(i+1,:) = tempVec(7);

    minThrottle(i) = value(rho1*exp(-Z0_k(i))*(1-(F*tempVec-Z0_k(i)) + ((F*tempVec-Z0_k(i))^2)/2));
    maxThrottle(i) = value(rho2*exp(-Z0_k(i))*(1-(F*tempVec-Z0_k(i))));

    massVec = [massVec; exp(value(F*tempVec))];
    
    thrustMag = [thrustMag; massVec(i+1)*value(norm(Eu*Gamma_k(i).mat*Eta,2))];
    capGamma = [capGamma; massVec(i+1)*value(e4'*Gamma_k(i).mat*Eta)];
    
    tempVec = value(Gamma_k(i).mat*Eta);

    uMatrix(i,:) = tempVec(1:3);
    ThrustMatrix(i) = tempVec(4);
    UValues = value(U(:,i));


end

maxThrottle(end) = maxThrottle(end-1);
minThrottle(end) = minThrottle(end-1);

thrustMag = [thrustMag; 0];
capGamma = [capGamma; 0];



%% Plotting
fprintf("\nPlotting\n");
close all; 
%Fig1
figure; hold on; grid on;
title('Horizontal Plane Transfer')
plot(PosMatrix(:,2), PosMatrix(:,3))
xlim([-1500 2000]);
ylim([0 3500]);
xlabel('East [m]');
ylabel('North [m]');


%Fig2
figure; hold on; grid on;
title('Vertical Plane Transfer')
plot(PosMatrix(:,3), PosMatrix(:,1))
xlim([0 3500]);
ylim([-500 3000]);
gammaConstraintX = [tan(gammaAngle)*5000, 0, tan(gammaAngle)*5000];
gammaConstraintZ = PosMatrix(end,3)+[-5000, 0, 5000];
plot(gammaConstraintZ, gammaConstraintX, '--');
legend({"Position", "Descent Angle Constraint"})
xlabel('North [m]');
ylabel('Height [m]');

%Fig3
figure; hold on; grid on;
title('Angle Above Surface')
angleAboveSurface = 180/pi*atan2(PosMatrix(:,1)-PosMatrix(end,1), sqrt((PosMatrix(:,2)-PosMatrix(end,2)).^2+(PosMatrix(:,3)-PosMatrix(end,3)).^2));
plot(time, angleAboveSurface);
plot(time, 0*time+4)
legend({"Angle", "Limit"})


%fig4
thrust100 = rho2/8*10;
figure; hold on; grid on;
title('Throttle Level')
plot(time,capGamma/thrust100)
plot(time,rho1/thrust100*ones(length(time),1));
plot(time,rho2/thrust100*ones(length(time),1));
legend({"Throttle", "Lower Limit" ,"Upper Limit"});
ylim([0 1])

%Fig5
figure; hold on; grid on;
title("Position Vs time")
plot(time,PosMatrix(:,1))
plot(time,PosMatrix(:,2))
plot(time,PosMatrix(:,3))
legend({"Height [m]", "East [m]" ,"North [m]"});
xlabel('Time [sec]');
ylabel('Distance from Target [m]');


%Fig6
figure; hold on; grid on;
title("Velocity Vs time")
plot(time,VelMatrix(:,1))
plot(time,VelMatrix(:,2))
plot(time,VelMatrix(:,3))
legend({"Vertical Vel [m/s]", "East Vel [m/s]" ,"North Vel [m/s]"});
xlabel('Time [sec]');
ylabel('Velocity Relative to Target [m]');

%Fig7
figure; hold on; grid on;
title("Controlled Acceleration Vs time")
plot(time,uMatrix(:,1))
plot(time,uMatrix(:,2))
plot(time,uMatrix(:,3))
legend({"X Accel", "Y Accel" ,"Z Accel"});

%Fig8
thrustAngle = 180/pi*atan2(sqrt(uMatrix(:,2).^2+uMatrix(:,3).^2), uMatrix(:,1));
rate = diff(thrustAngle)'./diff(time);

figure; hold on; grid on;
title("Thrust Angle and Rate Vs time")
plot(time, thrustAngle);
plot(time(1:end-1), rate);
legend({"Angle", "Rate"});


%% Continuous Time
deltaT_c = 0.1;
Time_c = 0:deltaT_c:finalTime;
indexRatio = deltaT/deltaT_c;

U_c = [];
for i = 1:N
 %   indecesLower = indexRatio*(i-1)+1;
 %   indecesUpper = indexRatio*i;
    U_c = [U_c, value(U(:,i)).*ones(7,indexRatio)];
end

U_c = [U_c, zeros(7,1)];


C = eye(7);
D = 0;
ContinuousModel = ss(A,B,C,D);

  
[XOut] =  lsim(ContinuousModel, U_c, Time_c, Y0);

%%
hold on
plot(Time_c, XOut(:,1), '--')
plot(Time_c, XOut(:,2), '--')
plot(Time_c, XOut(:,3), '--')
