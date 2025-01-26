%
% DISCLAIMER:
% This project is a recreation of the results in Blackmore, Acikmese, and Schaf's  research paper: “Minimum-Landing-Error Powered-Descent Guidance for  Mars Landing Using Convex Optimization."
% I claim no credit for any of the original ideas presented. I have only written the code as a project to learn more about their work. 

%% Add LMI solver paths
pathToYALMIP = '';
pathToMosek = '';

addpath(genpath(pathToYALMIP));
addpath(genpath(pathToMosek));

%%
clc; clear; close all

%Vehicle properties
rho1 = 4972;
rho2 = 13260;
mDry = 1505;

%Drag Addition
rho = .020;
V = 26.8224+100;
Density = 0.5*rho*V^2;
CA = 1.4;
S = pi/4*(2.65*(1905/608)^(1/3))^2;
Drag = Density*CA*S;
WindAzimuth = rand*2*pi;
DragAccel = Drag/mDry;

%Guidance parameters
guidanceSize = 200;
guidancIterations = 10;
finalTime = 78.4;
guidanceRunRate = finalTime/guidancIterations;

%Initial Conditions
r0 = [1500, 500, 2000]';
rDot0 = [-75, 0, 100]';
mWet = 1905;
alpha = 4.53e-4;
grav = [-3.7114, 0, 0]';
In = [r0; rDot0; mWet; finalTime];
Y0 = [r0; rDot0; log(mWet)];

%Continuous time state matrices
A = [zeros(3) eye(3) zeros(3,1);...
         zeros(3) zeros(3) zeros(3,1);...
         zeros(1,3) zeros(1,3) 0];
Bgrav = [zeros(3); eye(3); zeros(1,3)];

Bu = [zeros(3) zeros(3,1);
          eye(3) zeros(3,1);
          zeros(1,3) -alpha];
B = [Bgrav, Bu];
C = eye(7);
D = 0;
ContinuousModel = ss(A,B,C,D);



%% Looping through each guidance iteration
for i = 1:guidancIterations
    %Calculate the optimal trajectory for the "initial" condition
    [Uvals, deltaT] = calculateOptimalTrajectory(In);

    %Creating a continuous input matrix from the discrete results
    deltaT_c = 0.1;
    indexRatio = ceil(deltaT/deltaT_c);
    deltaT_c = deltaT/indexRatio;
    Time_c(i).vec = [0:deltaT_c:finalTime] + guidanceRunRate*(i-1);
    U_c = [];
    for j = 1:guidanceSize
        U_c = [U_c, Uvals(:,j).*ones(7,indexRatio)];
    end
    U_c = [U_c, zeros(7,1)];

    %Adding in random acceleration due to wind
    if i ==guidancIterations %No noise on final iteration
         randAccel(i).mat = [zeros(7, size(U_c,2))];
    else
        noiseCenterMag = DragAccel; %Random Acceleration between 0 and 0.4 m/s2
        noiseMag = normrnd(noiseCenterMag*ones(1,size(U_c,2)), 0.008);  %Random Acceleration centered around (noiseCenterMag);
        yDirectionNoise = cos(WindAzimuth)*noiseMag;
        zDirectionNoise = sin(WindAzimuth)*noiseMag;
        randAccel(i).mat = [zeros(4, size(U_c,2)); yDirectionNoise;zDirectionNoise ; zeros(1, size(U_c,2))];
    end

    %Running the continuous tome solution
    [XOut] =  lsim(ContinuousModel, U_c+randAccel(i).mat, Time_c(i).vec, Y0);
    PlannedXOut =  lsim(ContinuousModel, U_c, Time_c(i).vec, Y0);
    
    %Saving this iteration's results. 
    output(i).mat = XOut;
    plannedOutput(i).mat = PlannedXOut;

    %Preparing hte input vector of the next iteration
    stopIndex(i) = floor(length(Time_c(i).vec)*guidanceRunRate/finalTime);
    finalTime = finalTime-guidanceRunRate;
    Y0 = [XOut(stopIndex(i),1:3), XOut(stopIndex(i),4:6), XOut(stopIndex(i),7)];
    In = [XOut(stopIndex(i),1:3)'; XOut(stopIndex(i),4:6)'; exp(XOut(stopIndex(i),7)); finalTime];
end



%% Graph Plotting
close all;
figure; hold on; grid on;
figure; hold on; grid on;

CoolMap = flipud(winter(guidancIterations));
SummerMap= flipud(summer(guidancIterations));
AutumnMap= flipud(autumn(guidancIterations));
JetMap= flipud(jet(guidancIterations));
for i = 1:(guidancIterations-1)
    %Figure 1 Position states over time
    figure(1)
    plot(Time_c(i).vec(1:stopIndex(i)), output(i).mat(1:stopIndex(i),1), 'color' ,CoolMap(i,:))
    plot(Time_c(i).vec(1:stopIndex(i)), output(i).mat(1:stopIndex(i),2), 'color' ,SummerMap(i,:))
    plot(Time_c(i).vec(1:stopIndex(i)), output(i).mat(1:stopIndex(i),3), 'color' ,AutumnMap(i,:))
    plot(Time_c(i).vec(1:end), plannedOutput(i).mat(1:end,1), 'color' ,CoolMap(i,:),'LineStyle', '--')
    plot(Time_c(i).vec(1:end), plannedOutput(i).mat(1:end,2), 'color' ,SummerMap(i,:),'LineStyle', '--')
    plot(Time_c(i).vec(1:end), plannedOutput(i).mat(1:end,3), 'color' ,AutumnMap(i,:),'LineStyle', '--')

    %Figure 2 Velocity states over time
    figure(2)
    plot(Time_c(i).vec(1:stopIndex(i)), output(i).mat(1:stopIndex(i),4), 'color' ,CoolMap(i,:))
    plot(Time_c(i).vec(1:stopIndex(i)), output(i).mat(1:stopIndex(i),5), 'color' ,SummerMap(i,:))
    plot(Time_c(i).vec(1:stopIndex(i)), output(i).mat(1:stopIndex(i),6), 'color' ,AutumnMap(i,:))
    plot(Time_c(i).vec(1:end), plannedOutput(i).mat(1:end,4), 'color' ,CoolMap(i,:),'LineStyle', '--')
    plot(Time_c(i).vec(1:end), plannedOutput(i).mat(1:end,5), 'color' ,SummerMap(i,:),'LineStyle', '--')
    plot(Time_c(i).vec(1:end), plannedOutput(i).mat(1:end,6), 'color' ,AutumnMap(i,:),'LineStyle', '--')

end

%Plotting final position, with no wind disturbance for final iteration
figure(1)
plot(Time_c(guidancIterations).vec(1:stopIndex(guidancIterations)), output(guidancIterations).mat(1:stopIndex(guidancIterations),1), 'color' ,CoolMap(guidancIterations,:))
plot(Time_c(guidancIterations).vec(1:stopIndex(guidancIterations)), output(guidancIterations).mat(1:stopIndex(guidancIterations),2), 'color' ,SummerMap(guidancIterations,:))
plot(Time_c(guidancIterations).vec(1:stopIndex(guidancIterations)), output(guidancIterations).mat(1:stopIndex(guidancIterations),3), 'color' ,AutumnMap(guidancIterations,:))

%Plotting final velocity, with no wind disturbance for final iteration
figure(2)
plot(Time_c(guidancIterations).vec(1:stopIndex(guidancIterations)), output(guidancIterations).mat(1:stopIndex(guidancIterations),4), 'color' ,CoolMap(guidancIterations,:))
plot(Time_c(guidancIterations).vec(1:stopIndex(guidancIterations)), output(guidancIterations).mat(1:stopIndex(guidancIterations),5), 'color' ,SummerMap(guidancIterations,:))
plot(Time_c(guidancIterations).vec(1:stopIndex(guidancIterations)), output(guidancIterations).mat(1:stopIndex(guidancIterations),6), 'color' ,AutumnMap(guidancIterations,:))


figure; hold on; grid on;
AutumnMap= flipud(turbo(guidancIterations));
for i = 1:guidancIterations

    figure(3);
    plot3(output(i).mat(1:end,2), output(i).mat(1:end,3), output(i).mat(1:end,1), 'color', AutumnMap(i,:))
    plot3(plannedOutput(i).mat(1:end,2), plannedOutput(i).mat(1:end,3), plannedOutput(i).mat(1:end,1), 'color',AutumnMap(i,:) ,'LineStyle', '--')

end


figure(1);
title('Position states over time')
L1 = plot(nan, nan, 'k');
L2 = plot(nan, nan, 'k--');
legend([L1, L2], {"Position Without Updating again", 'Planned Position'})


figure(2);
title('Velocity states over time')
L1 = plot(nan, nan, 'k');
L2 = plot(nan, nan, 'k--');
legend([L1, L2], {"Velocity Without Updating again", 'Planned Velocity'})

figure(3);
xlabel('East [m]')
ylabel('North [m]');
zlabel('Height [m]');
title('Lander Trajectory');
L1 = plot(nan, nan, 'k');
L2 = plot(nan, nan, 'k--');
legend([L1, L2], {"Trajectory Followed Without Updating again", 'Planned Trajectory'})





