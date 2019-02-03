%% 1 Electron Modelling

global C

C.q_0 = 1.60217653e-19;             % electron charge
C.m_0 = 9.10938215e-31;             % electron mass
C.mn = 0.26*C.m_0;                  % Effective Electron Mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.T = 300;                          % Kelvin

vth = sqrt(2*C.kb*C.T/C.mn); %Thermal velocity
MTBC = 0.2e-12; %Mean time between colissions (s)
MFP = vth*MTBC; %Mean free path (m)

numPart = 10000; %Number of particles
xlim = 200e-9;
ylim = 100e-9;
dt = ylim/vth/100; %Scale time

%Random starting positions
x=rand(1,numPart)*xlim;
y=rand(1,numPart)*ylim;

%Random angle
hAngle = 360; %highest angle
lAngle = 0; %lowest angle
angle = (hAngle-lAngle).*rand(1,numPart) + lAngle; %Random angle within range

%Component velocities travelling at random angle
vx = cos(angle)*vth;
vy = sin(angle)*vth;

%Subset of particles to be plotted
numPartPlot = 7;
subset = randi([1,numPart],numPartPlot,1).';

maxTime = 1000;
for time=1:maxTime

    %y boundaries
    yBoundTop = y >= ylim;
    y(yBoundTop) = ylim;
    yBoundBottom = y<=0;
    y(yBoundBottom) = 0;
    yBound = yBoundTop | yBoundBottom;
    vy(yBound) = -1.*vy(yBound); %Reverse velocity
    
    %Updating y position
    yPrev = y;
    y = y + vy*dt; 
    
    %x boundaries
    rightBound = (x>=xlim & vx>=0); %Positive xvelocities reaching right boundary
    x(rightBound) = 0; %Relocate particle to left side
    leftBound = (x<=0 & vx<=0); %Negative xvelocities reaching left boundary
    x(leftBound) = xlim; %Relocate particle to right side
    
    %Updating x position
    xPrev = x;
    x = x + vx*dt; 
    
    %Plotting
    for i=1:numPartPlot
        plot([xPrev(subset(i)) x(subset(i))],[yPrev(subset(i)) y(subset(i))])
    end
    title('Particles Trajectories at Vth')
    axis ([0 xlim 0 ylim])
    drawnow
    hold on

    %Semiconductor temperature
    v = sqrt(vx.^2+vy.^2);
    overallTemp = C.mn*sum(v.^2)/(2*C.kb);
    avgTemp(time) = overallTemp/numPart;  
end

%Plotting temperature
xtime = linspace(0,maxTime,maxTime);
figure(2)
plot(xtime,avgTemp)
axis ([0 maxTime 200 400])
title('Average Semiconductor Temperature Over Time')
xlabel('Time (s)')
ylabel('Temperature (K)')
hold on

% Figure 2 shows that the semiconductor temperature remains constant since 
% all particles are travelling at the same velocity.

%Displaying final calcs
fprintf('Part 1:',vth);
fprintf('\nThe thermal velocity is %d m/s.',vth);
fprintf('\nThe Mean Free Path is %d m.\n',MFP);




