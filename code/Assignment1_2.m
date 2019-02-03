%% 2 Collisions with Mean Free Path (MFP)

clear

global C

C.q_0 = 1.60217653e-19;             % electron charge
C.m_0 = 9.10938215e-31;             % electron mass
C.mn = 0.26*C.m_0;                  % Effective Electron Mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.T = 300;                          % Kelvin

vth = sqrt(2*C.kb*C.T/C.mn); %Thermal velocity
MTBC1 = 0.2e-12; %Mean time between colissions (s)

numPart = 10000; %Number of particles
xlim = 200e-9;
ylim = 100e-9;
bin = numPart/10;
dt = ylim/vth/100; %Scale time

%Random starting positions
x=rand(1,numPart)*xlim;
y=rand(1,numPart)*ylim;

%Random angle
hAngle = 360; %highest angle
lAngle = 0; %lowest angle
angle = (hAngle-lAngle).*rand(1,numPart) + lAngle; %Random angle within range

%Random MB velocity
MBfact = vth;
MBvx = randn(1,numPart)*MBfact;
MBvy = randn(1,numPart)*MBfact;
MBv = sqrt(MBvx.^2+MBvy.^2);

%Random MB velocities travelling at random angle
vx = MBvx.*cos(angle);
vy = MBvy.*sin(angle);

%Plotting velocity distributions
figure(3)
subplot(3,1,1);
hist(MBvx,bin)
title('Hysteresis Plot of Velocity x Component')
xlabel('Velocity (m/s)')

subplot(3,1,2);
hist(MBvy,bin)
title('Hysteresis Plot of Velocity y Component')
xlabel('Velocity (m/s)')

subplot(3,1,3);
hist(MBv,bin)
title('Hysteresis Plot of Velocity')
xlabel('Velocity (m/s)')

% Figure 3 shows the Maxwell-Boltzmann distribution for velocity.

%Scatter probability
Pscat = 1-exp(-dt/MTBC1);

%Initialize Mean Free Path vector
MFPs = zeros(1,numPart);

%Subset of particles to be plotted
numPartPlot = 7;
subset = randi([1,numPart],numPartPlot,1).';

maxTime = 500;
for time=1:maxTime
    %Scattering
    scatter = Pscat > rand(1,numPart); %Particles that will scatter
    vx(scatter) = 0;
    vy(scatter) = 0;    
    angle = (hAngle-lAngle).*rand(1,numPart) + lAngle; %Random angle within range
    
    %Random MB velocity
    MBvx = randn(1,numPart)*MBfact; %New x component velocity
    MBvy = randn(1,numPart)*MBfact; %New y component velocity
    MBv = sqrt(MBvx.^2+MBvy.^2); %New actual particle velcoity
   
    %Random MB velocities travelling at random angle
    vx(scatter) = MBvx(scatter).*cos(angle(scatter));
    vy(scatter) = MBvy(scatter).*sin(angle(scatter));

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
    figure(4)
    title('Particles Trajectories Using Maxwell-Boltzmann Distribution')
    for i=1:numPartPlot
        plot([xPrev(subset(i)) x(subset(i))],[yPrev(subset(i)) y(subset(i))])
    end
    axis ([0 xlim 0 ylim])
    drawnow
    hold on
    
    %Semiconductor temperature
    v = sqrt(vx.^2+vy.^2);
    overallTemp = C.mn*sum(v.^2)/(2*C.kb);
    avgTemp(time) = overallTemp/numPart;  
    
    %Mean Free Path
    MFPs(scatter) = 0;
    notScatter = ismissing(scatter,0); %Prep index for particles that have not scattered
    MFPs(notScatter) =  MFPs(notScatter) + v(notScatter)*dt; %Add distance travelled by particles if not scattered
    MFP = sum(MFPs)/numPart;
    
    %Mean Time Between Collisions
    MTBC = MFP*numPart/sum(v);
end

%Plotting Temperature
xtime = linspace(0,maxTime,maxTime);
figure(5)
plot(xtime,avgTemp)
axis ([0 maxTime 200 400])
title('Average Semiconductor Temperature Over Time')
xlabel('Time (s)')
ylabel('Temperature (K)')
hold on

% The average temperature of the system has an average of 300K given that
% the individual particles share the same average velocity of Vth centred
% on a Maxwell-Boltzmann distribution.

%Displaying final calcs
fprintf('Part 2:',vth);
fprintf('\nThe Mean free path is %d m/s.',MFP);
fprintf('\nThe Mean Time Between Collisions is %d m.\n',MTBC);
