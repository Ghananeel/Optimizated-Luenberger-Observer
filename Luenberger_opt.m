%% State estimation of inverted pendulum
% Author : Ghananeel Rotithor

% Simulation (Inverted Pendulum) parameters
M = .5;
m = 0.2;
b = 0.1;
I = 0.006;
g = 9.8;
l = 0.3;

p = I*(M+m)+M*m*l^2; %denominator for the A and B matrices

% State Transition Matrix
A = [0      1              0           0;
     0 -(I+m*l^2)*b/p  (m^2*g*l^2)/p   0;
     0      0              0           1;
     0 -(m*l*b)/p       m*g*l*(M+m)/p  0];
 
% Input Gain Matrix
B = [0; (I+m*l^2)/p; 0; m*l/p];

% Output Matrix and LQR parameters
C = [1 0 0 0; 0 0 1 0]; Q = ones(size(A,1),size(A,1)); R = 1; 

n = size(A,1); p = size(C,1);

% Optimize for P and W
cvx_begin sdp
    variable P(n,n) symmetric
    variable W(n,p)
    A'*P+P*A-W*C-C'*W'+0.2*P<=zeros(n,n);
    P >= 0.1*eye(n);
cvx_end

% Compute Observer Gain L
L = inv(P)*W;

dt = 0.1;
t = 0:dt:10;
Px = zeros(size(A,1),size(A,1),length(t));

% Create state and initialize
x = zeros(size(A,1),length(t));
x(:,1) = -0.1*rand(size(A,1),1);

% Create estimate and initialize
x_hat = zeros(size(A,1),length(t)); 
x_hat(:,1) = 0.1*rand(size(A,1),1);

for i = 2:length(t)
    % Forward Propagate Continuous Riccati Equation
    P_dot =  Px(:,:,i-1)*A + A'*Px(:,:,i-1) - Px(:,:,i-1)*(B*B')*Px(:,:,i-1) + Q;
    Px(:,:,i) = Px(:,:,i-1) + P_dot*dt;
    % Compute Optimal Input Gain and Input
    K = B'*Px(:,:,i);
    u(i-1) = -K*x(:,i-1);
    
    % Forward Propage Dynamical System with True State
    x_dot = A*x(:,i-1)+B*u(i-1);
    x(:,i) = x(:,i-1) + x_dot * dt;
    
    % Forward Propage Observer with Estimate
    xhat_dot = A*x_hat(:,i-1)+B*u(i-1)+L*(C*x(:,i-1)-C*x_hat(:,i-1));
    x_hat(:,i) = x_hat(:,i-1) + xhat_dot * dt;
    
    % Compute Lyapunov Funtion Value
    V(i) = (x(:,i)-x_hat(:,i))'*P*(x(:,i)-x_hat(:,i));
end
V(1) = (x(:,1)-x_hat(:,1))'*P*(x(:,1)-x_hat(:,1));

% Plotting script
figure
subplot(131)
plot(t,x','LineWidth',2); hold on;
plot(t,x_hat','--','LineWidth',2);
legend('$x$','$\dot{x}$','$\phi$','$\dot{\phi}$',...
        '$\hat{x}$','$\dot{\hat{x}}$','$\hat{\phi}$','$\dot{\hat{\phi}}$',...
        'Interpreter','latex','FontSize',18);
xlabel('Time(s)','Interpreter','latex','FontSize',18);
ylabel('State Values','Interpreter','latex','FontSize',18);

subplot(132)
plot(t,V,'LineWidth',2); 
xlabel('Time(s)','Interpreter','latex','FontSize',18);
ylabel('$e^{T}Pe$','Interpreter','latex','FontSize',18);
legend('$V(e)$','Interpreter','latex','FontSize',18);

subplot(133)
plot(t(1:end-1),u,'LineWidth',2); 
xlabel('Time(s)','Interpreter','latex','FontSize',18);
ylabel('Force $F$','Interpreter','latex','FontSize',18);
legend('$u$','Interpreter','latex','FontSize',18);


