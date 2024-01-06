clear;
close all;
clc;

alpha_u = 1;
alpha_u_dot = 1;
alpha_x = 1;
sigma_1 = 1;
sigma_2 = 1; %1.1148
lambda_2 = 1.0369;
p_control = sqrt((sqrt(sigma_1^2*alpha_u^2+4*sigma_2*alpha_x*alpha_u_dot)-sigma_1*alpha_u)/(2*sigma_2*alpha_u_dot));

L =[3,-1,0,0,0,0,0,0,0,-1,0,0,0,0,0,-1;
    -1,3,-1,0,0,0,0,-1,0,0,0,0,0,0,0,0;
    0,-1,4,-1,0,0,-1,0,0,0,0,0,0,0,-1,0;
    0,0,-1,3,-1,0,0,0,0,0,0,-1,0,0,0,0;
    0,0,0,-1,3,-1,0,0,0,0,0,0,0,0,0,-1;
    0,0,0,0,-1,3,-1,0,0,0,0,0,-1,0,0,0;
    0,0,-1,0,0,-1,4,-1,0,0,-1,0,0,0,0,0;
    0,-1,0,0,0,0,-1,4,-1,0,0,0,0,0,-1,0;
    0,0,0,0,0,0,0,-1,3,-1,0,0,0,-1,0,0;
    -1,0,0,0,0,0,0,0,-1,3,-1,0,0,0,0,0;
    0,0,0,0,0,0,-1,0,0,-1,4,-1,0,0,0,-1;
    0,0,0,-1,0,0,0,0,0,0,-1,3,-1,0,0,0;
    0,0,0,0,0,-1,0,0,0,0,0,-1,3,-1,0,0;
    0,0,0,0,0,0,0,0,-1,0,0,0,-1,3,-1,0;
    0,0,-1,0,0,0,0,-1,0,0,0,0,0,-1,4,-1;
    -1,0,0,0,-1,0,0,0,0,0,-1,0,0,0,-1,4];

I_n = [1, 0; 0, 1];
u = -p_control*(kron(L, I_n));

x1 = [7,7]'; x2 = [8,4]'; x3 = [3,2]'; x4 = [2,5]'; x5 = [-1,1]'; 
x6 = [-7,2]'; x7 = [-1,8]'; x8 = [-5,3]'; x9 = [-8,-2]'; x10 = [-1,-6]'; 
x11 = [-4,-4]'; x12 = [-5,-5]'; x13 = [1,-7]'; x14 = [3,-4]'; x15 = [3,-2]'; 
x16 = [5,-2]';

y = [;];
y(:,1) = [x1',x2',x3',x4',x5',x6',x7',x8',x9',x10',x11',x12',x13',x14',x15',x16']';
t = [];
t(1) = 0;
delta_t = 0.01;
i = 1;
E = 0;

while i < 1000
    E = E + p_control^2*((3*y(1,i)-y(31,i)-y(3,i)-y(19,i))^2+((3*y(2,i)-y(32,i)-y(4,i)-y(20,i))^2))*delta_t;
    y(:,i + 1) = y(:,i) + u*delta_t*y(:,i);
    t(i + 1) = t(i) + delta_t;
    i = i + 1;
    
end

figure(1)
plot(t,y(1,:),'color','#3396CF')
hold on;
plot(t,y(3,:),'Color','#FB8857')
plot(t,y(5,:),'Color','#C37F95')
plot(t,y(7,:),'Color','#79C4D1')

plot(t,y(9,:),'Color','#38375A')
plot(t,y(11,:),'Color','#48948B')
plot(t,y(13,:),'Color','#135A67')
plot(t,y(15,:),'Color','#69356E')
plot(t,y(17,:),'Color','#9C3F2F')
plot(t,y(19,:),'Color','#C4D9A7')
plot(t,y(21,:),'Color','#C630C9')
plot(t,y(23,:),'Color','#3855A5')
plot(t,y(25,:),'Color','#E22C90')
plot(t,y(27,:),'Color','#45287A')
plot(t,y(29,:),'Color','#D64673')
plot(t,y(31,:),'Color','red')
hold off;



% figure(1)
% plot(t,y(1,:),'color','#FD6D5A')
% hold on;
% plot(t,y(3,:),'Color',[0 0 0.80392])
% plot(t,y(5,:),'Color',[0.27843 0.23529 0.5451])
% plot(t,y(7,:),'Color',[0.51373 0.43529 1])
% 
% plot(t,y(9,:),'Color',[0.51373 0.5451 0.5451])
% plot(t,y(11,:),'Color',[1 0.64706 0])
% plot(t,y(13,:),'Color',[0.5451 0.4902 0.48235])
% plot(t,y(15,:),'Color',[0.5451 0.51373 0.52549])
% plot(t,y(17,:),'Color',[0 0 0])
% plot(t,y(19,:),'Color',[1 0.84314 0])
% plot(t,y(21,:),'Color',[0.62745 0.12549 0.94118])
% plot(t,y(23,:),'Color',[1 0 1])
% plot(t,y(25,:),'Color',[1 0 0])
% plot(t,y(27,:),'Color',[0 1 1])
% plot(t,y(29,:),'Color',[1 1 0])
% plot(t,y(31,:),'Color',[0 1 0])
% hold off;
grid;
h1 = legend({'$x_{1-1}$',...
    '$x_{2-1}$','$x_{3-1}$','$x_{4-1}$','$x_{5-1}$','$x_{6-1}$','$x_{7-1}$','$x_{8-1}$','$x_{9-1}$','$x_{10-1}$','$x_{11-1}$','$x_{12-1}$','$x_{13-1}$','$x_{14-1}$','$x_{15-1}$','$x_{16-1}$'},'Location','northeast','NumColumns',2,'Interpreter','latex');

%h1 = legend({'$x_{1-1}$','$x_{2-1}$','$x_{3-1}$','$x_{4-1}$','$x_{5-1}$','$x_{6-1}$','$x_{7-1}$','$x_{8-1}$','$x_{9-1}$','$x_{10-1}$','$x_{11-1}$','$x_{12-1}$','$x_{13-1}$','$x_{14-1}$','$x_{15-1}$','$x_{16-1}$'},'Location','northwest','Interpreter','latex');
set(h1,'FontName','Times New Roman','FontSize',11,'FontWeight','normal')
set(h1,'Color','none');
set(h1,'Box','off')
xlabel('Times(s)','Interpreter','latex');
ylabel('$x_{i1},i=1,...,16$','Interpreter','latex');
set(gca,'Fontsize',13)

figure(2)
plot(t,y(1,:),'color','#3396CF')
hold on;
plot(t,y(3,:),'Color','#FB8857')
plot(t,y(5,:),'Color','#C37F95')
plot(t,y(7,:),'Color','#79C4D1')

plot(t,y(9,:),'Color','#38375A')
plot(t,y(11,:),'Color','#48948B')
plot(t,y(13,:),'Color','#135A67')
plot(t,y(15,:),'Color','#69356E')
plot(t,y(17,:),'Color','#9C3F2F')
plot(t,y(19,:),'Color','#C4D9A7')
plot(t,y(21,:),'Color','#C630C9')
plot(t,y(23,:),'Color','#3855A5')
plot(t,y(25,:),'Color','#E22C90')
plot(t,y(27,:),'Color','#45287A')
plot(t,y(29,:),'Color','#D64673')
plot(t,y(31,:),'Color','red')
hold off;
grid;
h2 = legend({'$x_{1-2}$','$x_{2-2}$','$x_{3-2}$','$x_{4-2}$','$x_{5-2}$','$x_{6-2}$','$x_{7-2}$','$x_{8-2}$','$x_{9-2}$','$x_{10-2}$','$x_{11-2}$','$x_{12-2}$','$x_{13-2}$','$x_{14-2}$','$x_{15-2}$','$x_{16-2}$'},'Location','northeast','NumColumns',2,'Interpreter','latex');
set(h2,'FontName','Times New Roman','FontSize',11,'FontWeight','normal')
set(h2,'Color','none');
set(h2,'Box','off')
xlabel('Times(s)','Interpreter','latex');
ylabel('$x_{i2},i=1,...,16$','Interpreter','latex');
set(gca,'Fontsize',13)







t_f = (1/(2*lambda_2*p_control))*log(15*y(:,1)'*kron(L,I_n)*y(:,1)/0.01);
delta_max = [8,8]'; delta_min=[-9,-7]';
C2 = delta_max - delta_min;
C11= x1 - delta_max;C12 = x2 - delta_max;C13 = x3 - delta_max;C14 = x4 - delta_max;C15 = (x5 - delta_max)/8;C16 = x6 - delta_max;C17 = x7 - delta_max;C18 = x8 - delta_max;C19 = x9 - delta_max;C110 = x10 - delta_max;C111 = x11 - delta_max;C112 = x12 - delta_max;C113 = x13 - delta_max;C114 = x14 - delta_max;C115 = x15 - delta_max;C116 = x16 - delta_max;
C11'*C11/(2*15*p_control)*(1-exp(-2*15*p_control*t_f));
2*C11'*x1/(15*p_control)*(1-exp(-3*p_control*t_f));
x1'*x1*t_f;
Jz1 = C11'*C11/(2*15*p_control)*(1-exp(-2*15*p_control*t_f))+2*C11'*C2/(15*p_control)*(1-exp(-3*p_control*t_f))+C2'*C2*t_f/2;
Jz2 = C12'*C12/(2*15*p_control)*(1-exp(-2*15*p_control*t_f))+2*C12'*C2/(15*p_control)*(1-exp(-15*p_control*t_f))+C2'*C2*t_f;
Jz3 = C13'*C13/(2*15*p_control)*(1-exp(-2*15*p_control*t_f))+2*C13'*C2/(15*p_control)*(1-exp(-15*p_control*t_f))+C2'*C2*t_f;
Jz4 = C14'*C14/(2*15*p_control)*(1-exp(-2*15*p_control*t_f))+2*C14'*C2/(15*p_control)*(1-exp(-15*p_control*t_f))+C2'*C2*t_f;
Jz5 = C15'*C15/(2*3*p_control)*(1-exp(-2*15*p_control*t_f))/64+2*C15'*x5/(3*p_control)*(1-exp(-15*p_control*t_f))+x5'*x5/8*t_f;
Jz6 = C16'*C16/(2*15*p_control)*(1-exp(-2*15*p_control*t_f))+2*C16'*C2/(15*p_control)*(1-exp(-15*p_control*t_f))+C2'*C2*t_f;
Jz7 = C17'*C17/(2*15*p_control)*(1-exp(-2*15*p_control*t_f))+2*C17'*C2/(15*p_control)*(1-exp(-15*p_control*t_f))+C2'*C2*t_f;
Jz8 = C18'*C18/(2*15*p_control)*(1-exp(-2*15*p_control*t_f))+2*C18'*C2/(15*p_control)*(1-exp(-15*p_control*t_f))+C2'*C2*t_f;
Jz9 = C19'*C19/(2*15*p_control)*(1-exp(-2*15*p_control*t_f))+2*C19'*C2/(15*p_control)*(1-exp(-15*p_control*t_f))+C2'*C2*t_f;
Jz10 = C110'*C110/(2*15*p_control)*(1-exp(-2*15*p_control*t_f))+2*C110'*C2/(15*p_control)*(1-exp(-15*p_control*t_f))+C2'*C2*t_f;
Jz11 = C111'*C111/(2*15*p_control)*(1-exp(-2*15*p_control*t_f))+2*C111'*C2/(15*p_control)*(1-exp(-15*p_control*t_f))+C2'*C2*t_f;
Jz12 = C112'*C112/(2*15*p_control)*(1-exp(-2*15*p_control*t_f))+2*C112'*C2/(15*p_control)*(1-exp(-15*p_control*t_f))+C2'*C2*t_f;
Jz13 = C113'*C113/(2*15*p_control)*(1-exp(-2*15*p_control*t_f))+2*C113'*C2/(15*p_control)*(1-exp(-15*p_control*t_f))+C2'*C2*t_f;
Jz14 = C114'*C114/(2*15*p_control)*(1-exp(-2*15*p_control*t_f))+2*C114'*C2/(15*p_control)*(1-exp(-15*p_control*t_f))+C2'*C2*t_f;
Jz15 = C115'*C115/(2*15*p_control)*(1-exp(-2*15*p_control*t_f))+2*C115'*C2/(15*p_control)*(1-exp(-15*p_control*t_f))+C2'*C2*t_f;
Jz16 = C116'*C116/(2*15*p_control)*(1-exp(-2*15*p_control*t_f))+2*C116'*C2/(15*p_control)*(1-exp(-15*p_control*t_f))+C2'*C2*t_f;

JE1 = p_control^2*9*Jz1+ p_control^2*9*(C2'*C2)*t_f/2
JE2 = p_control^2*225*Jz2+ p_control^2*225*(C2'*C2)*t_f;
JE3 = p_control^2*225*Jz3+ p_control^2*225*(C2'*C2)*t_f;
JE4 = p_control^2*225*Jz4+ p_control^2*225*(C2'*C2)*t_f;


JE5 = p_control^2*9*Jz5+ p_control^2*9*(x5'*x5)*t_f;


JE6 = p_control^2*225*Jz6+ p_control^2*225*(C2'*C2)*t_f;
JE7 = p_control^2*225*Jz7+ p_control^2*225*(C2'*C2)*t_f;
JE8 = p_control^2*225*Jz8+ p_control^2*225*(C2'*C2)*t_f;
JE9 = p_control^2*225*Jz9+ p_control^2*225*(C2'*C2)*t_f;
JE10 = p_control^2*225*Jz10+ p_control^2*225*(C2'*C2)*t_f;
JE11 = p_control^2*225*Jz11+ p_control^2*225*(C2'*C2)*t_f;
JE12 = p_control^2*225*Jz12+ p_control^2*225*(C2'*C2)*t_f;
JE13 = p_control^2*225*Jz13+ p_control^2*225*(C2'*C2)*t_f;
JE14 = p_control^2*225*Jz14+ p_control^2*225*(C2'*C2)*t_f;
JE15 = p_control^2*225*Jz15+ p_control^2*225*(C2'*C2)*t_f;
JE16 = p_control^2*225*Jz16+ p_control^2*225*(C2'*C2)*t_f;
E

figure(5)

y_5 = [9.087,6.41;22.6989,16.05;13.7423,9.7;4.3457,3.06;8.0836,5.7;8.3482,5.9];
b = bar(y_5)
grid;
ax = gca;
ax.XTickLabel = {'A','B','C','D','E','F'}
h = legend('$\overline{t_f}$','$t_f$','Interpreter','latex');
ylabel('Times(s)','Interpreter','latex');
set(h,'FontName','Times New Roman','FontSize',11,'FontWeight','normal')
set(h,'Color','none');
set(h,'Box','off')
set(gca,'Fontsize',14)


figure(6)
y = [5.7,1.97;5.4,1.57;5.62,1.79;6.12,2.3;5.8,2.02;5.84,2];
b = bar(y)
grid;
ax = gca;
ax.XTickLabel = {'A','B','C','D','E','F'}

h = legend(' $\log_{10}{\overline{J_{E_1}}}$','$\log_{10}{J_{E_1}}$','Interpreter','latex','Location','northwest');
ylabel('Energy','Interpreter','latex');
set(gca,'Fontsize',14)
set(h,'FontName','Times New Roman','FontSize',11,'FontWeight','normal')
set(h,'Color','none');
set(h,'Box','off')





alpha_u_s1= 1;
alpha_u_dot_s1 = 1;
alpha_x_s1 = 1;
sigma_1_s1(1) = 0;
sigma_2_s1 = 1; %1.1148
lambda_2 = 1.0369;


L =[3,-1,0,0,0,0,0,0,0,-1,0,0,0,0,0,-1;
    -1,3,-1,0,0,0,0,-1,0,0,0,0,0,0,0,0;
    0,-1,4,-1,0,0,-1,0,0,0,0,0,0,0,-1,0;
    0,0,-1,3,-1,0,0,0,0,0,0,-1,0,0,0,0;
    0,0,0,-1,3,-1,0,0,0,0,0,0,0,0,0,-1;
    0,0,0,0,-1,3,-1,0,0,0,0,0,-1,0,0,0;
    0,0,-1,0,0,-1,4,-1,0,0,-1,0,0,0,0,0;
    0,-1,0,0,0,0,-1,4,-1,0,0,0,0,0,-1,0;
    0,0,0,0,0,0,0,-1,3,-1,0,0,0,-1,0,0;
    -1,0,0,0,0,0,0,0,-1,3,-1,0,0,0,0,0;
    0,0,0,0,0,0,-1,0,0,-1,4,-1,0,0,0,-1;
    0,0,0,-1,0,0,0,0,0,0,-1,3,-1,0,0,0;
    0,0,0,0,0,-1,0,0,0,0,0,-1,3,-1,0,0;
    0,0,0,0,0,0,0,0,-1,0,0,0,-1,3,-1,0;
    0,0,-1,0,0,0,0,-1,0,0,0,0,0,-1,4,-1;
    -1,0,0,0,-1,0,0,0,0,0,-1,0,0,0,-1,4];

I_n = [1, 0; 0, 1];


x1 = [7,7]'; x2 = [8,4]'; x3 = [3,2]'; x4 = [2,5]'; x5 = [-1,1]'; 
x6 = [-7,2]'; x7 = [-1,8]'; x8 = [-5,3]'; x9 = [-9,-2]'; x10 = [-1,-6]'; 
x11 = [-4,-4]'; x12 = [-5,-5]'; x13 = [1,-7]'; x14 = [3,-4]'; x15 = [4,-2]'; 
x16 = [5,-2]';

y_s1 = [;];
y_s1(:,1) = [x1',x2',x3',x4',x5',x6',x7',x8',x9',x10',x11',x12',x13',x14',x15',x16']';
t_s1 = [];
t_s1(1) = 0;
delta_t = 0.01;
i = 1;
t_f = [];
j = 1;
p_control_s1 = [];
while i < 103 
    p_control_s1(i) = sqrt((sqrt(sigma_1_s1(i)^2*alpha_u_s1^2+4*sigma_2_s1*alpha_x_s1*alpha_u_dot_s1)-sigma_1_s1(i)*alpha_u_s1)/(2*sigma_2_s1*alpha_u_dot_s1)); 
    sigma_1_s1(i+1)= sigma_1_s1(i) + 0.01;
    i = i + 1;
end
p_control_s2 = [];
alpha_u_s2= 1;
alpha_u_dot_s2 = 1;
alpha_x_s2 = 1;
sigma_1_s2 = 1;
sigma_2_s2(1) = 0; %1.1148u
while j < 111
    p_control_s2(j) = sqrt((sqrt(sigma_1_s2^2*alpha_u_s2^2+4*sigma_2_s2(j)*alpha_x_s2*alpha_u_dot_s2)-sigma_1_s2*alpha_u_s2)/(2*sigma_2_s2(j)*alpha_u_dot_s2)); 
    sigma_2_s2(j+1)= sigma_2_s2(j) + 0.01;
    j = j + 1;
end

figure(3)
plot(sigma_1_s1(1,1:102),p_control_s1(1,:),sigma_2_s2(1,1:110),p_control_s2(1,:));
grid;
xlabel('$\sigma_1,\sigma_2$','Interpreter','latex');ylabel('$p^{*}$','Interpreter','latex');
h = legend('$\sigma_1$','$\sigma_2$','Interpreter','latex');
set(gca,'Fontsize',14)
set(h,'FontName','Times New Roman','FontSize',11,'FontWeight','normal')
set(h,'Color','none');
set(h,'Box','off')











alpha_u_u(1)= 0;
alpha_u_dot_u = 1;
alpha_x_u = 1;
sigma_1_u = 1;
sigma_2_u = 1; %1.1148
lambda_2 = 1.0369;


L =[3,-1,0,0,0,0,0,0,0,-1,0,0,0,0,0,-1;
    -1,3,-1,0,0,0,0,-1,0,0,0,0,0,0,0,0;
    0,-1,4,-1,0,0,-1,0,0,0,0,0,0,0,-1,0;
    0,0,-1,3,-1,0,0,0,0,0,0,-1,0,0,0,0;
    0,0,0,-1,3,-1,0,0,0,0,0,0,0,0,0,-1;
    0,0,0,0,-1,3,-1,0,0,0,0,0,-1,0,0,0;
    0,0,-1,0,0,-1,4,-1,0,0,-1,0,0,0,0,0;
    0,-1,0,0,0,0,-1,4,-1,0,0,0,0,0,-1,0;
    0,0,0,0,0,0,0,-1,3,-1,0,0,0,-1,0,0;
    -1,0,0,0,0,0,0,0,-1,3,-1,0,0,0,0,0;
    0,0,0,0,0,0,-1,0,0,-1,4,-1,0,0,0,-1;
    0,0,0,-1,0,0,0,0,0,0,-1,3,-1,0,0,0;
    0,0,0,0,0,-1,0,0,0,0,0,-1,3,-1,0,0;
    0,0,0,0,0,0,0,0,-1,0,0,0,-1,3,-1,0;
    0,0,-1,0,0,0,0,-1,0,0,0,0,0,-1,4,-1;
    -1,0,0,0,-1,0,0,0,0,0,-1,0,0,0,-1,4];

I_n = [1, 0; 0, 1];


x1 = [7,7]'; x2 = [8,4]'; x3 = [3,2]'; x4 = [2,5]'; x5 = [-1,1]'; 
x6 = [-7,2]'; x7 = [-1,8]'; x8 = [-5,3]'; x9 = [-9,-2]'; x10 = [-1,-6]'; 
x11 = [-4,-4]'; x12 = [-5,-5]'; x13 = [1,-7]'; x14 = [3,-4]'; x15 = [4,-2]'; 
x16 = [5,-2]';

y_u = [;];
y_u(:,1) = [x1',x2',x3',x4',x5',x6',x7',x8',x9',x10',x11',x12',x13',x14',x15',x16']';
t = [];
t(1) = 0;
delta_t = 0.01;
i = 1;
t_f = [];
j = 1;
p_control_u = [];
while j < 10000 
    p_control_u(j) = sqrt((sqrt(sigma_1_u^2*alpha_u_u(j)^2+4*sigma_2_u*alpha_x_u*alpha_u_dot_u)-sigma_1_u*alpha_u_u(j))/(2*sigma_2_u*alpha_u_dot_u)); 
    alpha_u_u(j+1)= alpha_u_u(j) + 0.01;
    j = j + 1;
end
p_control_dot = [];
alpha_u_dot= 1;
alpha_u_dot_dot(1) = 0;
alpha_x_dot = 1;
sigma_1_dot = 1;
sigma_2_dot = 1; %1.1148u
while i < 10000
    p_control_dot(i) = sqrt((sqrt(sigma_1_dot^2*alpha_u_dot^2+4*sigma_2_dot*alpha_x_dot*alpha_u_dot_dot(i))-sigma_1_dot*alpha_u_dot)/(2*sigma_2_dot*alpha_u_dot_dot(i))); 
    alpha_u_dot_dot(i+1)= alpha_u_dot_dot(i) + 0.01;
    i = i + 1;
end
k = 1;
p_control_x = [];
alpha_u_x= 1;
alpha_u_dot_x = 1;
alpha_x_x(k) = 0;
sigma_1_x = 1;
sigma_2_x = 1; %1.1148u
while k < 10000
    p_control_x(k) = sqrt((sqrt(sigma_1_x^2*alpha_u_x^2+4*sigma_2_x*alpha_x_x(k)*alpha_u_dot_x)-sigma_1_x*alpha_u_x)/(2*sigma_2_x*alpha_u_dot_x)); 
    alpha_x_x(k+1)= alpha_x_x(k) + 0.01;
    k = k + 1;
end


figure(4)
plot(alpha_u_u(1,1:9999),p_control_u(1,:),alpha_u_dot_dot(1,1:9999),p_control_dot(1,:),alpha_x_x(1,1:9999),p_control_x(1,:));
grid;
xlabel('$\alpha_u,\alpha_{\dot{u}},\alpha_x$','Interpreter','latex');ylabel('$p^{*}$','Interpreter','latex');
h = legend('$\alpha_u$','$\alpha_{\dot{u}}$','$\alpha_x$','Interpreter','latex');

set(gca,'Fontsize',14)
set(h,'FontName','Times New Roman','FontSize',11,'FontWeight','normal')
set(h,'Color','none');
set(h,'Box','off')









