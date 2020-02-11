%--------------- ELLIPTIC Mesh Generation---------
%-------------------------------------------------
clear all
clc
close all
%---------------parameters--------------------
nx=24;
ny=18;
maxit=1800;
show=1;            % 1 for yes o for no  to disp solution while solving
Ermax=10^-5;

%---------initializing------------------

f = @(x) x.^2;

c1=[-1+zeros(1,ny);linspace(0,4,ny)];        %left
c2=[linspace(-1,0,nx);zeros(1,nx)];          %bottom
c3=[linspace(0,2,ny);f(linspace(0,2,ny))];   %right
c4=[linspace(-1,2,nx);4+zeros(1,nx)];        %top
alpha=zeros(nx,ny);
beta=zeros(nx,ny);
gamma=zeros(nx,ny);

X=zeros(nx,ny);

X(1,:)=c3(1,:);           %right
X(nx,:)=c1(1,:);          %left
X(:,ny)=flipud(c4(1,:)');  %top
X(:,1)=flipud(c2(1,:)');   %bottom
Y=zeros(nx,ny);

 Y(1,:)=c3(2,:);           %right
 Y(nx,:)=c1(2,:);          %left
 Y(:,ny)=flipud(c4(2,:)'); %top
 Y(:,1)=flipud(c2(2,:)');  %bottom

newX=X; newY=Y;
Er1=zeros(1,maxit);
Er2=zeros(1,maxit);
%------------------------------------
%---calculating by iterations--------
%------------------------------------

for t=1:maxit
    
i=2:nx-1;
j=2:ny-1;

alpha(i,j)=(1/4)*((X(i,j+1)-X(i,j-1)).^2 + (Y(i,j+1) - Y(i,j-1)).^2);

beta(i,j)=(1/16)*((X(i+1,j)-X(i-1,j)).*(X(i,j+1)-X(i,j-1))+...
(Y(i+1,j) - Y(i-1,j)).*(Y(i,j+1) - Y(i,j-1)));

gamma(i,j)=(1/4)*((X(i+1,j)-X(i-1,j)).^2 + (Y(i+1,j) - Y(i-1,j)).^2);

newX(i,j)=((-0.5)./(alpha(i,j)+gamma(i,j)+10^-9)).*(2*beta(i,j).*(X(i+1,j+1)-X(i-1,j+1)-X(i+1,j-1) + X(i-1,j-1))...
-alpha(i,j).*(X(i+1,j)+X(i-1,j))-gamma(i,j).*(X(i,j+1)+X(i,j-1)));

newY(i,j)=((-0.5)./(alpha(i,j)+gamma(i,j)+10^-9)).*(2*beta(i,j).*(Y(i+1,j+1)-Y(i-1,j+1)-Y(i+1,j-1) + Y(i-1,j-1))...
-alpha(i,j).*(Y(i+1,j)+Y(i-1,j))-gamma(i,j).*(Y(i,j+1)+Y(i,j-1)));

Er1(1,t)=max(max(abs(newX-X)));
Er2(1,t)=max(max(abs(newY-Y)));

% Neuman BC
newY(nx,:)= newY(nx-1,:);     %left

%-------- dont use  -----------
% just an example to show that b.c are tricky

% newY(1,:)= newY(2,:);         %rigth

%--------------------------
X=newX;
Y=newY;
if Er1(t)<Ermax &&Er2(t)<Ermax
  break
end
if show==1
  if ceil(t/10)*10==t
    clf 
    hold on
    axis equal
    for m=1:nx
    plot(X(m,:),Y(m,:),'b');
    end
    for m=1:ny
    plot(X(:,m),Y(:,m),'Color',[0 0 0]);
    end
    pause(0.001)
  end
end
end
if t==maxit
warning('convergence not reached')
end
clf
hold on
axis equal
for m=1:nx
plot(X(m,:),Y(m,:),'b');
end
for m=1:ny
plot(X(:,m),Y(:,m),'Color',[0 0 0]);
end
