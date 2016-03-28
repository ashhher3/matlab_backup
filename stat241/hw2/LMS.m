%% import data and store it in matrices of the right shape
data=importdata('hw2.data');
x=data(:,1:2);
y=data(:,3);

%% calculate theta*
theta=(x'*x)\(x'*y);

%% compute eigenvalues and eigenvectors
[V,D]=eig(x'*x);
lambdaMax=max(max(D));

%% plot contours
thetaStep=-0.4:0.001:1.2;

J=nan(length(thetaStep),1);
for i=1:length(thetaStep)
    for j=1:length(thetaStep)
        J(j,i)=(y-x*[thetaStep(i); thetaStep(j)])'*(y-x*[thetaStep(i); thetaStep(j)]);
    end
end

figure(1)
clf
hold on
contour(thetaStep,thetaStep,J,20);
xlabel('theta(1)')
ylabel('theta(2)')

%% LMS paths

thetaT=nan(2,3,31);
thetaT(:,:,1) = zeros(2,3);
rho=[2/lambdaMax; 1/(2*lambdaMax); 1/(8*lambdaMax)];

for i=1:length(y)
    thetaT(:,1,i+1) = thetaT(:,1,i) + rho(1)*(y(i) - thetaT(:,1,i)'*x(i,:)')*x(i,:)';
    thetaT(:,2,i+1) = thetaT(:,2,i) + rho(2)*(y(i) - thetaT(:,2,i)'*x(i,:)')*x(i,:)';
    thetaT(:,3,i+1) = thetaT(:,3,i) + rho(3)*(y(i) - thetaT(:,3,i)'*x(i,:)')*x(i,:)';
end
 
plot(squeeze(thetaT(1,1,:)),squeeze(thetaT(2,1,:)),'r-o');
plot(squeeze(thetaT(1,2,:)),squeeze(thetaT(2,2,:)),'g-o');
plot(squeeze(thetaT(1,3,:)),squeeze(thetaT(2,3,:)),'c-o');