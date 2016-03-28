function Vhat=estimateValue(alpha,beta,V)
n=length(V);

Vhat=V;
Vhat(1,:)=[1 1 1 1 1];

for i=2:n
    Vhat(i,1)=(1-alpha)*Vhat(i-1,1)+alpha*V(i-1,1);
    Vhat(i,2)=(1-alpha)*Vhat(i-1,2)+alpha*V(i-1,2);
    Vhat(i,3)=(1-beta)*Vhat(i-1,3)+beta*V(i-1,3);
    Vhat(i,4)=(1-beta)*Vhat(i-1,4)+beta*V(i-1,4);
end

end