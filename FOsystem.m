%fractional order system to create measurements for kalman filter
%Adrian J Guel C 30/01/2020
function [t,y]=FOsystem(ks,b,N,Q,R,T,t,L,u)
    order=2;
    m=1;
    A=[1,T;-(ks*T/m),(1-b*T/m)];
    B=[0;T/m];
    Ad=A-eye(order,order);
    C=[1,0];
    y=model(Ad,B,C,Q,R,u,N,t,order);

    function y=model(Ad,B,C,Q,R,u,n,t,order)
        x=zeros(order,length(t));
        y=zeros(1,length(t));
        x(:,1)=[1;0];
        for k=2:length(t)
            x(:,k)=Ad*x(:,k-1)+B*u(k-1)-memo(k,x,n,order)+mvnrnd([0,0],Q,1)';
            y(:,k)=C*x(:,k-1)+mvnrnd(0,R,1);
        end
    end

    function r=memo(k,x,n,order)
        r=zeros(order,1);
        if k>=L
            for j=1:L-1
                r=r+((-1)^j)*OrderY(j,n,order)*x(:,k-j);
            end
        else
            for j=1:k-1
                r=r+((-1)^j)*OrderY(j,n,order)*x(:,k-j);
            end
        end
    end

    function b=OrderY(j,n,order)
        b=zeros(order,order);
        if j==0
            b=eye(order,order);
        else
            for i=1:order
                b(i,i)=gamma(n(i)+1)/(gamma(j+1)*gamma(n(i)-j+1));
            end
        end
    end

end
