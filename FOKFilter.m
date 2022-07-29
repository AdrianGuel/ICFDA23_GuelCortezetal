function varargout= FOKFilter(t,u,y,x_i,Q,R,T,L)
%fractional order kalman filter
%Adrian J Guel C 28/01/2020

n=5;
P_k=zeros(n,n,length(t));
P_k(:,:,1)=500*eye(n,n);
H=[1,0,0,0,0];
x_e=zeros(n,length(t));
x_e(3:5,1)=x_i;
alpha=0.03;

for k=2:length(t)
    x_model=f(x_e(:,k-1),u(k-1),T)-memo(k,x_e,[x_e(n,k-1);x_e(n,k-1);1;1;1],n);
    P_k(:,:,k)=(Fk(x_e(:,k-1),T)+OrderY(1,[x_e(n,k-1);x_e(n,k-1);1;1;1],n))*P_k(:,:,k-1)*(Fk(x_e(:,k-1),T)+OrderY(1,[x_e(n,k-1);x_e(n,k-1);1;1;1],n))'+Q+memoP(k,P_k,[x_e(n,k-1);x_e(n,k-1);1;1;1],n);
    K_k=P_k(:,:,k)*H'*(H*P_k(:,:,k)*H'+R)^(-1);
    x_e(:,k)=x_model+K_k*(y(k)-x_model(1));
    P_k(:,:,k)=(eye(n,n)-K_k*H)*P_k(:,:,k);
    Q=(1-alpha)*Q+alpha*K_k*((y(:,k)-H*x_e(:,k))*(y(:,k)-H*x_e(:,k))')*K_k';
end

[te,ye]=FOsystem(x_e(3,end),x_e(4,end),[x_e(5,end);x_e(5,end)],0*eye(2,2),0,T,t,L,u);
RMSE = norm(y(1,:) - ye(1,:),2)^2+0.1*norm(x_e(3:5,end),2)^2; 
%RMSE = norm(y(1,:) - ye(1,:),2)^2;

    if abs(nargout)==1
        varargout={RMSE};
    else
        varargout={RMSE,te,ye,x_e};
    end

function y=Fk(x,T)
    y=[0,T,0,0,0;
        -T*x(3),-T*x(4),-T*x(1),-T*x(2),0;
        0,0,0,0,0;
        0,0,0,0,0;
        0,0,0,0,0];
end

function y=f(x,F,T)
    y=[T*x(2);
        -T*x(3)*x(1)-T*x(4)*x(2)+T*F;
        0;
        0;
        0];
end


function r=memo(k,x,n,order)
    r=zeros(order,1);
    %L=200;
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

function r=memoP(k,P,N,order)
    r=zeros(order,order);
    %L=200;
    if k>=L
        for j=2:L-1
            r=r+OrderY(j,N,order)*P(:,:,k-j)*OrderY(j,N,order)';
        end
    else
        for j=2:k-1
            r=r+OrderY(j,N,order)*P(:,:,k-j)*OrderY(j,N,order)';
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
