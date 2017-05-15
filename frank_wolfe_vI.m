function [Z,E,D]=frank_wolfe_vI(X, niter)
[m,n]=size(X);
S = 1:m*n;
xij = X(S);

%S %index for know values
%xij= %known values
f=@(zij) 1/2*sum((zij-xij).^2); 
fGrad=zeros(m,n);
Z=zeros(m,n);
delta=10*trace(X'*X); %nuclear normf
L=1; %Lipschitz constant
dm=2*delta;
L_bar=L;
dm_bar=dm;
gamma1=0.8; %0<=gamma1<=gamma2<=1
gamma2=0.8;

fGrad(S)=Z(S)-xij;
[u,~,v]=svds(fGrad,1);

Z=-delta*u*v';
U=u;V=-v;D=delta;
B=max(f(0)+trace(fGrad'*Z),0);

cnt=0;

for k=1:niter
    fGrad(S)=Z(S)-xij;
    
    [u,~,v]=svds(fGrad,1);
    
    %In-face
    r=length(D);
    if r==min(m,n) %int(B)
        Z_hat=delta*u*v';
        %binary search for alpha_stop
        alpha_stop=binary_search(U,D,V,u,v,delta);
        
    else
        G=1/2*(V'*fGrad'*U+U'*fGrad*V);
        [ug,~]=eigs(G,1);
        ug=ug/norm(ug);
        Z_hat=U*delta*(ug*ug')*V';
        alpha_stop=(delta*ug'*D^-1*ug-1)^-1;
    end
    d=Z-Z_hat;
    
    if trace(fGrad'*d)<0 | 1==1
        
        Z_B=Z+alpha_stop*d;

        beta=min(-trace(fGrad'*d)/(L_bar*norm(d,'fro')^2),alpha_stop);
        Z_A=Z+beta*d; %beta [0,alpha_stop]

        if 1/(f(Z_B(S))-B)>=1/(f(Z(S))-B)+gamma1/(2*L_bar*dm_bar^2)
            Z=Z_B;
            if r==min(m,n)
                [U,D,V]=svd_update(U,(1+alpha_stop)*D,V,-alpha_stop*delta*u,v);
            else
                [R,D,R]=svd_thin((1+alpha_stop)*D+alpha_stop*delta*(ug*ug'));
                U=U*R;
                V=V*R;
            end

        elseif 1/(f(Z_A(S))-B)>=1/(f(Z(S))-B)+gamma2/(2*L_bar*dm_bar^2)
            Z=Z_A;
            if r==min(m,n)
                [U,D,V]=svd_update(U,(1+beta)*D,V,-beta*delta*u,v);
            else
                [R,D,R]=svd_thin((1+beta)*D+beta*delta*(ug*ug'));
                U=U*R;
                V=V*R;
            end

        else
            cnt=cnt+1;
            Z_tilda=-delta*u*v';
            Bw=f(Z(S))+trace(fGrad'*(Z_tilda-Z));
            B=max(B,Bw); 
            alpha=min(trace(fGrad'*(Z-Z_tilda))/(L_bar*norm(Z-Z_tilda,'fro')^2),1);
            Z=Z+alpha*(Z_tilda-Z);
            [U,D,V]=svd_update(U,(1-alpha)*D,V,-alpha*delta*u,v);

        end
    
    else
        cnt=cnt+1;
        Z_tilda=-delta*u*v';
        Bw=f(Z(S))+trace(fGrad'*(Z_tilda-Z));
        B=max(B,Bw); 
        alpha=min(trace(fGrad'*(Z-Z_tilda))/(L_bar*norm(Z-Z_tilda,'fro')^2),1);
        Z=Z+alpha*(Z_tilda-Z);
        [U,D,V]=svd_update(U,(1-alpha)*D,V,-alpha*delta*u,v);
    end
    
end

% rank_est=rank(Z) %estimation rank
err_est=sum((X(:)-Z(:)).^2)/(size(X,1)*size(X,2)) %error on all values
err_est_percentage=sqrt(sum((X(S)-Z(S)).^2)/(norm(X(S))^2))*100 %error on known valuse

E=X-Z;
end
function [u,d,v]=svd_thin(x)

[u,d,v]=svd(x);
idx=find(diag(d)>1e-6);
u=u(:,idx);
v=v(:,idx);
d=d(idx,idx);

end

function [u,s,v]=svd_update(U,S,V,a,b)

%USV'+ab'
r=length(S);
m=U'*a;p=a-U*m;
p_norm=norm(p);P=p/p_norm;
n=V'*b;q=b-V*n;
q_norm=norm(q);Q=q/q_norm;
K=[S,zeros(r,1);zeros(1,r),0]+[m;p_norm]*[n;q_norm]';
[uk,s,vk]=svd_thin(K);
u=[U P]*uk;
v=[V Q]*vk;

end
   
function alpha_stop=binary_search(U,D,V,u,v,delta)

    alpha_stop=0;alpha_max=1;
    for i=1:10
        alpha=(alpha_max-alpha_stop)/2;
        [~,s,~]=svd_update(U,(1+alpha)*D,V,-alpha*delta*u,v);
        if sum(diag(s))<=delta
            alpha_stop=alpha;
        else
            alpha_max=alpha;
        end
    end

end