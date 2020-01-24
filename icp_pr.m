function [RR,TT,data] = icp_pr(model,data)
% ICP (modified iterative closest point) algorithm

%% General variables
thd=1e-5;  % threshold to stop icp iterations
minIter=5;   % min number of icp iterations
maxIter=6;  % max number of icp iterations (although 100 is reasonable)

%% Specific Variables
% Size of model points and data points
m=size(model,1);
M=size(model,2);
N=size(data,2);

% Create closest point search structure
DT=[];
indx=zeros(N,1);
vi=ones(N,1);

% Initiate transformation
RR=eye(m);
TT=zeros(m,1);

%% Start the ICP algorithm
res=9e99;
for iter=1:maxIter
    oldres=res;
    % Find closest model points to data points
    for i=1:N
        minval=9e99;
        for j=1:M
            val=norm(data(:,i)-model(:,j));
            if val<minval
                minval=val;
                vi(i)=j;
                indx(i)=val;
            end
        end
    end
    
    res=mean(indx.^2);
    med=mean(data,2);
    mem=mean(model(:,vi),2);
    C=data*model(:,vi)'-(N*med)*mem';
    [U,~,V]=svd(C);
    Rot_i=V*U'; % New rotation
    if det(Rot_i)<0
        V(:,end)=-V(:,end);
        Rot_i=V*U';
    end
    Trans_i=mem-Rot_i*med; %New Translation
    
    % Apply the new transformation
    data=Rot_i*data;                       
    for i=1:m
        data(i,:)=data(i,:)+Trans_i(i);      
    end
    
    %Update the transformation
    RR=Rot_i*RR;                            
    TT=Rot_i*TT+Trans_i;                        
    
    if iter >= minIter
        if abs(oldres-res) < thd
            break
        end
    end
    
end

