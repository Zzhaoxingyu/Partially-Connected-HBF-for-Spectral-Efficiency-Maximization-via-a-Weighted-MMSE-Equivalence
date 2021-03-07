function a=annealing_algorithm(a1,a2,b1,b2)
%模拟退火算法实现搜索函数最大值
%Anneal.m
N=20;%粒子数量
temp=20;%初始温度
T=200;%迭代次数
k=0.1;%温度位移系数
kt=0.05;%温度概率系数
de=0.99;%温度降低速率
minx=0;
maxx=2*pi;%区间
location=2*pi*rand(1,N);%粒子初始位置
present_value=equation(location,a1,a2,b1,b2);%粒子当前解
%---------------------------
for t=1:T %
    dx_av=k*temp;%当前温度下粒子平均移动距离
    probability=exp(-1/(kt*temp));
    %disp(probability);
    temp=temp*de;%温度变化
    for p=1:N
        dx=0.5*dx_av*randn+dx_av;%以平均移动距离为中心正态分布，
        if rand>0.5 %0.5的概率为-
            dx=-dx;
        end
        local=location(p)+dx;
        if (local<maxx)&&(local>minx)%判断是否越界
            local_value=equation(local,a1,a2,b1,b2);
            if local_value>present_value(p)
                location(p)=local;
                present_value(p)=local_value;
            else if rand<probability
                    location(p)=local;
                    present_value(p)=local_value;
                end
            end
        end
    end
end
a = location(N);
end
