%利润矩阵
r=[4 3 4 4 5 6;
   3 4 5 3 4 5;
   5 3 4 5 5 4;
   3 3 4 4 6 6;
   3 3 3 4 5 7];
%布料单价
a=[6 6 7 8 9 10];

num_shop=5;%车间数
num_cloth=6;%布的种类
%目标函数
r=r';
c=-r(:);%取负数转化为最小化问题
num_vars=num_shop*num_cloth;%决策变量数
%加工能力约束，一个车间不超过10000
A=zeros(num_shop,num_vars);
b=zeros(num_shop,1);
for i = 1:num_shop
    A(i,(i-1)*num_material+1:i*num_material)=1;
    b(i)=10000;
end
%资金能力约束
a_expanded=zeros(1, num_vars);
for i=1:num_shop
    a_expanded((i-1)*num_material+1:i*num_material)=a;
end
A = [A;a_expanded];
b = [b;400000];
%生产要求约束
A_low=-eye(num_vars);
b_low=-ones(num_vars, 1)*1000;
A=[A; A_low];
b=[b; b_low];
%无等式约束
%初始化分支定界参数
best_obj=-inf;%最优目标函数值
best_x=zeros(num_vars,1);%最优决策变量
nodes_explored=0;%遍历节点数
lb=ones(num_vars,1)*1000;%下界
ub=ones(num_vars,1)*10000;%上界
%创建根节点
root_node.x=[];
root_node.lb=lb;
root_node.ub=ub;
root_node.is_integer=false;
root_node.obj=inf;

%创建节点队列
node_queue = {root_node};
%检查是否为整数解的函数,返回是否为全整数解以及非整数解的变量下标
function [is_integer,fractional_idx]=check_integer(x)
        tol=1e-8;
        fractional_idx=find(abs(x-round(x))>tol);
        is_integer=isempty(fractional_idx);
end
%分支定界法
while ~isempty(node_queue)
    %每次都遍历队列的下一个节点并且遍历完后出队
    current_node=node_queue{1};
    node_queue(1)=[];
    nodes_explored=nodes_explored+1;
    
    %求解当前节点的线性规划问题
    [x,obj,status]=linprog(c,A,b,[],[],current_node.lb,current_node.ub);
    
    %如果问题不可行，跳过此节点
    if status~=1
        continue;
    end
    x_mar=[];
    for i=1:5
        x_mar=[x_mar,x(6*i-5:6*i)];
    end
    x_mar=x_mar';
    disp(x_mar);
    %如果目标值比当前最优值还差，剪枝
    if obj>=-best_obj
        continue;
    end 
    %检查解是否为整数解
    [is_integer,fractional_idx]=check_integer(x);
    if is_integer
        %是整数解，更新最优解
        if -obj>best_obj
            best_obj=-obj;
            best_x=x;
        end
    else
        %不是整数解，进行分支
        branch_var=fractional_idx(1);  %选择第一个非整数变量进行分支
        
        %创建左子节点（添加上限约束）
        left_node=current_node;
        left_node.ub(branch_var)=floor(x(branch_var));
        left_node.is_integer=false;
        left_node.obj=inf;
        
        %创建右子节点（添加下限约束）
        right_node = current_node;
        right_node.lb(branch_var)=ceil(x(branch_var));
        right_node.is_integer=false;
        right_node.obj=inf;
        
        %将子节点加入队列
        node_queue=[node_queue, {left_node}, {right_node}];
    end
end

x_matrix=[];
    for i=1:5
        x_matrix=[x_matrix,best_x(6*i-5:6*i)];
    end
x_matrix=x_matrix';
fprintf('每个车间每种布料的分配量:\n');
disp(x_matrix);
fprintf('最大总利润:%.2f元\n',best_profit);
fprintf('遍历的节点数:%d\n',nodes_explored);
sum=0;
r=r';
for i=1:5
    for j=1:6
        sum=sum+x_matrix(i,j)*r(i,j);
    end
end
disp(sum);
