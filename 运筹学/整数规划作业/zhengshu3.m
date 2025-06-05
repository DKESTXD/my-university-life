a=[3,5,7,4,2,1,9,5,4,2];%重量
c=[4,6,7,3,1,1,10,6,3,3];%价值
b=25;%最大重量
best_x=[];% 最优解
best_value=0;%最优值
best_weight=0;%最优值对应的总重量
lb=0;%隐枚举的下界
for i=1:1024%共有2^10种可能
    tool=dec2bin(i-1,10);%用二进制表示组合
    sum_weight=0;%总重量
    sum_value=0;%总价值
    for j=1:length(tool)%先求总价值   
        sum_value=sum_value+(tool(j)-48)*c(j);
    end
    if sum_value>lb%看是否满足隐枚举条件
        for j=1:length(tool)
            sum_weight=sum_weight+(tool(j)-48)*a(j);
        end
        if sum_weight<=b%再判断是否满足重量约束
            lb=sum_value;%更新下界
            best_x=tool;
            best_value=sum_value;
            best_weight=sum_weight;
        end
    end
end
disp(best_x);
disp(best_value);
disp(best_weight);
