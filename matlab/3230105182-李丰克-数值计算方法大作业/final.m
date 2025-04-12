strinput="K4Fe(CN)6,KMnO4,H2SO4";
stroutput="CO2,KNO3,H2O,K2SO4,MnSO4,Fe2(SO4)3";
cheminput=strsplit(strinput,",");
chemoutput=strsplit(stroutput,",");
function [element,element_num]=unpack_element(str)
    element_temp=[];
    num_temp=[];
    element=[];
    element_num=[];
    temp1=strrep(str,"·","");
    temp=char(temp1);
    for j=1:length(temp)%遍历每一个字符
        if (double(temp(j))>=65)&&(double(temp(j))<=90)%只有这个字符是大写字母时才会停下
            if j==length(temp)%判断是否是最后一个字符
                element_temp=[element_temp,string(temp(j:j))];
                num_temp=[num_temp,1];
            else
                if (double(temp(j+1))>=97)&&(double(temp(j+1))<=122)%如果下一个是小写字母
                    element_temp=[element_temp,string(temp(j:j+1))];%这两个字母组合为元素名称
                    if j+1==length(temp)%判断接下来还有没有字符
                        num_temp=[num_temp,1];
                    else
                        if (double(temp(j+2))>=65)&&(double(temp(j+2))<=90)%如果再下一个是大写字母
                            num_temp=[num_temp,1];%此元素下标为1
                        end
                        if (double(temp(j+2))>=48)&&(double(temp(j+2))<=57)%如果再下一个是数字
                              k=j+2;
                              num=[];
                              while k<=length(temp)&&(double(temp(k))>=48)&&(double(temp(k))<=57)%把下面所有数字位组合成字符串
                                num=[num,temp(k)];
                                k=k+1;
                              end
                              n=str2num(num);%字符串转数字即为下标
                              num_temp=[num_temp,n];
                        end
                    end
                end
    
                if (double(temp(j+1))>=65)&&(double(temp(j+1))<=90)%如果下一个是大写字母
                    element_temp=[element_temp,string(temp(j:j))];
                    num_temp=[num_temp,1];
                end
                if (double(temp(j+1))>=48)&&(double(temp(j+1))<=57)%如果下一个是数字
                    element_temp=[element_temp,string(temp(j:j))];
                    k=j+1;
                    num=[];
                    while k<=length(temp)&&((double(temp(k))>=48)&&(double(temp(k))<=57))%把下面所有数字位组合成字符串
                            num=[num,temp(k)];
                            k=k+1;
                    end
                    n=str2num(num);%字符串转数字即为下标
                    num_temp=[num_temp,n];
                end
            end
        end
    end
    if ~isempty(element_temp)
    element=[element,element_temp(1)];
    element_num=[element_num,num_temp(1)];
    if length(element_temp)>1
        for i=2:length(element_temp)%把相同元素合并起来
            flag=0;
            for j=1:length(element)%用element_temp里的值遍历element里的值，有就相加数，没有就添上
                if element_temp(i)==element(j)
                    element_num(j)=element_num(j)+num_temp(i);
                    flag=1;
                    break;
                end
            end
            if flag==0
                element=[element,element_temp(i)];
                element_num=[element_num,num_temp(i)];
            end
        end
    end
    end
end

function [element,element_num]=unpack_complex_chem(str)
    left=[];
    right=[];
    temp1=strrep(str,"·","");
    temp=char(temp1);
    element_with=[];
    element_with_temp=[];%存放括号里面的字符串
    num_with_temp=[];%存放每个括号的下标
    num_with_length=[];%存放每个括号下标数字长度，没有则为0
    for i=1:length(temp)%首先找到所有括号的位置
        if temp(i)=='('
            left=[left,i];
        end
        if temp(i)==')'
            right=[right,i];
        end
    end
    if ~isempty(left)%如果存在括号
        for i=1:length(left)
            element_with_temp=[element_with_temp,string(temp(left(i)+1:right(i)-1))];%把括号内部的字符串提取出来
            if right(i)==length(temp)
                num_with_temp=[num_with_temp,1];
                num_with_length=[num_with_length,0];
            else 
                k=right(i)+1;
                num=[];
                while k<=length(temp)&&((double(temp(k))>=48)&&(double(temp(k))<=57))%对括号后的数字进行处理
                    num=[num,temp(k)];
                    k=k+1;
                end
                if isempty(num)%如果没有数字，下标1，长度0
                    num_with_temp=[num_with_temp,1];
                    num_with_length=[num_with_length,length(num)];
                else
                    n=str2num(num);%字符串转数字即为下标
                    num_with_temp=[num_with_temp,n];
                    num_with_length=[num_with_length,length(num)];
                end
            end
        end
        str_temp=temp;
        for i=1:length(left)
            substring=temp(left(i):(right(i)+num_with_length(i)));
            str_temp=strrep(str_temp,substring,"");
        end
        element_with=[];
        element_with_num=[];
        for i=1:length(element_with_temp)
            [a,b]=unpack_element(element_with_temp(i));
            b=b.*num_with_temp(i);
            element_with=[element_with,a];
            element_with_num=[element_with_num,b];
        end
        [a,b]=unpack_element(str_temp);
        element_with=[element_with,a];
        element_with_num=[element_with_num,b];
    end
    if isempty(left)
        [element_with,element_with_num]=unpack_element(temp);%如果没有括号
        element=[];
        element_num=[];
    end
    if ~isempty(element_with)
    element=[];
    element_num=[];
    element=[element,element_with(1)];
    element_num=[element_num,element_with_num(1)];
    if length(element_with)>1
        for i=2:length(element_with)%把相同元素合并起来
            flag=0;
            for j=1:length(element)%用element_temp里的值遍历element里的值，有就相加数，没有就添上
                if element_with(i)==element(j)
                    element_num(j)=element_num(j)+element_with_num(i);
                    flag=1;
                    break;
                end
            end
            if flag==0
                element=[element,element_with(i)];
                element_num=[element_num,element_with_num(i)];
            end
        end
    end
    end
end

ele_matrix_input={};%定义元胞数组
elenum_matrix_input={};
for i=1:length(cheminput)
    [a,b]=unpack_complex_chem(cheminput(i));
    ele_matrix_input{end+1}=a;
    elenum_matrix_input{end+1}=b;
end
ele_matrix_output={};%定义元胞数组
elenum_matrix_output={};
for i=1:length(chemoutput)
    [a,b]=unpack_complex_chem(chemoutput(i));
    ele_matrix_output{end+1}=a;
    elenum_matrix_output{end+1}=b;
end
disp(ele_matrix_input);
disp(elenum_matrix_input);
disp(ele_matrix_output);
disp(elenum_matrix_output);

a_marix=[];
a_marix=[a_marix,ele_matrix_input{1}(1)];
for i=1:length(ele_matrix_input)
    for j=1:length(ele_matrix_input{i})
        flag=0;
        for k=1:length(a_marix)
            if a_marix(k)==ele_matrix_input{i}(j)
                flag=1;
                break;
            end
        end
        if flag==0
            a_marix=[a_marix,ele_matrix_input{i}(j)];
        end
    end
end
disp(a_marix);
a=length(a_marix);%参与反应的元素个数，根据元素守恒列方程式，为方程式个数，即矩阵行数
b=length(cheminput)+length(chemoutput);%化学式个数，即未知量个数，为矩阵的列数
disp("列数");
disp(b);
disp("行数");
disp(a);
A=zeros(a,b);
for i=1:length(a_marix)
    for j=1:length(cheminput)
        flag=0;
        for k=1:length(ele_matrix_input{j})
            if ele_matrix_input{j}(k)==a_marix(i)
                A(i,j)=elenum_matrix_input{j}(k);
                flag=1;
                break;
            end
        end
        if flag==0
            A(i,j)=0;
        end
    end
    for j=length(cheminput)+1:b
        flag=0;
        l=j-length(cheminput);
        for k=1:length(ele_matrix_output{l})
            if ele_matrix_output{l}(k)==a_marix(i)
                A(i,j)=-elenum_matrix_output{l}(k);
                flag=1;
                break;
            end
        end
        if flag==0
            A(i,j)=0;
        end
    end
end
disp(A);
disp("秩数")
disp(rank(A));
delete_num=a-rank(A);
A_temp=A;
A_temp(1:delete_num,:)=[];
while rank(A_temp)~=rank(A)
    A_temp=A;
    randomint=randperm(a,delete_num);
    A_temp(randomint,:)=[];
end
give_num=b-rank(A);
B=[];
for i=1:rank(A)
    sum=0;
    for j=1:give_num
        sum=sum+A_temp(i,j);
    end
    B=[B,-sum];
end
A_temp_temp=A_temp;
A_temp(:,1:give_num)=[];


disp("系数矩阵")
disp(A_temp);
disp(B);
%LU分解法
function X=LU_solve(A,b)
    [m,n]=size(A);
    L=eye(n);
    U=A;
    %LU分解
    for k=1:n
        for i=k+1:n
            L(i,k)=U(i,k)/U(k,k);
            U(i,k:n)=U(i,k:n)-L(i,k)*U(k,k:n);
        end
    end
    %前向代入求解Ly=b
    y=zeros(n,1);
    for i=1:n
        y(i)=b(i)-L(i,1:i-1)*y(1:i-1);
    end
    %后向代入求解Ux=y
    X=zeros(n,1);
    for i=n:-1:1
        X(i)=(y(i)-U(i,i+1:n)*X(i+1:n))/U(i,i);
    end
end
%高斯消元法
function x=Gauss(a,b)
    [m,n]=size(a);
    for k=1:n-1
        for i=k+1:n
            factor=a(i,k)/a(k,k);
            for j=k+1:n
                a(i,j)=a(i,j)-factor*a(k,j);
            end
            b(i)=b(i)-factor*b(k);
        end
    end
    x(n)=b(n)/a(n,n);
    for i=n-1:-1:1
        sum=b(i);
       for j=i+1:n
           sum=sum-a(i,j)*x(j);
       end
       x(i)=sum/a(i,i);
    end
end
%雅各比迭代法
function [x,k]=jacobi(a,b,x0,epsilon)
    n=size(a,1);
    x=zeros(1,n);
    flag=0;
    k=0;
    while 1
    k=k+1;
    for i=1:n
        sum=0;
        for j=1:n
            if i~=j
                sum=sum+a(i,j)*x0(j);
            end
        end
        x(i)=(b(i)-sum)/a(i,i);
        if abs((x(i)-x0(i))/x(i))<epsilon
            flag=1;
        end    
        x0(i)=x(i);
    end    
    if flag==1
        break;
    end  
    if k==1000
        break;
    end
    end
end
%高斯-赛德尔迭代法
function [x,k]=GS(a,b,x0,epsilon)
    n=size(a,1);
    x=zeros(1,n);
    flag=0;
    k=0;
    while 1
    k=k+1;
    for i=1:n
        sum1=0;
        sum2=0;
        if i~=1
            for j=1:i-1
                sum1=sum1+a(i,j)*x(j);
            end
        end
        for j=i+1:n
            sum2=sum2+a(i,j)*x0(j);
        end
        x(i)=(b(i)-sum1-sum2)/a(i,i);
        if abs((x(i)-x0(i))/x(i))<epsilon
            flag=1;
        end    
        x0(i)=x(i);
    end    
    if flag==1
        break;
    end   
    if k==1000
        break;
    end
end
end
function x=ColMajorGauss(a,b)
    [m,n]=size(a);
    for k=1:n-1
        col=a(k:n,k);%取出第k列
        col=abs(col);%取绝对值
        [max_value,max_index]=max(col);%获得最大值的索引并转为真实索引
        if max_index+k-1 ~=k%对比是否为当前行，如果不是则换行
            a([k,max_index+k-1], :) = a([max_index+k-1,k],:);
            b([k,max_index+k-1])=b([max_index+k-1,k]);%常数项也要交换
        end
        for i=k+1:n%选好主元后再进行高斯消元
            factor=a(i,k)/a(k,k);
            for j=k:n
                a(i,j)=a(i,j)-factor*a(k,j);
            end
            b(i)=b(i)-factor*b(k);
        end
    end
    x(n)=b(n)/a(n,n);
    for i=n-1:-1:1
        sum=b(i);
       for j=i+1:n
           sum=sum-a(i,j)*x(j);
       end
       x(i)=sum/a(i,i);
    end
end
x0=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
format rat;
X=ColMajorGauss(A_temp,B);

disp(X);
D_array=[];
for i=1:length(X)
    [I,D]=numden(sym(X(i)));
    D=eval(D);
    D_array=[D_array,D];
end
disp("分数矩阵");
disp(D_array);
tool=D_array(1);
for i=2:length(X)
    if ~isinteger(tool/D_array(i))
        tool=tool*D_array(i)/gcd(tool,D_array(i));
    end
end
X=X*tool;
disp(X);
disp("最小公倍数为")
disp(tool);
fprintf("%d",tool);
fprintf(cheminput(1));
for i=2:length(cheminput)
    fprintf("+");
    fprintf("%d",int32(X(i-1)));
    fprintf(cheminput(i));
end
fprintf("=\n");
fprintf("%d",int32(X(length(cheminput))));
fprintf(chemoutput(1));
for i=2:length(chemoutput)
    fprintf("+");
    fprintf("%d",int32(X(length(cheminput)+i-1)));
    fprintf(chemoutput(i));
end
    