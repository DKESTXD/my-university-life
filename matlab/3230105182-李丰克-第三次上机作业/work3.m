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

a_test=[2 3 -1;1 -1 4;3 2 1];
b_test=[5 11 10];
x_test=Gauss(a_test,b_test);
for i=1:3
    fprintf("%.16f\n",x_test(i));
end


a=[1 2 4 8 16;1 3 8 27 81;1 4 16 64 256;1 5 25 125 625;1 6 36 216 1295];
b=[16 81 256 625 1296];
aplusb=[a,b'];
disp(rank(aplusb));
disp(rank(a));
c=Gauss(a,b);

[m,n]=size(a);

for i=1:n
    sum=0;
    for j=1:n
        sum=sum+a(i,j)*c(j);
    end
    fprintf("%.16f\n",sum);
end
fprintf("\n");

for i=1:n
    fprintf("%.16f\n",c(i));
end
fprintf("\n");

x=linsolve(a,b');
for i=1:n 
    sum=0;
    for j=1:n
        sum=sum+a(i,j)*x(j);
    end
    fprintf("%.16f\n",sum);
end
fprintf("\n");
for i=1:n
    fprintf("%.16f\n",x(i));
end
fprintf("\n");
c_4=round(c,4); 
for i=1:n 
    sum=0;
    for j=1:n
        sum=sum+a(i,j)*c_4(j);
    end
    fprintf("%.16f\n",sum);
end
fprintf("\n");
a_2=[1 2 4 8 16;1 3 8 27 81;1 4 16 64 256;1 5 25 125 625;1 6 36 216 1296];
for i=1:n 
    sum=0;
    for j=1:n
        sum=sum+a_2(i,j)*c(j);
    end
    fprintf("%.16f\n",sum);
end