//Date :13/01/2022
//Aim : To solve a second order ordinary differential equation with Dirichlet boundary conditions using shooting method
clc
clear
function Dy = f(x,y)
    Dy(1)=y(2)
    Dy(2)= (%pi^2)*y(1)-2*(%pi^2)*sin(%pi*x)
endfunction
x0=0 
y0=0 //initial value of y(0)
y_prime(1)=input("First approximation : ")
y_prime(2)=input("Second Approximation : ")
yn=input("Enter the given Boundary value of y for x(n) : ")
x=[0:0.01:1]//range of the independent variable (from initial to the boundary)
h=x(2)-x(1)
n=length(x)//corresponding values in the range of independent variable
Tol=input("Enter the tolerance : ")
exec('C:\Users\MANIDIPA BANERJEE\Desktop\MP III Scilab SEM - IV\Rk_4 Func.sce', -1)
Y(1,:)=RKfourth(x0,y0,y_prime(1),x,f,n,h)//solution with initial value and first approximation
Y(2,:)=RKfourth(x0,y0,y_prime(2),x,f,n,h)//solution with initial value and second approximation
err(1)=abs(Y(1,n)-yn)//error against first approximation
err(2)=abs(Y(2,n)-yn)//error against second approximation
i=2 //will be taking third approximation by taking i = 2 
if (err(1)<Tol||err(2)<Tol) then
    break;
else
    while(err(i)>Tol)
        y_prime(i+1)=y_prime(i)-(((y_prime(i)-y_prime(i-1))*(Y(i,n)-yn))/((Y(i,n)-yn)-(Y(i-1,n)-yn))) //secant method to upgrade the value of approximation // here for first iteration i would be second approximation and i-1 would be first approximation 
        Y(i+1,:)=RKfourth(x0,y0,y_prime(i+1),x,f,n,h)//solution using upgraded (i+1) approximation
        err(i+1)=abs(Y(i+1,n)-yn)//for next upgraded approximation
        i=i+1
    end
end
disp("approximation : ",y_prime)
disp("solution : ",Y(i,:))
//Plotting with true solution
plot(x,Y(i,:))
R=sin(%pi*x)//original solution
plot(x,R,"*g")
title("Solving Boundary valued Problem using Shooting Method")
title color Red
title fontsize 4
xlabel("X----->")
xlabel color magenta fontsize 4
ylabel("Y(X)---->")
ylabel color magenta fontsize 4
legend(["Using Shooting Method";"Exact Solution"])
xgrid


