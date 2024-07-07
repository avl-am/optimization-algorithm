% This function draw the benchmark functions

% 定义绘图函数，函数申明
function func_plot(func_name)
% 调用了名为 Get_F 的函数，并将返回的结果赋值给变量 lb、ub、dim 和 fobj
[lb,ub,dim,fobj]=Get_F(func_name);

% switch 语句；根据变量 func_name 的值进行分支处理

% 定义了x与y的取值范围；[开始：步长：结束]
switch func_name 
    case 'F1' 
        x=-100:2:100; y=x; %[-100,100] 
    case 'F2' 
        x=-100:2:100; y=x; %[-10,10]
    case 'F3' 
        x=-100:2:100; y=x; %[-100,100]
    case 'F4' 
        x=-100:2:100; y=x; %[-100,100]
    case 'F5' 
        x=-200:2:200; y=x; %[-5,5]
    case 'F6' 
        x=-100:2:100; y=x; %[-100,100]
    case 'F7' 
        x=-1:0.03:1;  y=x  %[-1,1]
    case 'F8' 
        x=-500:10:500;y=x; %[-500,500]
    case 'F9' 
        x=-5:0.1:5;   y=x; %[-5,5]    
    case 'F10' 
        x=-20:0.5:20; y=x;%[-500,500]
    case 'F11' 
        x=-500:10:500; y=x;%[-0.5,0.5]
    case 'F12' 
        x=-10:0.1:10; y=x;%[-pi,pi]
    case 'F13' 
        x=-5:0.08:5; y=x;%[-3,1]
    case 'F14' 
        x=-100:2:100; y=x;%[-100,100]
    case 'F15' 
        x=-5:0.1:5; y=x;%[-5,5]
    case 'F16' 
        x=-1:0.01:1; y=x;%[-5,5]
    case 'F17' 
        x=-5:0.1:5; y=x;%[-5,5]
    case 'F18' 
        x=-5:0.06:5; y=x;%[-5,5]
    case 'F19' 
        x=-5:0.1:5; y=x;%[-5,5]
    case 'F20' 
        x=-5:0.1:5; y=x;%[-5,5]        
    case 'F21' 
        x=-5:0.1:5; y=x;%[-5,5]
    case 'F22' 
        x=-5:0.1:5; y=x;%[-5,5]     
    case 'F23' 
        x=-5:0.1:5; y=x;%[-5,5]  
end    

%计算函数的数值
L=length(x);
%定义空的矩阵用来存储函数的数值
f=[];

% 使用两个嵌套的循环遍历变量 x 和 y 的所有元素。根据 func_name 的取值来计算对应函数的数值，并赋值给矩阵 f 的对应位置
for i=1:L
    for j=1:L
        if strcmp(func_name,'F15')==0 && strcmp(func_name,'F19')==0 && strcmp(func_name,'F20')==0 && strcmp(func_name,'F21')==0 && strcmp(func_name,'F22')==0 && strcmp(func_name,'F23')==0
            f(i,j)=fobj([x(i),y(j)]);
        end

        if strcmp(func_name,'F15')==1
            f(i,j)=fobj([x(i),y(j),0,0]);
        end

        if strcmp(func_name,'F19')==1
            f(i,j)=fobj([x(i),y(j),0]);
        end

        if strcmp(func_name,'F20')==1
            f(i,j)=fobj([x(i),y(j),0,0,0,0]);
        end   

        if strcmp(func_name,'F21')==1 || strcmp(func_name,'F22')==1 ||strcmp(func_name,'F23')==1
            f(i,j)=fobj([x(i),y(j),0,0]);
        end          
    end
end
% surfc 函数来绘制函数的三维图像；'LineStyle','none' 参数用于去除曲面的网格线。
surfc(x,y,f,'LineStyle','none');

end

