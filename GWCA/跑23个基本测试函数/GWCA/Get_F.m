% 函数细节

% 函数声明
function [LB,UB,Dim,F_obj] = Get_F(F)
    switch F
        case 'F1'
            F_obj = @F1;
            LB=-100;
            UB=100;
            Dim =10;    
        case 'F2'
            F_obj = @F2;
            LB=-10;
            UB=10;
            Dim = 10;        
        case 'F3'
            F_obj = @F3;
            LB=-100;
            UB=100;
            Dim = 10;        
        case 'F4'
            F_obj = @F4;
            LB=-100;
            UB=100;
            Dim = 10;        
        case 'F5'
            F_obj = @F5;
            LB=-30;
            UB=30;
            Dim = 10;        
        case 'F6'
            F_obj = @F6;
            LB=-100;
            UB=100;
            Dim = 10;        
        case 'F7'
            F_obj = @F7;
            LB=-1.28;
            UB=1.28;
            Dim = 10;        
        case 'F8'
            F_obj = @F8;
            LB=-500;
            UB=500;
            Dim = 10;        
        case 'F9'
            F_obj = @F9;
            LB=-5.12;
            UB=5.12;
            Dim = 10;        
        case 'F10'
            F_obj = @F10;
            LB=-32;
            UB=32;
            Dim = 10;        
        case 'F11'
            F_obj = @F11;
            LB=-600;
            UB=600;
            Dim = 10;        
        case 'F12'
            F_obj = @F12;
            LB=-50;
            UB=50;
            Dim = 10;        
        case 'F13'
            F_obj = @F13;
            LB=-50;
            UB=50;
            Dim = 10;        
        case 'F14'
            F_obj = @F14;
            LB=-65.536;
            UB=65.536;
            Dim=2;
        case 'F15'
            F_obj = @F15;
            LB=-5;
            UB=5;
            Dim=4;
        case 'F16'
            F_obj = @F16;
            LB=-5;
            UB=5;
            Dim=2;
        case 'F17'
            F_obj = @F17;
            LB=[-5,0];
            UB=[10,15];
            Dim=2;
        case 'F18'
            F_obj = @F18;
            LB=-2;
            UB=2;
            Dim=2;      
        case 'F19'
            F_obj = @F19;
            LB=0;
            UB=1;
            Dim=3;   
        case 'F20'
            F_obj = @F20;
            LB=0;
            UB=1;
            Dim=6;     
            
        case 'F21'
            F_obj = @F21;
            LB=0;
            UB=10;
            Dim=4;    
            
        case 'F22'
            F_obj = @F22;
            LB=0;
            UB=10;
            Dim=4;  
        case 'F23'
            F_obj = @F23;
            LB=0;
            UB=10;
            Dim=4;
    end
end

%% 函数

% F1

function f= F1(x)
    f = sum(x.^2);
end

% F2

function f = F2(x)
    f = sum(abs(x))+prod(abs(x));
end

% F3

function f = F3(x)
    dim=size(x,2);
    f = 0;
    for i=1:dim
        f = f+sum(x(1:i))^2;
    end
end

% F4
function f = F4(x)
    f = max(abs(x));
end

% F5
function f = F5(x)
    dim=size(x,2);
    f=sum(100*(x(2:dim)-(x(1:dim-1).^2)).^2+(x(1:dim-1)-1).^2);
end

% F6
function f = F6(x)
    f=sum(abs((x+.5)).^2);
end

% F7
function f = F7(x)
    dim=size(x,2);
    f=sum([1:dim].*(x.^4))+rand;
end

% F8
function f = F8(x)
    f=sum(-x.*sin(sqrt(abs(x))));
end

% F9
function f = F9(x)
    dim=size(x,2);
    f=sum(x.^2-10*cos(2*pi.*x))+10*dim;
end

% F10
function f = F10(x)
    dim=size(x,2);
    f=-20*exp(-.2*sqrt(sum(x.^2)/dim))-exp(sum(cos(2*pi.*x))/dim)+20+exp(1);
end

% F11
function f = F11(x)
    dim=size(x,2);
    f=sum(x.^2)/4000-prod(cos(x./sqrt([1:dim])))+1;
end

% F12
function f= F12(x)
    dim=size(x,2);
    f = (pi/dim)*(10*((sin(pi*(1+(x(1)+1)/4)))^2)+sum((((x(1:dim-1)+1)./4).^2).*...
    (1+10.*((sin(pi.*(1+(x(2:dim)+1)./4)))).^2))+((x(dim)+1)/4)^2)+sum(Ufun(x,10,100,4));
end

% F13
function f = F13(x)
    dim=size(x,2);
    f=.1*((sin(3*pi*x(1)))^2+sum((x(1:dim-1)-1).^2.*(1+(sin(3.*pi.*x(2:dim))).^2))+...
    ((x(dim)-1)^2)*(1+(sin(2*pi*x(dim)))^2))+sum(Ufun(x,5,100,4));
end

% F14

function f = F14(x)
    aS=[-32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32;,...
    -32 -32 -32 -32 -32 -16 -16 -16 -16 -16 0 0 0 0 0 16 16 16 16 16 32 32 32 32 32];
    for j=1:25
        bS(j)=sum((x'-aS(:,j)).^6);
    end
    f=(1/500+sum(1./([1:25]+bS))).^(-1);
end

% F15
function f = F15(x)
    aK=[.1957 .1947 .1735 .16 .0844 .0627 .0456 .0342 .0323 .0235 .0246];
    bK=[.25 .5 1 2 4 6 8 10 12 14 16];bK=1./bK;
    f=sum((aK-((x(1).*(bK.^2+x(2).*bK))./(bK.^2+x(3).*bK+x(4)))).^2);
end

% F16
function f = F16(x)
    f=4*(x(1)^2)-2.1*(x(1)^4)+(x(1)^6)/3+x(1)*x(2)-4*(x(2)^2)+4*(x(2)^4);
end

% F17
function f = F17(x)
    f=(x(2)-(x(1)^2)*5.1/(4*(pi^2))+5/pi*x(1)-6)^2+10*(1-1/(8*pi))*cos(x(1))+10;
end

% F18
function f = F18(x)
    f=(1+(x(1)+x(2)+1)^2*(19-14*x(1)+3*(x(1)^2)-14*x(2)+6*x(1)*x(2)+3*x(2)^2))*...
    (30+(2*x(1)-3*x(2))^2*(18-32*x(1)+12*(x(1)^2)+48*x(2)-36*x(1)*x(2)+27*(x(2)^2)));
end

% F19
function  f = F19(x)
    aH=[3 10 30;.1 10 35;3 10 30;.1 10 35];cH=[1 1.2 3 3.2];
    pH=[.3689 .117 .2673;.4699 .4387 .747;.1091 .8732 .5547;.03815 .5743 .8828];
    f=0;
    for i=1:4
        f = f - cH(i)*exp(-(sum(aH(i,:).*((x-pH(i,:)).^2))));
    end
end

% F20
function f = F20(x)
    aH=[10 3 17 3.5 1.7 8;.05 10 17 .1 8 14;3 3.5 1.7 10 17 8;17 8 .05 10 .1 14];
    cH=[1 1.2 3 3.2];
    pH=[.1312 .1696 .5569 .0124 .8283 .5886;.2329 .4135 .8307 .3736 .1004 .9991;...
    .2348 .1415 .3522 .2883 .3047 .6650;.4047 .8828 .8732 .5743 .1091 .0381];
    f=0;
    for i=1:4
        f=f-cH(i)*exp(-(sum(aH(i,:).*((x-pH(i,:)).^2))));
    end
end

% F21
function f = F21(x)
    aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
    cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];
    
    f=0;
    for i=1:5
        f=f-((x-aSH(i,:))*(x-aSH(i,:))'+cSH(i))^(-1);
    end
end

% F22
function f = F22(x)
    aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
    cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];
    
    f=0;
    for i=1:7
        f=f-((x-aSH(i,:))*(x-aSH(i,:))'+cSH(i))^(-1);
    end
end

% F23
function f = F23(x)
    aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
    cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];
    
    f=0;
    for i=1:10
        f=f-((x-aSH(i,:))*(x-aSH(i,:))'+cSH(i))^(-1);
    end
end

function f = Ufun(x,a,k,m)
    f = k.*((x-a).^m).*(x>a)+k.*((-x-a).^m).*(x<(-a));
end
