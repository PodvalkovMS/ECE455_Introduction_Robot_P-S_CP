f=dlmread('arm') 
n=f(1,1)
lamda=f(1,2)
LMAX=0
for i=2:length(f(:,1))
l(i-1)=f(i,1)
teta(i-1)=f(i,2)
LMAX=f(i,1)+LMAX
end

x1=dlmread('trajectory')
m=x1(1,1)
lamda_trajectory=x1(1,2)
for i=2:length(x1(:,1))
R=sqrt(x1(i,1)^2+x1(i,2)^2)
if R>LMAX
x(i-1)=x1(i,1)*LMAX/R
y(i-1)=x1(i,2)*LMAX/R
else if R<l(1)
x(i-1)=x1(i,1)*l(1)/R
y(i-1)=x1(i,2)*l(1)/R            
    else
    x(i-1)=x1(i,1)
    y(i-1)=x1(i,2)
    end    
end
end

deltaMax=l(1)/100
index=1
teta3=teta
Length_XDelta=1
while index<=m
while Length_XDelta>=deltaMax/100
xa=forward_kenimatic(teta3,l,n)
XDelta(1)=x(index)-xa(1)
XDelta(2)=y(index)-xa(2)
Length_XDelta=sqrt(XDelta(1)^2+XDelta(2)^2)
if Length_XDelta>deltaMax
XDelta=XDelta*deltaMax/Length_XDelta        
end 
Jocobi=find_revers_jacobi(find_real_jacobian(teta3,l,n), lamda)
Teta_Delat=Jocobi*XDelta'
teta3=teta3+Teta_Delat'
end
result(index,:)=teta3
index=index+1
Length_XDelta=1
end

val=writefile(result,'angles')

function val=writefile(Q,fName) % rewrite file in fName directory , write a label t and matrix Q
    
mat1 = Q;             
        %# a file name
fid = fopen(fName,'w');            
dlmwrite(fName,mat1,'-append','delimiter',' ');  %print the matrix'
val=1
end



function TETE=find_teta(teta,n)
TETE(1)=teta(1)
for i=2:n
  TETE(i)=teta(i)+TETE(i-1)    
end
end

function L_and_trig=find_l_and_trig(l,teta,n)
for i=1:n
  L_and_trig(1,i)=-l(i)*sin(teta(i))
  L_and_trig(2,i)=l(i)*cos(teta(i))
end
end

function J=find_jacobian(trig,n)
J(1,n)=trig(1,n)
J(2,n)=trig(2,n)
for i=1:n-1
J(1,n-i)=J(1,n-i+1)+trig(1,n-i)
J(2,n-i)=J(2,n-i+1)+trig(2,n-i)
end
end

function J=find_real_jacobian(teta,l,n)
TETE=find_teta(teta,n)
trig=find_l_and_trig(l,TETE,n)
J=find_jacobian(trig,n)
end


function coordinat=forward_kenimatic(teta,l,n)
TETE=find_teta(teta,n)
trig=find_l_and_trig(l,TETE,n)
trig2(1,:)=trig(2,:)
trig2(2,:)=-trig(1,:)
coordinat(1)=sum(trig2(1,:))
coordinat(2)=sum(trig2(2,:))
end

function Jlimda=find_revers_jacobi(J,lambda)
Jlimda=J'*inv(J*J'+lambda^2*eye(2,2))
end