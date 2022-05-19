% 1 part 

[k,t,x]=readfile('robot.key') 
sis=size(x.data)
i=1
z=1
m=sis(1)
while i<m
    for b=1:sis(2)
q(z,b)=x.data(i,b)
    end
    i=i+2
    z=z+1
end
z=1    
i=2
while i<m+1
    for b=1:sis(2)
q_vel(z,b)=x.data(i,b) % velocity
    end
    i=i+2
    z=z+1
end

q_length=size(q)
step_size=(k-1)/(t-1) %step u

for j=1:q_length(2) % for each joint

p=find_polynom(q(:,j),q_vel(:,j))
    
Q(:,j)=find_new_point(p,step_size,t)
end
%write in file
val=writefile(Q,'robot/robot.ang', t)

% 2 part
[k,t,obj]=readfile('object.key') 
% create a two matrix xyz and trans 
sis=size(obj.data)
step_size=(k-1)/(t-1)
m=sis(1)
n=sis(2)
for i=0:k-1
    X_obj(i+1)=obj.data((i*3+1),n)
    Y_obj(i+1)=obj.data((i*3+2),n)
    Z_obj(i+1)=obj.data((i*3+3),n)
end

X_vel=find_velocity(X_obj)
Y_vel=find_velocity(Y_obj)
Z_vel=find_velocity(Z_obj)

px=find_polynom(X_obj,X_vel)
py=find_polynom(Y_obj,Y_vel)
pz=find_polynom(Z_obj,Z_vel)

X=find_new_point(px,step_size, t)
Y=find_new_point(py,step_size, t)
Z=find_new_point(pz,step_size, t)

for index=1:k
for i=1:3
    for j=1:3
        rotationMatrix(i,j)=obj.data((i+(index-1)*3),j)
    end
end
quat = quaternion(rotationMatrix,'rotmat','frame')
vector=rotvec(quat)
for index2=1:3 
vector_matrix(index,index2)=vector(index2)
end
end

vector_vel(:,1)=find_velocity(vector_matrix(:,1))
vector_vel(:,2)=find_velocity(vector_matrix(:,2))
vector_vel(:,3)=find_velocity(vector_matrix(:,3))

p1=find_polynom(vector_matrix(:,1),vector_vel(:,1))
p2=find_polynom(vector_matrix(:,2),vector_vel(:,2))
p3=find_polynom(vector_matrix(:,3),vector_vel(:,3))

I=find_new_point(p1,step_size, t)
J=find_new_point(p2,step_size, t)
K=find_new_point(p3,step_size, t)


for loop=1:t 
  
quat1 = quaternion([I(loop) J(loop) K(loop)],'rotvec')
matrix=rotmat(quat1,'frame')

for i=1:3
    for j=1:3
        OBJ_MATRIX((loop-1)*3+i,j)=matrix(i,j)
    end
end

end

for loop=1:t 
OBJ_MATRIX((loop-1)*3+1,4)=X(loop)
OBJ_MATRIX((loop-1)*3+2,4)=Y(loop)
OBJ_MATRIX((loop-1)*3+3,4)=Z(loop)

end

val=writefile(OBJ_MATRIX,'object/object.traj', t)

function [k,t,x]=readfile(s) %read file input file output valeu of key frame and need one
delimiterIn = ' ';
headerlinesIn = 1;
x1=importdata(s )
k=x1(1,1)  % key frame
t=x1(1,2)   % how much frame I need
x=importdata(s, delimiterIn, 1 )
end

function val=writefile(Q,fName, t) % rewrite file in fName directory , write a label t and matrix Q
str = num2str(t);      
mat1 = Q;             
        %# a file name
fid = fopen(fName,'w');            
if fid ~= -1
  fprintf(fid,'%s\r\n',str);       %# print the string
  fclose(fid);                    
end
dlmwrite(fName,mat1,'-append','delimiter',' ');  %print the matrix'
val=1
end

function VEL=find_velocity(point) % find the velocity in D(i)=(P(i+1)-P(i-1))/2 in first and last point velocity iz zero
VEL(1)=0
for i=2:length(point)-1
VEL(i)=(point(i+1)-point(i-1))/2
end
VEL(length(point))=VEL(length(point)-1)+0,1
end
function massive=find_new_point(p,step_size,t ) % find a new point between key points use 2 point and polynome and step size
l=1
u=0
for loop=1:t % find t new point


    if u<1
        massive(loop)=u*u*u*p(l, 1)+u*u*p(l,2)+u*p(l,3)+p(l,4)
       
    else
        u=u-1
        l=l+1 % start with new polynom
         massive(loop)=u*u*u*p(l, 1)+u*u*p(l,2)+u*p(l,3)+p(l,4)
        
    end
   
    u=step_size+u
end
end

function p=find_polynom(point, point_vel) % find hermet polynom

h1=[2 -3 0  1]
h2=[-2 3 0  0]
h3=[1 -2 1  0]
h4=[1 -1 0  0]

for i=1:length(point)-1  %find polynome for one joint
p_Trah=h1*point(i)+h2*point(i+1)+h3*point_vel(i)+h4*point_vel(i+1)
    p(i,:)=p_Trah 
end
end