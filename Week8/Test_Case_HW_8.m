
N =3;
x = randi(9,N^2,1) ;
%x= ones(N^2,1); 

[ nzA, ir, ic ] = Create_Poisson_problem_nzA(N);
y = SparseMvMult(nzA,ir,ic,x);

%% Generate Matrix A and extract nzA, ir and ic from the whole Matrix
n=N;
N = N^2;

A= zeros(N,N) ;
for i = 1: N
        if i-n>0  
            A(i,i-n) = -1;
        end
        
        if (i-1 >0 && mod(i-1,n)>0)
            A(i,i-1)= -1;
        end

        A(i,i)= 4;

        if (i+1 <=N && mod(i,n)>0)
            A(i,(i+1)) = -1;
        end
        if (i+n) <=N 
            A(i,(i+n)) = -1;
        end
end

% % Since A is an SPD and A' = A hence row indices provided are the same as 
% % the column indices
S= sparse(A);
[ic_t,ir_temp,nzA_t] = find(S);   
ir_size = size(ir_temp) ; 
t = ir_temp(1,1);
ir_t = zeros(N+1,1);
ir_t(1,1)  = t;
ir_t(N+1,1) = size(nzA_t,1)+1; 
k =2 ;
for z= 2:ir_size 
    if ir_temp(z,1) ~= t 
        ir_t(k,1) = z;
        t = ir_temp(z,1);
        k = k+1;
    end

end

y_t =A*x;
%% Print Matrix A, nzA, ir, ic and y
A 
nzA 
ir 
ic 
x
y
%% Test Case 1 - Check the values of nzA
test_Case_1 = all(nzA == nzA_t);

if test_Case_1==1 
    disp('Test Case 1 - Validate the values of nzA - Passed!') 
else 
    disp('Test Case 1 - Validate the values of nzA - Failed!') 
end

%% Test Case 2 - Check the values of ic
test_Case_2 = all(ic == ic_t);
if test_Case_2==1 
    disp('Test Case 2 - Validate the values of ic - Passed!') 
else 
    disp('Test Case 2 - Validate the values of ic - Failed!') 
end
%% Test Case 3 - Check the values of ir
test_Case_3 = all(ir == ir_t);
if test_Case_3==1 
    disp('Test Case 3 - Validate the values of ir - Passed!') 
else 
    disp('Test Case 3 - Validate the values of ir - Failed!') 
end

%% Test Case 4 - Check the values of y
test_Case_4 = all(y== y_t);
if test_Case_3==1 
    disp('Test Case 4 - Validate the values of y - Passed!') 
else 
    disp('Test Case 4 - Validate the values of y - Failed!') 
end