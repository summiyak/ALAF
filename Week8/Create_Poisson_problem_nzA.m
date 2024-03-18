function [ nzA, ir, ic ] = Create_Poisson_problem_nzA( N )
n=N;
N = N^2;

size_nzA = (N + 2*2*(N-n));
nzA = zeros(size_nzA,1);
nzA_i = 1 ;  %index for iterating through the array nzA

ir = zeros(N,1);
ir(1)= 1;

ic = zeros(size_nzA,1);
ic_i =1 ; %index for iterating through the array ic

for i = 1: N
        if i-n>0  
            nzA(nzA_i) = -1;
            nzA_i = nzA_i+1;
            ic(ic_i) = i-n;
            ic_i = ic_i +1;
        end
        if (i-1 >0 && mod(i-1,n)>0)
            nzA(nzA_i) = -1;
            nzA_i = nzA_i+1;
            ic(ic_i) = i-1;
            ic_i = ic_i +1;
        end

        nzA(nzA_i) = 4;
        nzA_i = nzA_i+1;
        ic(ic_i) = i;
        ic_i = ic_i +1;

        if (i+1 <=N && mod(i,n)>0)  
            nzA(nzA_i) = -1;
            nzA_i = nzA_i+1;
            ic(ic_i) = i+1;
            ic_i = ic_i +1;
        end

        if (i+n) <=N
            nzA(nzA_i) = -1;
            nzA_i = nzA_i+1;
            ic(ic_i) = i+n;
            ic_i = ic_i +1;
        end
        ir(i+1)= nzA_i;
end