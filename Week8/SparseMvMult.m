
function y = SparseMvMult(nzA,ir,ic,x)
N = size(ir,1)-1 ;

r = 1;
ind_nextrow = ir(r+1);
y =zeros(N,1);


y_i = 1; 

for j=1:size(nzA,1) 
    if j ==ind_nextrow 
        
        y_i = y_i +1;
        r = r+1 ; 
        ind_nextrow = ir(r+1);
    end
    y(y_i) = y(y_i) + x(ic(j))*nzA(j,1) ;
    
end
