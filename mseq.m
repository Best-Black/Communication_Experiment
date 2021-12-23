function [seq]=mseq(coef)
m=length(coef);
len=2^m-1;  %得到最终生成的m序列的长度    
backQ=0; 
seq=zeros(1,len); 
registers = [1 zeros(1, m-2) 1]; 
for i=1:len
    seq(i)=registers(m);
    backQ = mod(sum(coef.*registers) , 2); 
    registers(2:length(registers)) = registers(1:length(registers)-1); 
    registers(1)=backQ; 
end
end