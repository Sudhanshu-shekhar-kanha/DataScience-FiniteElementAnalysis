function [ K ] = K_ele(B,C,t,AREA)

K=transpose(B)*C*B*t*AREA;

end
