function P=Ion_Pressure(C_H,C_S,alpha_B,c_p)
    P=sqrt((alpha_B*c_p)^2+4*(C_S+C_H)^2)-2*(C_S+C_H);
end

