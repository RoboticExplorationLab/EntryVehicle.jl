#Function returns the quaternion associated with a rotation matrix
using LinearAlgebra

#test part
θ = 0.0*pi/180

M = [-sin(θ) cos(θ) 0.0;
     0.0 0.0 1.0;
     cos(θ) sin(θ) 0.0]

Q = mat2quat(M)

#function
function mat2quat(R)
    r11, r12, r13 = R[1,1], R[1,2], R[1,3]
    r21, r22, r23 = R[2,1], R[2,2], R[2,3]
    r31, r32, r33 = R[3,1], R[3,2], R[3,3]
    N_0 = sqrt(1+r11+r22+r33)
    Q_0 = 0.5*[N_0; (r23-r32)/N_0; (r31-r13)/N_0; (r12-r21)/N_0]
    N_1 = sqrt(1+r11-r22-r33)
    Q_1 = 0.5*[(r23-r32)/N_1; N_1;(r12+r21)/N_1;(r31+r13)/N_1]
    N_2 = sqrt(1-r11+r22-r33)
    Q_2 = 0.5*[(r31-r13)/N_2;(r12+r21)/N_2; N_2;(r23+r32)/N_2]
    N_3 = sqrt(1-r11-r22+r33)
    Q_3 = 0.5*[(r12-r21)/N_3; (r31+r13)/N_3; (r23+r32)/N_3; N_3]
    if r22>=-r33 && r11>-r22 && r11>-r33
        Q = Q_0
    elseif r22<=-r33 && r11>r22 && r11>r33
        Q = Q_1
    elseif r22>=r33 && r11<r22 && r11<-r33
        Q = Q_2
    elseif r22<=r33 && r11<-r22 && r11<r33
        Q = Q_3
    end
    return Q
end
