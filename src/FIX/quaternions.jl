using LinearAlgebra

#File contains useful functions associated with quaternions

function qmult(q1,q2)
    #scalar first
    s1 = q1[1]
    s2 = q2[1]
    v1 = q1[2:4]
    v2 = q2[2:4]
    q3 = [s1*s2 - v1'*v2;
          s1*v2 + s2*v1 + cross(v1, v2)]
end

function qrot(q,x)
    xrot = x + 2.0*cross(q[2:4], cross(q[2:4],x) + q[1]*x)
end

function qconj(q)
    qc = [q[1]; -q[2:4]]
end

#FUNCTION BELOW SHOULD WORK AND BE MORE ROBUST ==> TO BE FIXED

function mat2quat_2(R)
    Q = zeros(4)
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



function mat2quat(M)
    #takes a rotation matrix as an input and returns a quaternion
    q1 = 0.5 * sqrt(1 + M[1,1] + M[2,2] + M[3,3])
    q2 = (M[2,3] - M[3, 2]) / (4 * q1)
    q3 = (M[3,1] - M[1,3]) / (4 * q1)
    q4 = (M[1,2] - M[2,1]) / (4 * q1)

    q = [q1, q2, q3, q4]
    q = q/norm(q)
    return q
end



function log_quat(q)
    if norm(q[2:4]) != 0.0
        Q = [0.0;acos(q[1]/norm(q))*(q[2:4]/(norm(q[2:4])))]
    else
        Q = [0.0; 0.0; 0.0; 0.0]
    end
    return Q
end

function exp_quat(q) #v is a rotation vector
    if norm(q[2:4]) != 0.0
        Q = exp(q[1])*[cos(norm(q[2:4])); (q[2:4]/norm(q[2:4]))*sin(norm(q[2:4]))]
    else
        Q = [1.0; 0.0; 0.0; 0.0]
    end
    return Q
end

function quat2euler_error(q, qref)
    #qref is the reference that we use
    #return vector euler error between q and qref
    V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    δq = qmult(qconj(qref), q)
    e = 2*V*log_quat(δq)
    return e
end

function euler2quat_error(e, qref)
    #e is the euler vector of the rotation diff between q and qref
    #this function returns q given the error and the reference
    V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    q = qmult(qref, exp_quat(V'*e/2))
    return q
end

function quat2euler(q)
    #return the euler vector associated with quaternion q
    V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    e = 2*V*log_quat(q)
    return e
end

function euler2quat(e)
    V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    q = exp_quat(V'*e/2)
    return q
end
