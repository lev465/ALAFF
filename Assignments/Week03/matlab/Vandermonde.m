A = [1 0; 0 1; 1 1]
b = [1; 1; 0]

A_H = transpose(A)
A_HA = A_H * A
A_HA_inv = inv(A_HA)

x_h = A_HA_inv*A_H*b

b_h = A*x_h

Q = [1/sqrt(2) -1/sqrt(6); 0 sqrt(2/3); 1/sqrt(2) 1/sqrt(6)]

Q_H = transpose(Q)

R = Q_H*A

R_H = transpose(R)

y = Q_H*b

x_h_QR = [1/3; 1/3]

R*x_h_QR

q_0 = [sqrt(2)/2; sqrt(2)/2]
q_1 = [-sqrt(2)/2; sqrt(2)/2]

(6/sqrt(2)) * q_0 - (2/sqrt(2)) * q_1

U = [1 2 3 4; 0 3 -3 4; 0 0 3 4; 0 0 0 4]
U_inv = inv(U)
U*U_inv
U_inv*U