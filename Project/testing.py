import numpy as np

q = 2. ; m = 1.
r_n = np.array([0.,0.,0.],dtype=float) ; v_n = np.array([3000.,0.,0.],dtype=float)
B = np.array([0,0,1.5],dtype=float) ; a_n = q/m * np.cross(v_n,B)
dt = 0.1

k_1v = a_n*dt
k_1x = v_n*dt

k_2v = np.cross(v_n+0.5*k_1v,B) * q/m * dt
k_2x = (v_n + 0.5*k_1v) * dt

k_3v = np.cross(v_n+0.5*k_2v,B) * q/m * dt
k_3x = (v_n + 0.5*k_2v) * dt

k_4v = np.cross(v_n+k_3v,B) * q/m * dt
k_4x = (v_n + k_3v) * dt

velocity = v_n + 1/6 * (k_1v + 2*k_2v + 2*k_3v + k_4v)
position = r_n + 1/6 * (k_1x + 2*k_2x + 2*k_3x + k_4x)

print(velocity)
print(position)