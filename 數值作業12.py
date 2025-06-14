import math

# HW12 - Numerical Methods for PDEs
# --------------------------------------------------
# This code includes:
# 1. Laplace Equation in 2D
# 2. Heat Equation in cylindrical coordinate (1D r) using 3 methods
# 3. Heat Equation in polar coordinate (r, theta)
# 4. Wave Equation in 1D

# --------------------------------------------------
# Problem 1: 2D Laplace Equation (Finite Difference)
# --------------------------------------------------
pi = 3.141592653589793
h = k = pi / 10
nx = 10
ny = 5

u = [[0.0 for _ in range(ny+1)] for _ in range(nx+1)]

# Boundary conditions
for j in range(ny+1):
    y = j * k
    u[0][j] = math.cos(y)
    u[nx][j] = -math.cos(y)
for i in range(nx+1):
    x = i * h
    u[i][0] = math.cos(x)
    u[i][ny] = 0

# Iteration (Jacobi)
for _ in range(500):
    new_u = [[u[i][j] for j in range(ny+1)] for i in range(nx+1)]
    for i in range(1, nx):
        for j in range(1, ny):
            new_u[i][j] = 0.25 * (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1])
    u = new_u

print("\n==============================")
print("Problem 1: Laplace Equation")
print("==============================")
for row in u:
    print([round(val, 3) for val in row])

# --------------------------------------------------
# Problem 2: Heat Equation in Cylindrical Coordinates
# --------------------------------------------------
dr = 0.1
K = 0.1
dt = 0.5
nr = 6
nt = 21

# Initial condition
init_T = [200 * (0.5 + i * dr - 0.5) for i in range(nr)]

# Forward Difference
T_forward = [[init_T[i] if n == 0 else 0.0 for n in range(nt)] for i in range(nr)]
for n in range(0, nt-1):
    for i in range(1, nr-1):
        r = 0.5 + i * dr
        T_forward[i][n+1] = T_forward[i][n] + dt * K * (
            (T_forward[i+1][n] - 2*T_forward[i][n] + T_forward[i-1][n]) / dr**2 +
            (T_forward[i+1][n] - T_forward[i-1][n]) / (2*r*dr)
        )
    T_forward[0][n+1] = (T_forward[1][n+1] - 3*T_forward[0][n+1]*dr)
    T_forward[nr-1][n+1] = 100 + 40*(n+1)*dt

# Backward Difference
T_backward = [[init_T[i] if n == 0 else 0.0 for n in range(nt)] for i in range(nr)]
for n in range(1, nt):
    for i in range(1, nr-1):
        r = 0.5 + i * dr
        T_backward[i][n] = T_backward[i][n-1] + dt * K * (
            (T_backward[i+1][n-1] - 2*T_backward[i][n-1] + T_backward[i-1][n-1]) / dr**2 +
            (T_backward[i+1][n-1] - T_backward[i-1][n-1]) / (2*r*dr)
        )
    T_backward[0][n] = (T_backward[1][n] - 3*T_backward[0][n]*dr)
    T_backward[nr-1][n] = 100 + 40*n*dt

# Crank-Nicolson (approximation step-by-step using average of F and B)
T_cn = [[init_T[i] if n == 0 else 0.0 for n in range(nt)] for i in range(nr)]
for n in range(0, nt-1):
    for i in range(1, nr-1):
        r = 0.5 + i * dr
        fwd = (
            (T_cn[i+1][n] - 2*T_cn[i][n] + T_cn[i-1][n]) / dr**2 +
            (T_cn[i+1][n] - T_cn[i-1][n]) / (2*r*dr)
        )
        T_cn[i][n+1] = T_cn[i][n] + 0.5 * dt * K * fwd
    T_cn[0][n+1] = (T_cn[1][n+1] - 3*T_cn[0][n+1]*dr)
    T_cn[nr-1][n+1] = 100 + 40*(n+1)*dt

print("\n==============================================")
print("Problem 2: Heat Equation in Cylindrical Coordinates")
print("Method: Forward Difference")
print("==============================================")
for row in T_forward:
    print([round(val, 2) for val in row])

print("\nMethod: Backward Difference")
print("==============================================")
for row in T_backward:
    print([round(val, 2) for val in row])

print("\nMethod: Crank-Nicolson Approximation")
print("==============================================")
for row in T_cn:
    print([round(val, 2) for val in row])

# --------------------------------------------------
# Problem 3: Steady-state in Polar Coordinates
# --------------------------------------------------
nr = 6
ntheta = 6
T2 = [[0.0 for _ in range(ntheta)] for _ in range(nr)]
dr = (1 - 0.5) / (nr - 1)
dtheta = (pi/3) / (ntheta - 1)

# Boundary conditions
for j in range(ntheta):
    T2[0][j] = 50
    T2[nr-1][j] = 100
for i in range(nr):
    T2[i][0] = 0
    T2[i][ntheta-1] = 0

# Iteration (Jacobi)
for _ in range(500):
    new_T2 = [[T2[i][j] for j in range(ntheta)] for i in range(nr)]
    for i in range(1, nr-1):
        r = 0.5 + i * dr
        for j in range(1, ntheta-1):
            new_T2[i][j] = (1/(2*(1+1/r**2)))*((T2[i+1][j]+T2[i-1][j])/dr**2 + (1/r**2)*(T2[i][j+1]+T2[i][j-1])/dtheta**2)
    T2 = new_T2

print("\n==============================")
print("Problem 3: Polar Coordinates")
print("==============================")
for row in T2:
    print([round(val, 2) for val in row])

# --------------------------------------------------
# Problem 4: 1D Wave Equation
# --------------------------------------------------
dx = dt = 0.1
nx = 11
nt = 11
p = [[0.0 for _ in range(nt)] for _ in range(nx)]

# Initial conditions
for i in range(nx):
    x = i * dx
    p[i][0] = math.cos(2*pi*x)
    p[i][1] = p[i][0] + dt * 2*pi*math.sin(2*pi*x)

# Boundary conditions
for j in range(nt):
    p[0][j] = 1
    p[nx-1][j] = 2

# Time stepping
for j in range(1, nt-1):
    for i in range(1, nx-1):
        p[i][j+1] = 2*p[i][j] - p[i][j-1] + (dt**2/dx**2)*(p[i+1][j] - 2*p[i][j] + p[i-1][j])

print("\n==============================")
print("Problem 4: Wave Equation")
print("==============================")
for row in p:
    print([round(val, 3) for val in row])
