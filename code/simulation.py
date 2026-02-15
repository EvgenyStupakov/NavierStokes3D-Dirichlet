import numpy as np

# Parameters
L, nx, ny, nz = 1.0, 5, 5, 5
nu, n_modes, dt, n_steps = 0.1, 2, 0.05, 6
x = np.linspace(0,L,nx); y = np.linspace(0,L,ny); z = np.linspace(0,L,nz)
X,Y,Z = np.meshgrid(x,y,z,indexing='ij')

def F_nml(t,n,m,l):
    return 5/(n*m*l) if t<0.1 else 0

for t_step in range(n_steps):
    t = t_step*dt
    u = np.zeros_like(X)
    v = np.zeros_like(X)
    w = np.zeros_like(X)
    for n in range(1,n_modes+1):
        for m in range(1,n_modes+1):
            for l in range(1,n_modes+1):
                lam = np.pi**2*(n**2+m**2+l**2)
                a = F_nml(t,n,m,l)/lam*(1-np.exp(-nu*lam*min(t,0.1))) \
                    * np.exp(-nu*lam*max(0,t-0.1))
                u += a*np.sin(n*np.pi*X)*np.sin(m*np.pi*Y)*np.sin(l*np.pi*Z)
                v += a*np.sin(m*np.pi*X)*np.sin(l*np.pi*Y)*np.sin(n*np.pi*Z)
                w += a*np.sin(l*np.pi*X)*np.sin(n*np.pi*Y)*np.sin(m*np.pi*Z)
    speed = np.sqrt(u**2+v**2+w**2)
    print(f"\n--- t={t:.2f} ---\n", speed.round(2))

