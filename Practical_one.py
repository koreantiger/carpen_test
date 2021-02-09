import numpy as np
import matplotlib.pyplot as plt

nw = 3.5
no = 2
mu_w = 1e-3
mu_o = 2e-3
krwe = 0.5
kroe = 0.95
S_wc = 0.15
S_or = 0.15
phi = 0.3
u_t = 1e-4
k = 1e-12

# 1
S_init = S_wc
S_inj = 1 - S_or


# 2
ns = 100
sw = np.linspace(S_init,S_inj,ns)
sw_eff = (sw - S_wc)/(1 - S_wc - S_or)
so = 1 - sw
so_eff =(so-S_or)/(1-S_wc-S_or)
#3
fig = plt.figure(1)
Krw = krwe * (sw_eff)**nw
Kro = kroe * (so_eff)**no
lambda_w = k * Krw / mu_w
lambda_o = k * Kro / mu_o
f_w = lambda_w / (lambda_w + lambda_o)
f_o = lambda_o / (lambda_w + lambda_o)
axes1 = fig.add_subplot(1, 3, 1)
axes1.plot(sw, Krw)
axes1.plot(sw, Kro)
axes1.set_title('Relative permeability')
axes1.set_xlabel('S')
axes1.set_ylabel('K')
axes2 = fig.add_subplot(1, 3, 2)
axes2.plot(sw, lambda_w)
axes2.plot(sw, lambda_o)
axes2.set_xlabel('S')
axes2.set_ylabel(' Lambda')
axes2.set_title('Phase Mobility')
axes3 = fig.add_subplot(1, 3, 3)
axes3.plot(sw,f_w)
axes3.plot(sw,f_o)
axes3.set_xlabel('S')
axes3.set_ylabel('f')
axes3.set_title('fractional flow')
plt.show()

# four numerically
dfw = np.zeros(ns)
rr = np.arange(1, ns-1)
dfw[rr] = (f_w[rr+1] - f_w[rr-1])/(sw[rr+1] - sw[rr-1])
dfw[[0, ns-1]] = 2 * dfw[[1, ns -2]] - dfw[[2, ns - 3]]
plt.figure(2)
plt.plot(sw,dfw)
plt.show()


# six jump velocity

vw_jump = (f_w[1:] - f_w[0])/(sw[1:]-S_wc)


# shock saturation
#[shock_velocity, shock_index] = np.max(vw_jump)
shock_velocity = np.max(vw_jump)
shock_index = np.where(vw_jump == np.max(vw_jump))
shock_index = int(shock_index[0])
print(shock_velocity)

# plot fractional flow derivative and the jump velocity
plt.figure(3)
plt.plot(sw, dfw)
plt.plot(sw[1:], vw_jump)
plt.show()

# nine
plt.figure(4)
t = 10
x = u_t / phi * dfw * t
x1 = u_t / phi * dfw * 20

#plt.show()

# ten
sw2 = sw[shock_index:]
dfw2 = dfw[shock_index:]

# Eleven

xx = np.array([S_init, S_init])
xx1 = np.array([np.max(dfw), dfw2[0]])
sw3 = np.concatenate((xx,sw2))
print(sw3)
dfw3 = np.concatenate((xx1, dfw2))
x3 = u_t / phi * dfw3 * t
x3_20 = u_t / phi * dfw3 * 20
# plot(x, sw, x1, sw, x3, sw3, x3_20, sw3)
#plt.figure(5)
plt.plot(x, sw, x1, sw, x3, sw3, x3_20, sw3)
plt.xlabel("Dimensionless position")
plt.ylabel("Water Saturation")
plt.show()

# twelve
"""
#13
Wid = np.zeros(len(dfw2)+2)
Wid[2:] = 1 / dfw2

print(len(Wid))
L=50
t2 = L*phi/u_t/dfw2

plt.figure(6)
plt.plot(t2, sw2[1:])
#plt.show()


# 14
fw2 = np.zeros(len(f_w[shock_index:])+1)
fw2[0] = int(f_w[shock_index])
fw2[1:] = f_w[shock_index:]
print(len(fw2))
ss = np.zeros(len(sw2)+ 1)
ss[0] = S_init
ss[1:] = sw2
print(len(ss))
Npd = (1 - fw2) * Wid + (ss - S_wc)
plt.figure(7)
plt.plot(Wid, Npd)
plt.xlabel('Water PV injected')
plt.ylabel('Oil PV injected')
plt.show()

# 15

Swmed= np.zeros(len(dfw2)+1)
Swmed[0] = S_init
x = (1 - f_w[shock_index:-1]) / dfw2
print(len(sw2))
print(len(dfw2))
print(len(f_w[shock_index:-1]))
Swmed[1:] = (1 - f_w[shock_index:-1]) / dfw2 + sw2[1:]
plt.figure(7)
plt.plot(Wid[1:], Swmed)
plt.xlabel('Water PV injected')
plt.ylabel('Average Sw')
plt.show()

"""

#13
Wid = np.zeros(len(dfw2)+1)
Wid[1:] = 1 / dfw2
L=50
t2 = L*phi/u_t/dfw2
plt.figure(6)
plt.plot(t2[:-1], sw2[:-1])
plt.show()

fw2 = np.zeros(len(f_w[shock_index:])+1)
fw2[0] = int(f_w[shock_index])
fw2[1:] = f_w[shock_index:]
print(len(fw2))
ss = np.zeros(len(sw2)+ 1)
ss[0] = S_init
ss[1:] = sw2
print(len(ss))
Npd = (1 - fw2) * Wid + (ss - S_wc)
plt.figure(7)
plt.plot(Wid[:-1], Npd[:-1])
plt.xlabel('Water PV injected')
plt.ylabel('Oil PV injected')
plt.show()

# 15

Swmed= np.zeros(len(dfw2)+1)
Swmed[0] = S_init
x = (1 - f_w[shock_index:-1]) / dfw2
print(len(sw2))
print(len(dfw2))
print(len(f_w[shock_index:-1]))
Swmed[1:] = (1 - f_w[shock_index:-1]) / dfw2 + sw2[1:]
plt.figure(7)
plt.plot(Wid[1:], Swmed)
plt.xlabel('Water PV injected')
plt.ylabel('Average Sw')
plt.show()