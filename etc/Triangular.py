import random
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize, integrate
from scipy.linalg import eig, eigvalsh_tridiagonal
from math import *
from sklearn.metrics.pairwise import cosine_similarity


class Triangular:
    def __init__(self, Js, M):
        self.Js = Js
        self.M = M
        self.dim = (1 << M)
        self.dim2 = self.dim*2

    def product_D(self, temperature, v):
        w1 = np.array([1.0, np.exp(-2.0*self.Js[0]/temperature)])
        for s in range(self.dim):
            for i in range(self.M):
                j = (i+1) % self.M
                sigma_i = (s >> i) & 1
                sigma_j = (s >> j) & 1
                v[s] = v[s] * w1[sigma_i ^ sigma_j]

    def product_U(self, temperature, v, vtmp):  # v, vtmpはサイズ2倍にしておく
        w2 = np.array([1.0, np.exp(-2.0*self.Js[1]/temperature)])
        w3 = np.array([1.0, np.exp(-2.0*self.Js[2]/temperature)])

        # U1: dim expansion
        vtmp[:] = 0.0
        for s in range(self.dim):  # before
            sm, sp = s << 1, (s << 1) | 1
            sigma_d, sigma_u = s & 1, (s >> (self.M-1)) & 1
            vtmp[sp] += v[s]*w2[1 ^ sigma_d]*w3[1 ^ sigma_u]
            vtmp[sm] += v[s]*w2[0 ^ sigma_d]*w3[0 ^ sigma_u]
        for s in range(self.dim2):
            v[s] = vtmp[s]

        # U2-UM
        for i in range(1, self.M):
            vtmp[:] = 0.0
            for s in range(self.dim2):  # before
                sm, sp = s & ~(1 << i), s | (1 << i)
                sigma_d, sigma_u = (s >> i) & 1, (s >> (i+1)) & 1
                vtmp[sp] += v[s]*w2[1 ^ sigma_u]*w3[1 ^ sigma_d]
                vtmp[sm] += v[s]*w2[0 ^ sigma_u]*w3[0 ^ sigma_d]
            for s in range(self.dim2):
                v[s] = vtmp[s]

        # U(M+1): dim shrink
        vtmp[:] = 0.0
        for s in range(self.dim):  # after
            sm, sp = s, s | self.dim
            vtmp[s] += v[sm]
            vtmp[s] += v[sp]
        for s in range(self.dim):
            v[s] = vtmp[s]

    def product_U_T(self, temperature, v, vtmp):
        w2 = np.array([1.0, np.exp(-2.0*self.Js[1]/temperature)])
        w3 = np.array([1.0, np.exp(-2.0*self.Js[2]/temperature)])

        vtmp[:] = 0.0
        for s in range(self.dim):
            sm, sp = s, s | self.dim
            vtmp[sm] += v[s]
            vtmp[sp] += v[s]
        for s in range(self.dim2):
            v[s] = vtmp[s]

        for i in reversed(range(1, self.M)):
            vtmp[:] = 0.0
            for s in range(self.dim2):
                sm, sp = s & ~(1 << i), s | (1 << i)
                sigma_d, sigma_u = (s >> i) & 1, (s >> (i+1)) & 1
                vtmp[s] += v[sp]*w2[1 ^ sigma_u]*w3[1 ^ sigma_d]
                vtmp[s] += v[sm]*w2[0 ^ sigma_u]*w3[0 ^ sigma_d]
            for s in range(self.dim2):
                v[s] = vtmp[s]

        vtmp[:] = 0.0
        for s in range(self.dim):
            sm, sp = s << 1, (s << 1) | 1
            sigma_d, sigma_u = s & 1, (s >> (self.M-1)) & 1
            vtmp[s] += v[sp]*w2[1 ^ sigma_d]*w3[1 ^ sigma_u]
            vtmp[s] += v[sm]*w2[0 ^ sigma_d]*w3[0 ^ sigma_u]
        for s in range(self.dim):
            v[s] = vtmp[s]

    def product(self, temperature, v, vtmp, transpose=False):
        if transpose:
            self.product_D(temperature, v)
            self.product_U_T(temperature, v, vtmp)
        else:
            self.product_U(temperature, v, vtmp)
            self.product_D(temperature, v)

    def power_iter_1(self, temperature, v0, eps=1e-8):
        # lambda_1 & v1R
        vo = v0
        vn = np.zeros_like(vo)
        vtmp1, vtmp2 = np.zeros((2, self.dim2))
        lmd_o, lmd_n = 0.0, 1.0
        while True:
            vtmp1[:] = 0.0
            vtmp1[:self.dim] = vo
            self.product(temperature, vtmp1, vtmp2)
            vn[:] = vtmp1[:self.dim]
            lmd_n = np.dot(vn, vn)/np.dot(vn, vo)
            vo[:] = vn/np.linalg.norm(vn)
            if np.abs(lmd_n-lmd_o) < eps:
                break
            lmd_o = lmd_n
        v1R, lmd1 = vo.copy(), lmd_n

        # v1L
        vo = v0
        vn = np.zeros_like(vo)
        vtmp1, vtmp2 = np.zeros((2, self.dim2))
        lmd_o, lmd_n = 0.0, 1.0
        while True:
            vtmp1[:] = 0.0
            vtmp1[:self.dim] = vo
            self.product(temperature, vtmp1, vtmp2, transpose=True)
            vn[:] = vtmp1[:self.dim]
            lmd_n = np.dot(vn, vn)/np.dot(vn, vo)
            vo[:] = vn/np.linalg.norm(vn)
            if np.abs(lmd_n-lmd_o) < eps:
                break
            lmd_o = lmd_n
        v1L = vo

        return lmd1, v1R, v1L

    def decrease_term(self, lmd, vR, vL, v, vtmp):
        vL_v = np.dot(vL, v)
        vtmp[:self.dim] = lmd * vL_v * vR

    def power_iter_2(self, temperature, v0, lmd1, v1R, v1L, eps=1e-8):
        vo = v0
        vn = np.zeros_like(vo)
        vtmp1, vtmp2 = np.zeros((2, self.dim2))
        lmd_o, lmd_n = 0.0, 1.0
        cnt = 0
        while True:
            vtmp1[:] = 0.0
            vtmp1[:self.dim] = vo
            self.product(temperature, vtmp1, vtmp2)
            vn[:] = vtmp1[:self.dim]
            self.decrease_term(lmd1, v1R, v1L, vo, vtmp1)
            vn[:] = vn - vtmp1[:self.dim]
            lmd_n = np.dot(vn, vn)/np.dot(vn, vo)
            vo[:] = vn/np.linalg.norm(vn)
            if np.abs(lmd_n-lmd_o) < eps:
                break
            lmd_o = lmd_n
            cnt += 1
        return lmd_n, vo

    def calc_xi(self, temperature):
        v0 = np.random.rand(self.dim)
        lmd1, v1R, v1L = self.power_iter_1(temperature, v0)
        v0 = np.random.rand(self.dim)
        lmd2, v2 = self.power_iter_2(temperature, v0, lmd1, v1R, v1L)
        col = lmd2/lmd1
        return 1.0/np.abs(np.log(col))

# for confirmation
    def fill_T(self, temperature, T, transpose=False):
        v, vtmp = np.zeros((2, self.dim2))
        for i in range(self.dim):
            v[:] = 0.0
            v[i] = 1.0
            self.product(temperature, v, vtmp, transpose)
            T[:, i] = v[:self.dim]

    def calc_T_bruteforce(self, temperature):
        w = self.Js.copy()
        w = w/temperature
        w = np.vstack([w, -w]).T
        w = np.exp(w)

        T = np.ones((self.dim, self.dim))
        for s1 in range(self.dim):
            for s2 in range(self.dim):
                for i in range(self.M):
                    j = (i-1) % self.M
                    sigma = (s1 >> i) & 1
                    sigma_s = np.array([s1 >> j, s2 >> i, s2 >> j]) & 1
                    antiferro = sigma ^ sigma_s
                    for t in range(3):
                        T[s1, s2] *= w[t, antiferro[t]]
        return T

    def diag_bruteforce(self, temperature):
        T = np.zeros((self.dim, self.dim))
        self.fill_T(temperature, T)
        lmd, v = eig(T)
        lmd = np.array(lmd)
        lmd = np.real(lmd)
        v = np.array(v).T
        v = np.real(v)
        arg = np.argsort(lmd)[::-1]
        return lmd[arg], v[arg]


M = 8
Js_raw = np.array([1.0, 1.0, 1.0])
Js = Js_raw/2.0/np.sum(Js_raw)
temperature = 1.0
t = Triangular(Js, M)

v0 = np.random.rand(t.dim)
lmd1, v1R, v1L = t.power_iter_1(temperature, v0)
v0 = np.random.rand(t.dim)
lmd2, v2 = t.power_iter_2(temperature, v0, lmd1, v1R, v1L)

lmd_s, v_s = t.diag_bruteforce(temperature)


np.set_printoptions(precision=4)


print("J = (%d, %d, %d), M =" % tuple(np.rint(Js_raw)), M)

print("\nλの計算値: %.8f, %.8f" % (lmd1, lmd2))
print("λ1の相対誤差: %.2e %%" % ((lmd1-lmd_s[0])/lmd_s[0]*100))
print("λ2の相対誤差: %.2e %%" % ((lmd2-lmd_s[1])/lmd_s[1]*100))
print("\nv1のCos類似度: ", np.dot(v1R, v_s[0]) /
      np.linalg.norm(v1R)/np.linalg.norm(v_s[0]))
print("v2のCos類似度: ", np.dot(v2, v_s[1]) /
      np.linalg.norm(v2)/np.linalg.norm(v_s[1]))
print("\nv1（左右）のCos類似度: ", np.dot(v1R, v1L) /
      np.linalg.norm(v1R)/np.linalg.norm(v1L))
