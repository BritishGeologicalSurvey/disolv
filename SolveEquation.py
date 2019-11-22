# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 09:12:32 2019

@author: sarcol
"""
import sys
import numpy as np
from scipy.integrate import odeint
from scipy.sparse import diags
from scipy.optimize import minimize


def WriteMatrix(N_nodes, inflows, outflows, z, alpha, Cc, Nflows, Dd):

    """
    Creates matrices to solve
    """

    if abs(np.sum(inflows[:, 1]) - np.sum(outflows[:, 1]))/np.sum(
            inflows[:, 1]) > 0.001:
        sys.exit("Inflow does not equal outflow")

# ---------Put inflows and outflows in arrays------------

    qout = np.zeros(N_nodes)
    qin = np.zeros(N_nodes)

    for i in range(Nflows):

        nodein = int(inflows[i, 0]/z)
        nodeout = int(outflows[i, 0]/z)

        qin[nodein] = inflows[i, 1]
        qout[nodeout] = outflows[i, 1]

# --------------Calculate flows at every node------------

    Qculm = 0
    Flow_front = np.zeros(N_nodes)
    Flow_back = np.zeros(N_nodes)

    for i in range(N_nodes - 1):
        Flow = Qculm + (qin[i] - qout[i])
        Flow_front[i] = Flow
        Flow_back[i + 1] = Flow
        Qculm = Flow

    qin = qin/z
    qout = qout/z

    Dfront = abs(Flow_front * alpha) + Dd
    Dback = abs(Flow_back * alpha) + Dd

# -------------------Write function---------------------

    p1 = Dfront + Dback
    p2 = Dfront
    p3 = Dback
    p4 = Flow_front
    p5 = Flow_back
    p6 = qin * Cc
    p7 = qout

    Eqs = np.zeros([N_nodes, 3])

    for i in range(N_nodes - 2):

        j = i + 1

# ------In order ci , ci+1, ci-1
        Eqs[j, :] = np.array([-p1[j]/z**2 - p4[j]/(2*z) + p5[j]/(2*z) - p7[j],
                        p2[j]/z**2 - p4[j]/(2*z), p3[j]/z**2 + p5[j]/(2*z)])

# --First and last lines
    Eqs[0, :] = [-p1[0]/z**2 - p4[0]/(2*z) + p5[0]/(2*z) - p7[0],
                    (p2[0] + p3[0])/z**2 - (p4[0] - p5[0])/(2*z), 0]
    Eqs[-1, :] = [-p1[-1]/z**2 - p4[-1]/(2*z) + p5[-1]/(2*z) - p7[-1], 0,
                    (p2[-1]+p3[-1])/z**2 - (p4[-1] - p5[-1])/(2*z)]

# --Eqs are in order ci , ci+1, ci-1
    Matrix = diags([Eqs[1:, 2], Eqs[:, 0], Eqs[:-1, 1]], [-1, 0, 1], shape=(
            N_nodes, N_nodes)).todense()

    return Matrix, p6


def forward(N_nodes, inflows, outflows, z, alpha, Cc, Nflows, in_con, t, A,
            Dd):
    """
    Forward model
    """

    Matrix, p6 = WriteMatrix(N_nodes, inflows, outflows, z, alpha, Cc, Nflows,
                             Dd)[0], WriteMatrix(N_nodes, inflows, outflows, z,
                             alpha, Cc, Nflows, Dd)[1]

# --Make function
    def dX_dt(sm, t):
        return 1/A * np.squeeze(np.asarray(np.dot(Matrix, sm) + p6))

# --Solver
    wsol = odeint(dX_dt, in_con, t, atol=1e-08, rtol=1e-06, mxstep=5000000)

    return wsol


def inverse(N_nodes, inflows, outflows, z, alpha, Cc, Nflows, in_con, t, A,
            Obs, Bound, FracBounds, method, Dd, minimise_param):
    """
    Inverse model
    """
    def CalculateRMSE(invec):

        alpha = invec[0]
        variable_inflow = invec[1:int(1 + Nflows)]
        variable_outflow = invec[int(1 + Nflows):-int(1 + Nflows)]

        inflows[:, 1] = variable_inflow/np.sum(variable_inflow) * invec[-1]
        outflows[:, 1] = variable_outflow/np.sum(variable_outflow) * invec[-1]

        Matrix, p6 = WriteMatrix(N_nodes, inflows, outflows, z, alpha, Cc,
                                 Nflows, Dd)[0], WriteMatrix(N_nodes, inflows,
                                 outflows, z, alpha, Cc, Nflows, Dd)[1]

# ------Make function
        def dX_dt(sm, t):
            return 1/A * np.squeeze(np.asarray(np.dot(Matrix, sm) + p6))

# ------Solver
        wsol = odeint(dX_dt, in_con, t, atol=1e-08, rtol=1e-06, mxstep=5000000)

# ------Calculate root mean square error
        r = np.sqrt(np.nanmean(((wsol[1:, :].transpose() - Obs)**2)))

        return r

# Need to make sure the optimization preserves total inflow = total outflows

    inflow_frac = inflows[:, 1]/np.sum(inflows[:, 1])
    outflow_frac = outflows[:, 1]/np.sum(outflows[:, 1])

    Ini_tot_flow = np.sum(inflows[:, 1])

    FracLocations = inflows[:, 0] + outflows[:, 0]

# --create parameter array input
    param = [alpha] + inflow_frac.tolist() + outflow_frac.tolist() +\
            FracLocations.tolist() + [Ini_tot_flow]
# --create bounds input
    a = ()
    b = ()
    for i in range(Nflows):
        if inflows[i, 0] == 0:
            a = a + ((0, 0), )
        else:
            a = a + ((0.01, 1), )
        if outflows[i, 0] == 0:
            b = b + ((0, 0), )
        else:
            b = b + ((0.01, 1), )

    parambounds = ((Bound[0], Bound[1]), ) + a + b + FracBounds +\
                    ((Bound[2], Bound[3]), )

    return minimize(CalculateRMSE, param, method=method, bounds=parambounds, **minimise_param)
