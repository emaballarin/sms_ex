#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# IMPORTS:
from typing import Union, Tuple, List, Optional
import math
import matplotlib.pyplot as plt
import numpy as np

# TYPE HINTS:
realnum = Union[int, float]
trinumt = Tuple[realnum, realnum, realnum]


# FUNCTIONS
def sir_log(
    logger: List[List[Optional[realnum]]],
    t: realnum,
    s: realnum,
    i: realnum,
    r: realnum,
) -> None:
    logger[0].append(t)
    logger[1].append(s)
    logger[2].append(i)
    logger[3].append(r)


def dtdoctor(alpha, beta, npop, s, i) -> None:
    dt_sysa = 1.0 / (beta * i)
    dt_sysb = (1.0 / (beta * s - alpha)) * (-1 + npop / i)
    dt_sugg = min(dt_sysa, dt_sysb)
    print(
        "Your timestep should be less than ",
        dt_sugg,
        " or bad things will happen! Look...",
    )


def sir_ode_base(
    x: trinumt,
    tend: realnum = 1.0,
    deltat: realnum = 0.001,
    beta: realnum = 0.1,
    alpha: realnum = 0.01,
    stoch_frac: realnum = 0,
    log: bool = False,
    call_dtdoctor: bool = True,
) -> Tuple[trinumt, List[List[Optional[realnum]]], List[List[Optional[realnum]]]]:

    # Get args and unpack
    s: realnum
    i: realnum
    r: realnum
    s, i, r = x

    # Value checks
    if not (s >= 0 and i >= 0 and r >= 0):
        raise ValueError("ERROR: Values of S, I, R must be non-negative!")
    if stoch_frac < 0 or stoch_frac > 1:
        raise ValueError("ERROR: stoch_frac must be a fraction (i.e. between 0 and 1)!")

    # Internal initialization
    if log:
        logger: List[List[Optional[realnum]]] = [[], [], [], []]
        logger_discretized: List[List[Optional[realnum]]] = [[], [], [], []]
    else:
        logger: List[List[Optional[realnum]]] = [[None], [None], [None], [None]]
        logger_discretized: List[List[Optional[realnum]]] = logger
    npop: realnum = s + i + r
    t: realnum = 0.0
    niter: int = math.ceil(tend / deltat)

    if log:
        sir_log(logger, t, s, i, r)
        sir_log(logger_discretized, t, round(s), round(i), round(r))

    # Repeatedly integrate:
    for itr in range(niter):

        scopy: Optional[realnum]
        icopy: Optional[realnum]

        if call_dtdoctor:
            # Anamnesis for dt Doctor :)
            scopy, icopy = s, i
        else:
            # Make the analyzer happy
            scopy, icopy = None, None

        s, i = (
            s + deltat * (-beta * s * i),
            i + deltat * (beta * s * i - alpha * i),
        )
        r: realnum = npop - s - i
        # Sanity check
        if not (
            (s >= 0 and i >= 0 and r >= 0) and (s <= npop and i <= npop and r <= npop)
        ):
            if call_dtdoctor:
                dtdoctor(alpha, beta, npop, scopy, icopy)
            raise RuntimeError(
                "ERROR: Physical constraint violation. Try lowering the timestep!"
            )

        # Stochastic noise addition, estimated as Gaussian (mean +/- 3-sigma) ~ (-inf, +inf)
        sinoise = np.random.randn() * min(s, i) * stoch_frac / 3.0
        irnoise = np.random.randn() * min(i, r) * stoch_frac / 3.0
        s, i, r = s + sinoise, i - sinoise + irnoise, r - irnoise

        t: realnum = (1 + itr) * deltat

        if log:
            sir_log(logger, t, s, i, r)
            sir_log(logger_discretized, t, round(s), round(i), round(r))

    return (s, i, r), logger, logger_discretized


def sir_ode(
    x: trinumt,
    tend: realnum = 1.0,
    deltat: realnum = 0.001,
    beta: realnum = 0.1,
    alpha: realnum = 0.01,
    stoch_frac: realnum = 0,
    call_dtdoctor: bool = True,
) -> trinumt:
    ret, _, _ = sir_ode_base(
        x, tend, deltat, beta, alpha, stoch_frac, False, call_dtdoctor
    )
    return ret


def sir_plot(
    x: trinumt,
    tend: realnum = 1.0,
    deltat: realnum = 0.001,
    beta: realnum = 0.1,
    alpha: realnum = 0.01,
    stoch_frac: realnum = 0,
    call_dtdoctor: bool = True,
    discretize: bool = True,
) -> None:
    _, cont, disc = sir_ode_base(
        x, tend, deltat, beta, alpha, stoch_frac, True, call_dtdoctor
    )
    if discretize:
        data = disc
    else:
        data = cont

    x = np.array(data[0])
    ys = np.array(data[1])
    yi = np.array(data[2])
    yr = np.array(data[3])

    plt.plot(x, ys, "b")
    plt.plot(x, yi, "r")
    plt.plot(x, yr, "g")
    plt.show()

    return None


# TESTS:
sir_plot((1000, 2, 0), 500, 0.008, stoch_frac=0.05, discretize=True, call_dtdoctor=True)
