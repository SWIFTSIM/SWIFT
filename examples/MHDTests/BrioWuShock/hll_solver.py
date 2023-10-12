#!/usr/bin/env python3

"""
hll_solver.py

HLL Riemann solver for hydro and MHD problems.
In standalone mode, this script runs 1 of 3 pre-configured tests and makes a
plot.
As a module, it exposes the solve_MHD_Riemann_problem() function, which can
be used to import the solution to a Riemann problem into another script.
"""

import numpy as np

# data types for various types of vectors
# this allows us to use structured arrays that can be indexed by name rather
# than number
state_dtype = [
    ("rho", np.float64),
    ("rho_u", np.float64),
    ("rho_v", np.float64),
    ("rho_w", np.float64),
    ("Bx", np.float64),
    ("By", np.float64),
    ("Bz", np.float64),
    ("e", np.float64),
]
eigen_dtype = [
    ("l1", np.float64),
    ("l2", np.float64),
    ("l3", np.float64),
    ("l4", np.float64),
    ("l5", np.float64),
    ("l6", np.float64),
    ("l7", np.float64),
]
prim_dtype = [
    ("rho", np.float64),
    ("u", np.float64),
    ("v", np.float64),
    ("w", np.float64),
    ("p", np.float64),
    ("Bx", np.float64),
    ("By", np.float64),
    ("Bz", np.float64),
]


def get_state(gamma, primitives):
    """
    Convert a dictionary of primitive variables into a state vector
    Only useful to generate initial conditions, where providing an input
    dictionary actually makes sense.
    """
    state = np.zeros(1, dtype=state_dtype)
    state["rho"] = primitives["rho"]
    state["rho_u"] = primitives["rho"] * primitives["u"]
    state["rho_v"] = primitives["rho"] * primitives["v"]
    state["rho_w"] = primitives["rho"] * primitives["w"]
    state["Bx"] = primitives["Bx"]
    state["By"] = primitives["By"]
    state["Bz"] = primitives["Bz"]
    v2 = primitives["u"] ** 2 + primitives["v"] ** 2 + primitives["w"] ** 2
    B2 = primitives["Bx"] ** 2 + primitives["By"] ** 2 + primitives["Bz"] ** 2
    state["e"] = (
        primitives["p"] / (gamma - 1.0) + 0.5 * primitives["rho"] * v2 + 0.5 * B2
    )
    return state


def get_states(gamma, primitives):
    """
    Convert a structured array of primitive variables into a state vector
    """
    state = np.zeros(primitives.shape, dtype=state_dtype)
    state["rho"] = primitives["rho"]
    state["rho_u"] = primitives["rho"] * primitives["u"]
    state["rho_v"] = primitives["rho"] * primitives["v"]
    state["rho_w"] = primitives["rho"] * primitives["w"]
    state["Bx"] = primitives["Bx"]
    state["By"] = primitives["By"]
    state["Bz"] = primitives["Bz"]
    v2 = primitives["u"] ** 2 + primitives["v"] ** 2 + primitives["w"] ** 2
    B2 = primitives["Bx"] ** 2 + primitives["By"] ** 2 + primitives["Bz"] ** 2
    state["e"] = (
        primitives["p"] / (gamma - 1.0) + 0.5 * primitives["rho"] * v2 + 0.5 * B2
    )
    return state


def get_primitives(gamma, state):
    """
    Convert an array of state vectors into an array of primitive variables.
    """
    prim = np.zeros(state.shape, dtype=prim_dtype)
    prim["rho"] = state["rho"]
    prim["u"] = state["rho_u"] / state["rho"]
    prim["v"] = state["rho_v"] / state["rho"]
    prim["w"] = state["rho_w"] / state["rho"]
    B2 = state["Bx"] ** 2 + state["By"] ** 2 + state["Bz"] ** 2
    rhov2 = state["rho_u"] ** 2 + state["rho_v"] ** 2 + state["rho_w"] ** 2
    prim["p"] = (gamma - 1.0) * (state["e"] - 0.5 * rhov2 / state["rho"] - 0.5 * B2)
    prim["Bx"] = state["Bx"]
    prim["By"] = state["By"]
    prim["Bz"] = state["Bz"]
    return prim


def get_eigenvalues(gamma, state):
    """
    Get the eigenvalues for the given array of state vectors.
    """
    u = state["rho_u"] / state["rho"]
    B2 = state["Bx"] ** 2 + state["By"] ** 2 + state["Bz"] ** 2
    rhov2 = state["rho_u"] ** 2 + state["rho_v"] ** 2 + state["rho_w"] ** 2
    p = (gamma - 1.0) * (state["e"] - 0.5 * rhov2 / state["rho"] - 0.5 * B2)
    if (p < 0.0).any():
        print(f"state: {state}")
        raise RuntimeWarning("Pressure is negative!")
    c_a = np.abs(state["Bx"]) / np.sqrt(state["rho"])
    try:
        c_f = np.sqrt(
            0.5
            * (
                gamma * p
                + B2
                + np.sqrt((gamma * p + B2) ** 2 - 4.0 * gamma * p * state["Bx"] ** 2)
            )
            / state["rho"]
        )
    except RuntimeWarning as w:
        print(f"p: {p}, B2: {B2}, Bx: {state['Bx']}, rho: {state['rho']}")
        raise w
    try:
        c_s = np.sqrt(
            0.5
            * (
                gamma * p
                + B2
                - np.sqrt((gamma * p + B2) ** 2 - 4.0 * gamma * p * state["Bx"] ** 2)
            )
            / state["rho"]
        )
    except RuntimeWarning as w:
        print(f"p: {p}, B2: {B2}, Bx: {state['Bx']}, rho: {state['rho']}")
        raise w
    eigenvalues = np.zeros(state.shape, dtype=eigen_dtype)
    eigenvalues["l1"] = u - c_f
    eigenvalues["l2"] = u - c_a
    eigenvalues["l3"] = u - c_s
    eigenvalues["l4"] = u
    eigenvalues["l5"] = u + c_s
    eigenvalues["l6"] = u + c_a
    eigenvalues["l7"] = u + c_f
    assert (eigenvalues["l1"] <= eigenvalues["l2"]).all()
    assert (eigenvalues["l2"] <= eigenvalues["l3"]).all()
    assert (eigenvalues["l3"] <= eigenvalues["l4"]).all()
    assert (eigenvalues["l4"] <= eigenvalues["l5"]).all()
    assert (eigenvalues["l5"] <= eigenvalues["l6"]).all()
    assert (eigenvalues["l6"] <= eigenvalues["l7"]).all()
    return eigenvalues


def get_minmax(eigenvalues):
    """
    Get the minimum and maximum eigenvalues in the given array of eigenvalues.
    """
    eigenview = eigenvalues.view(dtype=np.float64).reshape((*eigenvalues.shape, -1))
    lmin = np.nanmin(eigenview, axis=-1, keepdims=True)
    lmax = np.nanmax(eigenview, axis=-1, keepdims=True)
    return lmin, lmax


def get_flux(gamma, state):
    """
    Get the flux vector(s) for the given array of state vectors.
    """
    u = state["rho_u"] / state["rho"]
    v = state["rho_v"] / state["rho"]
    w = state["rho_w"] / state["rho"]
    B2 = state["Bx"] ** 2 + state["By"] ** 2 + state["Bz"] ** 2
    rhov2 = state["rho_u"] ** 2 + state["rho_v"] ** 2 + state["rho_w"] ** 2
    p = (gamma - 1.0) * (state["e"] - 0.5 * rhov2 / state["rho"] - 0.5 * B2)
    pT = p + 0.5 * B2
    flux = np.zeros(state.shape, dtype=state_dtype)
    flux["rho"] = state["rho_u"]
    flux["rho_u"] = state["rho_u"] * u + pT - state["Bx"] ** 2
    flux["rho_v"] = state["rho_v"] * u - state["Bx"] * state["By"]
    flux["rho_w"] = state["rho_w"] * u - state["Bx"] * state["Bz"]
    flux["By"] = state["By"] * u - state["Bx"] * v
    flux["Bz"] = state["Bz"] * u - state["Bx"] * w
    flux["e"] = (state["e"] + pT) * u - state["Bx"] * (
        u * state["Bx"] + v * state["By"] + w * state["Bz"]
    )
    return flux


class RiemannSolver:
    """
    General interface for Riemann solvers.
    """

    def __init__(self, gamma):
        """
        Constructor.
        """
        self._gamma = gamma


class HLLRiemannSolver(RiemannSolver):
    """
    HLL solver.
    """

    def __init__(self, gamma):
        """
        Constructor.
        """
        super().__init__(gamma)

    def solve(self, left_state, right_state):
        """
        Solve the Riemann problem with given left and right state(s) using the
        HLL approximate solver for the flux.
        """
        Leigen = get_eigenvalues(self._gamma, left_state)
        Reigen = get_eigenvalues(self._gamma, right_state)
        Llmin, Llmax = get_minmax(Leigen)
        Rlmin, Rlmax = get_minmax(Reigen)
        SL = np.minimum(Llmin, Rlmin)
        SR = np.maximum(Llmax, Rlmax)
        vmax = np.maximum(np.abs(SL), np.abs(SR)).max()
        SL = np.minimum(SL, 0.0)
        SR = np.maximum(SR, 0.0)
        FL = get_flux(self._gamma, left_state)
        FR = get_flux(self._gamma, right_state)
        Fstar = np.zeros(FL.shape, dtype=state_dtype)
        Fstarview = Fstar.view(np.float64).reshape((*Fstar.shape, -1))
        FLview = FL.view(np.float64).reshape((*FL.shape, -1))
        FRview = FR.view(np.float64).reshape((*FR.shape, -1))
        UL = left_state.view(np.float64).reshape((*left_state.shape, -1))
        UR = right_state.view(np.float64).reshape((*right_state.shape, -1))
        Fstarview[:] = (SR * FLview - SL * FRview + SR * SL * (UR - UL)) / (SR - SL)
        return Fstar, vmax


def solve_MHD_Riemann_problem(riemann_problem_dict, verbose=True):
    """
    Solve the MHD Riemann problem with the given specification.

    The input dictionary should contain:
     - left_state/right_state: input left and right state, given as dictionaries
        containing the primitive variables rho, u, v, w, p, Bx, By and Bz.
     - gamma: the adiabatic index
     - ncell: the number of cells to use
     - xsize: the size of the spatial domain
     - time: the desired time of the solution
     - dt: the time step to use (fixed)

    Optionally, the Riemann problem can be (explicitly) solved in co-moving
    coordinates. This requires the following additional parameters:
     - comoving: True for co-moving integration
     - a_from_t: a function that returns the scale factor corresponding to a
        given time since the start of the calculation. This depends on the
        cosmology.
     - comoving_term: a function of the signature
        f(a, X)
       that returns A in the following equation:
        dV/dt = A*V
       if V depends on the scale factor with a power X, e.g. for V=rho,
        rho' = rho*a^3
       and X = 3.
       this term expresses the change in variable V because of the expansion of
       the Universe and depends on the cosmology.
    Note that we always express the variables in a frame that moves with the
    expansion of the Universe, so we (probably?) do not explicitly need to take
    into account the Hubble flow.
    We also assume for convenience that the time variable is the re-scaled time
    for which dt' = dt/a^2.
    """

    # initial conditions and problem setup
    gamma = riemann_problem_dict["gamma"]
    IC_L = get_state(gamma, riemann_problem_dict["left_state"])
    IC_R = get_state(gamma, riemann_problem_dict["right_state"])

    Ncell = riemann_problem_dict["ncell"]
    xmin = -0.5 * riemann_problem_dict["xsize"]
    xmax = 0.5 * riemann_problem_dict["xsize"]
    dx = (xmax - xmin) / Ncell
    xs = np.linspace(0.0, riemann_problem_dict["xsize"], Ncell)
    xs -= xs.mean()
    dt = riemann_problem_dict["dt"]
    tend = riemann_problem_dict["time"]
    Nstep = int(tend / dt)
    solver = HLLRiemannSolver(gamma)

    comoving = (
        riemann_problem_dict["comoving"]
        if "comoving" in riemann_problem_dict
        else False
    )
    if comoving:
        a_from_t = riemann_problem_dict["a_from_t"]
        comoving_term = riemann_problem_dict["comoving_term"]
        H_0 = riemann_problem_dict["H_0"]

    # initialise state vectors
    states = np.zeros(Ncell, dtype=state_dtype)
    states[: Ncell // 2] = IC_L
    states[Ncell // 2 :] = IC_R

    if comoving:
        # convert from comoving to physical coordinates
        prim = get_primitives(gamma, states)
        a = a_from_t(0.0)
        prim["rho"] /= a ** 3
        prim["u"] /= a
        prim["v"] /= a
        prim["w"] /= a
        """
        prim["u"] -= H_0 * xs * a
        prim["v"] -= H_0 * xs * a
        prim["w"] -= H_0 * xs * a
        """
        prim["p"] /= a ** (3 * gamma)
        prim["Bx"] /= a ** 2
        prim["By"] /= a ** 2
        prim["Bz"] /= a ** 2
        states = get_states(gamma, prim)

    # time integration loop
    warning = None
    for istep in range(Nstep):

        left = states[:-1]
        right = states[1:]

        try:
            flux, vmax = solver.solve(left, right)
        except RuntimeWarning as w:
            print("Problem!")
            warning = w
            break

        if comoving:
            # get the cell size at the current cosmological time
            a = a_from_t(istep * dt)
            stepdx = dx * a
            # use actual time rather than rescaled time for the integration
            stepdt = dt * a ** 2
        else:
            stepdx = dx
            stepdt = dt

        dtmax = 0.8 * stepdx / vmax
        if comoving:
            # convert the maximum dt to rescaled time
            dtmax /= a ** 2
        if verbose:
            print(f"step {istep+1} out of {Nstep}, dtmax: {dtmax}", end="\r")
        if dtmax < dt:
            print("Timestep too large!")
            break

        stateview = states.view(np.float64).reshape((*states.shape, -1))
        fluxview = flux.view(np.float64).reshape((*flux.shape, -1))
        stateview[1:-1] += (stepdt / stepdx) * fluxview[:-1]
        stateview[1:-1] -= (stepdt / stepdx) * fluxview[1:]

        if comoving:
            # add cosmological expansion terms
            # we could use exp() for an exact solution to the first order ODE,
            # but since dt is very small anyway, that makes little difference
            prim = get_primitives(gamma, states)
            prim["rho"] *= 1.0 + dt * comoving_term(a, 3)
            prim["u"] *= 1.0 + dt * comoving_term(a, 1)
            """
            prim["u"] -= H_0 * stepdx
            """
            prim["v"] *= 1.0 + dt * comoving_term(a, 1)
            prim["w"] *= 1.0 + dt * comoving_term(a, 1)
            prim["p"] *= 1.0 + dt * comoving_term(a, 3 * gamma)
            prim["Bx"] *= 1.0 + dt * comoving_term(a, 2)
            prim["By"] *= 1.0 + dt * comoving_term(a, 2)
            prim["Bz"] *= 1.0 + dt * comoving_term(a, 2)
            states = get_states(gamma, prim)

    if verbose:
        print(f"\nDone {Nstep} steps.")

    primitives = get_primitives(gamma, states)
    if comoving:
        # convert back to co-moving variables
        a = a_from_t(Nstep * dt)
        primitives["rho"] *= a ** 3
        primitives["u"] *= a
        primitives["v"] *= a
        primitives["w"] *= a
        """
        primitives["u"] += H_0 * xs * a
        primitives["v"] += H_0 * xs * a
        primitives["w"] += H_0 * xs * a
        """
        primitives["p"] *= a ** (3 * gamma)
        primitives["Bx"] *= a ** 2
        primitives["By"] *= a ** 2
        primitives["Bz"] *= a ** 2

    return xs, primitives


def plot_states(x, primitives, t, outputfile):
    """
    Plot the given primitive variable states to a file with a given name.
    """

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as pl

    fig, ax = pl.subplots(3, 3, sharex=True)

    ax[0][0].plot(x, primitives["rho"])
    ax[0][0].set_ylabel("rho")
    ax[0][1].plot(x, primitives["p"])
    ax[0][1].set_ylabel("p")
    ax[0][2].axis("off")
    ax[1][0].plot(x, primitives["u"])
    ax[1][0].set_ylabel("u")
    ax[1][1].plot(x, primitives["v"])
    ax[1][1].set_ylabel("v")
    ax[1][2].plot(x, primitives["w"])
    ax[1][2].set_ylabel("w")
    ax[2][0].plot(x, primitives["Bx"])
    ax[2][0].set_ylabel("Bx")
    ax[2][1].plot(x, primitives["By"])
    ax[2][1].set_ylabel("By")
    ax[2][2].plot(x, primitives["Bz"])
    ax[2][2].set_ylabel("Bz")

    ax[0][1].set_title(f"t={t}")
    ax[2][0].set_xlabel("x")
    ax[2][1].set_xlabel("x")
    ax[2][2].set_xlabel("x")

    #    ax[2][0].set_xlim(-0.05, 0.05)
    pl.tight_layout()
    pl.savefig(outputfile, dpi=300)


if __name__ == "__main__":
    """
    Standalone mode.

    Runs a test of choice (options: DaiWoodward or Sod) and creates an image
    with the given name.

    Optional parameters allow changing the time, time step size, size of the
    spatial domain, number of cells and adiabatic index.

    This assumes ordinary (non co-moving) integration.
    """

    import argparse

    argumentparser = argparse.ArgumentParser()
    argumentparser.add_argument("test", choices=["DaiWoodward", "Sod", "BrioWu"])
    argumentparser.add_argument("output")
    argumentparser.add_argument("--time", "-t", type=float, default=0.2)
    argumentparser.add_argument("--dt", "-d", type=float, default=1.0e-4)
    argumentparser.add_argument("--xsize", "-x", type=float, default=1.0)
    argumentparser.add_argument("--ncell", "-n", type=int, default=800)
    argumentparser.add_argument("--gamma", "-g", type=float, default=5.0 / 3.0)
    args = argumentparser.parse_args()

    gamma = 5.0 / 3.0
    if args.test == "DaiWoodward":
        left_state = {
            "rho": 1.08,
            "p": 0.95,
            "u": 1.2,
            "v": 0.01,
            "w": 0.5,
            "Bx": 4.0 / np.sqrt(4.0 * np.pi),
            "By": 3.6 / np.sqrt(4.0 * np.pi),
            "Bz": 2.0 / np.sqrt(4.0 * np.pi),
        }
        right_state = {
            "rho": 1.0,
            "p": 1.0,
            "u": 0.0,
            "v": 0.0,
            "w": 0.0,
            "Bx": 4.0 / np.sqrt(4.0 * np.pi),
            "By": 4.0 / np.sqrt(4.0 * np.pi),
            "Bz": 2.0 / np.sqrt(4.0 * np.pi),
        }
    elif args.test == "Sod":
        left_state = {
            "rho": 1.0,
            "p": 1.0,
            "u": 0.0,
            "v": 0.0,
            "w": 0.0,
            "Bx": 0.0,
            "By": 0.0,
            "Bz": 0.0,
        }
        right_state = {
            "rho": 0.125,
            "p": 0.1,
            "u": 0.0,
            "v": 0.0,
            "w": 0.0,
            "Bx": 0.0,
            "By": 0.0,
            "Bz": 0.0,
        }
    elif args.test == "BrioWu":
        left_state = {
            "rho": 1.0,
            "p": 5.0 / 3.0 - 1.0,
            "u": 0.0,
            "v": 0.0,
            "w": 0.0,
            "Bx": 0.75,
            "By": 1.0,
            "Bz": 0.0,
        }
        right_state = {
            "rho": 0.125,
            "p": 0.1 * (5.0 / 3.0 - 1.0),
            "u": 0.0,
            "v": 0.0,
            "w": 0.0,
            "Bx": 0.75,
            "By": -1.0,
            "Bz": 0.0,
        }

    riemann_problem = {
        "left_state": left_state,
        "right_state": right_state,
        "time": args.time,
        "dt": args.dt,
        "xsize": args.xsize,
        "gamma": args.gamma,
        "ncell": args.ncell,
        "solver": "HLL",
    }

    x, solution = solve_MHD_Riemann_problem(riemann_problem)

    plot_states(x, solution, args.time, args.output)
