"""Integer time-line helpers, a line by line port of src/timeline.h.

The figures of the time-stepping documentation import these functions so that
the plots use exactly the definitions of the code.
"""

NUM_TIME_BINS = 56
MAX_NR_TIMESTEPS = 1 << (NUM_TIME_BINS + 1)

# Okabe-Ito colour-blind safe palette.
PALETTE = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
}


def _c_div(a, b):
    """Integer division that truncates towards zero, as in C."""
    q = abs(a) // abs(b)
    return q if (a >= 0) == (b >= 0) else -q


def get_integer_timestep(bin_):
    """Length in integer ticks of the time-step of a time bin."""
    if bin_ <= 0:
        return 0
    return 1 << (bin_ + 1)


def get_time_bin(time_step):
    """Time bin of an integer time-step: floor(log2(time_step)) - 1."""
    return time_step.bit_length() - 2


def get_timestep(bin_, time_base):
    """Physical length of the time-step of a time bin."""
    return get_integer_timestep(bin_) * time_base


def get_integer_time_begin(ti_current, bin_):
    """Start of the step of a bin, as seen from ti_current (C semantics)."""
    dti = get_integer_timestep(bin_)
    if dti == 0:
        return 0
    return dti * _c_div(ti_current - 1, dti)


def get_integer_time_end(ti_current, bin_):
    """End of the step of a bin, as seen from ti_current."""
    dti = get_integer_timestep(bin_)
    if dti == 0:
        return 0
    mod = ti_current % dti
    if mod == 0:
        return ti_current
    return ti_current - mod + dti


def get_max_active_bin(time):
    """Highest time bin that is active at a point of the time-line."""
    if time == 0:
        return NUM_TIME_BINS
    bin_ = 1
    while not ((1 << (bin_ + 1)) & time):
        bin_ += 1
    return bin_


def get_min_active_bin(ti_current, ti_old):
    """Lowest active bin, from the length of the last global step."""
    return get_max_active_bin(ti_current - ti_old)


def apply_growth_rules(new_dti, old_bin, ti_end, min_ngb_bin=NUM_TIME_BINS,
                       delta_bin=2):
    """The alignment and growth rules of make_integer_timestep().

    Takes the integer step that the physics asks for and returns the integer
    step the particle is allowed to use. Port of src/timestep.h.
    """
    new_bin = get_time_bin(new_dti)
    new_bin = min(new_bin, min_ngb_bin + delta_bin)
    new_dti = get_integer_timestep(new_bin)
    current_dti = get_integer_timestep(old_bin)
    if old_bin > 0:
        new_dti = min(new_dti, 2 * current_dti)
    dti_timeline = MAX_NR_TIMESTEPS
    while new_dti < dti_timeline:
        dti_timeline //= 2
    new_dti = dti_timeline
    if new_dti > current_dti:
        if (MAX_NR_TIMESTEPS - ti_end) % new_dti > 0:
            new_dti = current_dti
    return new_dti
