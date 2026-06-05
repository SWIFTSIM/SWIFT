###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Katy Proctor (katy.proctor@fysik.su.se)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################

import re
from collections import Counter

# Matches the step header lines, e.g. "       7   4.272461e-05 ..."
STEP_LINE_RE = re.compile(r"^\s+(\d+)\s+\S+\s+\S+\s+\S+\s+\S+", re.MULTILINE)
PID_RE = re.compile(r"\[PID[ij](\d+)\]")


def split_into_timesteps(log_content):
    """
    Split the log into per-timestep blocks.
    """
    matches = list(STEP_LINE_RE.finditer(log_content))
    blocks = []
    for i, m in enumerate(matches):
        step_num = int(m.group(1))
        start = m.end()  # content starts after this step's header line
        end = matches[i + 1].start() if i + 1 < len(matches) else len(log_content)
        block = log_content[start:end]
        blocks.append((step_num, block))
    return blocks


def check_timestep(step_num, block):
    """
    Within a single timestep, check that no PID appears more than twice.
    """
    pids = PID_RE.findall(block)
    if not pids:
        return True  # no interactions this step

    counts = Counter(pids)
    more_than_two = {pid: c for pid, c in counts.items() if c > 2}

    total_unique = len(counts)
    exactly_two = sum(1 for c in counts.values() if c == 2)

    print(f"\n--- Step {step_num} ---")
    print(f"  Total unique PIDs: {total_unique}")
    print(f"  PIDs appearing exactly twice (OK): {exactly_two}")
    print(f"  PIDs appearing more than twice (VIOLATION): {len(more_than_two)}")

    if more_than_two:
        print("  Violating PIDs:")
        for pid, count in sorted(more_than_two.items(), key=lambda x: -x[1]):
            print(f"    PID {pid}: {count} times")
        return False

    return True


def count_pids(log_content):
    blocks = split_into_timesteps(log_content)

    if not blocks:
        print("No timestep headers found in log.")
        return

    steps_with_interactions = 0
    violation_steps = []

    for step_num, block in blocks:
        pids = PID_RE.findall(block)
        if not pids:
            continue  # skip quiet steps
        steps_with_interactions += 1
        ok = check_timestep(step_num, block)
        if not ok:
            violation_steps.append(step_num)

    print(f"\n{'='*50}")
    print(f"Steps with interactions: {steps_with_interactions}")
    print(f"Steps with violations:   {len(violation_steps)}")
    if violation_steps:
        print(f"Violation steps: {violation_steps}")
    else:
        print("No violations found — every particle scattered at most once per step.")


with open("output.log", "r") as f:
    log_content = f.read()

count_pids(log_content)
