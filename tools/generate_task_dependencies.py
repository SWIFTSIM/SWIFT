#!/usr/bin/env python3
import sys
import re
import typing
from task import tasks as ref_tasks


def get_string_between(s: str, beg: str, end: str):
    """
    Match the smallest string contained within
    the substring beg and end. Returns None if nothing
    is found
    """
    regex = beg + ".*?" + end
    m = re.search(regex, s)

    if m is None:
        return None

    m = m.group(0)
    return m[len(beg):-len(end)]


class Task:
    """
    Class representing a single task along with its links.
    """
    def __init__(self, s: str):
        """
        Initialize the task from a line in the input file.
        """

        self.name = s[:s.find("[")]
        self.iact_type = self.read_attribute("type", s)
        self.active = self.read_attribute("active", s)
        self.level = self.read_attribute("level", s)
        self.updates = self.read_attribute("updates", s)
        self.policy = self.read_attribute("policy", s)
        self.cell_contains = self.read_attribute("cell_contains", s)
        self.act_on = self.read_attribute("act_on", s)
        self.link = []

    def read_attribute(self, attribute_name: str, s: str):
        """
        Read a single attribute in the definition line.
        """
        attribute_name += "="
        attr = get_string_between(s, attribute_name, ",")
        if attr is not None:
            return attr

        attr = get_string_between(s, attribute_name, "]")
        return attr

    def add_link(self, name: str):
        """
        Add a link for this task.
        """
        self.link.append(name)

    def print_task(self):
        """
        Print the task.
        """
        print("{}: iact={}, active={}, level={}, link={}".format(
            self.name, self.iact_type, self.active,
            self.level, self.link))

    def write_if_condition(self, f: typing.TextIO, level: bool = True,
                           policy: bool = True, task_type: bool = False,
                           task_subtype: bool = False):
        """
        Write the if condition for this task.
        The conditions can be disabled if needed.

        Returns
        -------

        written: bool
            Did this function write something (e.g. empty if condition)?
        """
        if_condition = ""
        if self.level and level:
            if_condition += "(c->%s == c) && " % self.level
        if self.policy and policy:
            if_condition += "(e->policy && engine_policy_%s) && " % self.policy
        if task_type:
            if_condition += "(t->type == task_type_%s) && " % \
                ref_tasks[self.name]["type"]
        if task_subtype:
            if_condition += "(t->subtype == task_subtype_%s) && " % \
                ref_tasks[self.name]["subtype"]

        if len(if_condition) == 0:
            return False

        # Remove ' &&'
        if_condition = if_condition[:-4]

        if_condition = "if (%s) {\n" % if_condition
        f.write(if_condition)

        return True

    def write_maketask_definitions(self, f: typing.TextIO):
        """
        Write the code corresponding to the definition of the task
        (e.g. scheduler_addtask).
        It creates the tasks acting on a single cell at a time and
        does not care about dependencies.
        This will replace the work done in engine_make_hierarchical_tasks.
        """
        is_implicit = int(self.iact_type == "implicit")
        if self.iact_type != "single" and not is_implicit:
            print("Skipping", self.name)
            return

        # Compute the scheduler_addtask
        task_type = ref_tasks[self.name]["type"]
        creation = "c->{name} = scheduler_addtask(s, task_type_{task_type},"
        creation += " task_subtype_none, 0, {implicit}, c, NULL)"
        creation = creation.format(
            name=ref_tasks[self.name]["variable"], task_type=task_type,
            implicit=is_implicit)

        # Write the code inside a file
        if_done = self.write_if_condition(f)
        code = "\t %s;\n" % creation
        if if_done:
            code += "};\n"
        f.write(code)

    def write_maketask_single_single(
            self, f: typing.TextIO, link: str, tasks: dict):
        # Generate part of the variable name
        self_level = ""
        if self.level:
            self_level = "%s->" % self.level

        # Grab data for the link
        link = tasks[link]
        link_level = ""
        if link.level:
            link_level = "%s->" % link.level

        # Write the if condition
        write_policy = self.policy != link.policy
        if_done2 = link.write_if_condition(
            f, level=False, policy=write_policy)

        # Generate the scheduler_addunlock
        unlock = "\t\tscheduler_addunlock(s, ci->{level1}{name1},"
        unlock += " ci->{level2}{name2});\n"
        unlock = unlock.format(
            level1=self_level, name1=ref_tasks[self.name]["variable"],
            level2=link_level, name2=ref_tasks[link.name]["variable"]
        )

        # Close parenthesis
        if if_done2:
            unlock += "\t};\n"
        f.write(unlock)

    def write_maketask_iact_iact(
            self, f: typing.TextIO, link: str, tasks: dict):
        raise Exception("It should never happens")

    def write_maketask_iact_single(
            self, f: typing.TextIO, link: str, tasks: dict):
        return

    def write_maketask_single_iact(
            self, f: typing.TextIO, link: str, tasks: dict):
        return

    def write_maketask_extra_loop(
            self, f: typing.TextIO, tasks: dict):
        """
        Write the code corresponding to the link from this task
        (e.g. scheduler_addunlock)
        and defines it if it is a task acting on two cells.
        In order to have enough information for the links,
        a dictionary containing all the tasks is required.
        This will replace the work done in
        engine_make_extra_hydroloop_tasks.
        """

        # Write header
        f.write("// %s\n" % self.name)
        if len(self.link) == 0:
            return

        # Write the initial if condition
        if_done = self.write_if_condition(
            f, level=False, task_type=True)

        iact1 = self.iact_type

        # loop over all the links
        for link in self.link:
            if "self" in link or "pair" in link:
                print("Skipping self/pair")
                continue

            iact2 = tasks[link].iact_type

            # Deal with the case where both tasks are not iact
            if iact1 != "iact" and iact2 != "iact":
                self.write_maketask_single_single(f, link, tasks)
            # Deal with the case where both are iact
            elif iact1 == "iact" and iact2 == "iact":
                self.write_maketask_iact_iact(f, link, tasks)
            elif iact1 == "iact":
                self.write_maketask_iact_single(f, link, tasks)
            else:
                self.write_maketask_single_iact(f, link, tasks)

        # Close parenthesis
        if if_done:
            f.write("};\n")


class Reader:
    """
    Class dealing with the task system.
    """
    def __init__(self, filename: str):
        """
        Initialize the reader from a .task file
        """
        self.filename = filename
        self.tasks = {}
        self.name = None
        self.read()

    def read(self):
        """
        Read the .task file
        """
        with open(self.filename, "r") as f:
            for line in f:
                line = line.rstrip()
                self.read_line(line)

    def read_line(self, line: str):
        """
        Read a single line in the .task file
        """
        if "->" in line:
            self.read_link(line)
        elif "[" in line and "]" in line:
            self.read_definition(line)
        elif "label" in line:
            s = line.find('"')
            e = line.rfind('"')
            self.name = line[s+1:e]

    def read_link(self, line: str):
        """
        Read a line containing a link
        """
        names = line.split("->")
        if names[0] not in self.tasks:
            raise Exception(
                "Trying to link {} without any definition".format(
                    names[0]))
        self.tasks[names[0]].add_link(names[1])

    def read_definition(self, line: str):
        """
        Read a line containing a definition.
        """
        t = Task(line)
        self.tasks[t.name] = t

    def print_reader(self):
        """
        Print the task system
        """
        print("Simulation type:", self.name)
        for t in self.tasks:
            self.tasks[t].print_task()

    def write_maketask(self, filename: str):
        """
        Write the code corresponding to the task system into a file.
        """
        with open(filename, "w") as f:
            f.write("// Hierarchical taks\n")
            for name in self.tasks:
                self.tasks[name].write_maketask_definitions(f)

            f.write("// Dependencies\n")
            for name in self.tasks:
                self.tasks[name].write_maketask_extra_loop(f, self.tasks)


if __name__ == "__main__":
    reader = Reader(sys.argv[-1])

    reader.write_maketask("test.c")
