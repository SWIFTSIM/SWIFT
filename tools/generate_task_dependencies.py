#!/usr/bin/env python3
import sys
import re


def get_string_between(s, beg, end):
    regex = beg + ".*?" + end
    m = re.search(regex, s)

    if m is None:
        return None

    m = m.group(0)
    return m[len(beg):-len(end)]


class Task:
    def __init__(self, s):

        self.name = s[:s.find("[")]
        self.iact_type = self.read_attribute("type", s)
        self.active = self.read_attribute("active", s)
        self.level = self.read_attribute("level", s)
        self.updates = self.read_attribute("updates", s)
        self.policy = self.read_attribute("policy", s)
        self.cell_contains = self.read_attribute("cell_contains", s)
        self.act_on = self.read_attribute("act_on", s)
        self.link = []

    def read_attribute(self, attribute_name, s):
        attribute_name += "="
        attr = get_string_between(s, attribute_name, ",")
        if attr is not None:
            return attr

        attr = get_string_between(s, attribute_name, "]")
        return attr

    def add_link(self, name):
        self.link.append(name)

    def print_task(self):
        print("{}: iact={}, active={}, level={}, link={}".format(
            self.name, self.iact_type, self.active,
            self.level, self.link))


class Reader:
    def __init__(self, filename):
        self.filename = filename
        self.tasks = {}
        self.name = None
        self.read()

    def read(self):
        with open(self.filename, "r") as f:
            for line in f:
                line = line.rstrip()
                self.read_line(line)

    def read_line(self, line):
        if "->" in line:
            self.read_link(line)
        elif "[" in line and "]" in line:
            self.read_definition(line)
        elif "label" in line:
            s = line.find('"')
            e = line.rfind('"')
            self.name = line[s+1:e]

    def read_link(self, line):
        names = line.split("->")
        if names[0] not in self.tasks:
            raise Exception(
                "Trying to link {} without any definition".format(
                    names[0]))
        self.tasks[names[0]].add_link(names[1])

    def read_definition(self, line):
        t = Task(line)
        self.tasks[t.name] = t

    def print_reader(self):
        print("Simulation type:", self.name)
        for t in self.tasks:
            self.tasks[t].print_task()


if __name__ == "__main__":
    reader = Reader(sys.argv[-1])

    reader.print_reader()

