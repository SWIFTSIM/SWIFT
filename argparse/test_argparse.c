#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "argparse.h"

static const char *const usages[] = {
    "test_argparse [options] [[--] args]",
    "test_argparse [options]",
    NULL,
};

#define PERM_READ (1 << 0)
#define PERM_WRITE (1 << 1)
#define PERM_EXEC (1 << 2)

struct stuff {
  const char *path[10];
  int npath;
};

static int callback(struct argparse *self, const struct argparse_option *opt) {
  printf("Called back... %s\n", *(char **)opt->value);
  struct stuff *data = (struct stuff *)opt->data;
  data->path[data->npath] = *(char **)opt->value;
  data->npath++;
  return 1;
}

int main(int argc, const char **argv) {
  int force = 0;
  int self_gravity = 0;
  int int_num = 0;
  float flt_num = 0.f;
  struct stuff data;
  data.npath = 0;
  data.path[0] = NULL;
  const char *buffer;
  int perms = 0;
  int npath;

  struct argparse_option options[] = {
      OPT_HELP(),
      OPT_GROUP("Basic options"),
      OPT_BOOLEAN('f', "force", &force, "force to do", NULL, 0, 0),
      OPT_BOOLEAN(0, "self-gravity", &self_gravity, "use self gravity", NULL, 0,
                  0),
      OPT_STRING('P', "path", &buffer, "path to read", &callback,
                 (intptr_t)&data, 0),
      OPT_INTEGER('i', "int", &int_num, "selected integer", NULL, 0, 0),
      OPT_FLOAT('s', "float", &flt_num, "selected float", NULL, 0, 0),
      OPT_END(),
  };

  struct argparse argparse;
  argparse_init(&argparse, options, usages, 0);
  argparse_describe(
      &argparse,
      "\nA brief description of what the program does and how it works.",
      "\nAdditional description of the program after the description of the "
      "arguments.");
  argc = argparse_parse(&argparse, argc, argv);
  if (force != 0) printf("force: %d\n", force);
  if (self_gravity != 0) printf("self_gravity: %d\n", self_gravity);
  if (data.npath > 0) {
    for (int i = 0; i < data.npath; i++) printf("path: %s\n", data.path[i]);
  }
  if (int_num != 0) printf("int_num: %d\n", int_num);
  if (flt_num != 0) printf("flt_num: %g\n", flt_num);
  if (argc != 0) {
    printf("argc: %d\n", argc);
    int i;
    for (i = 0; i < argc; i++) {
      printf("argv[%d]: %s\n", i, *(argv + i));
    }
  }
  if (perms) {
    printf("perms: %d\n", perms);
  }
  return 0;
}
