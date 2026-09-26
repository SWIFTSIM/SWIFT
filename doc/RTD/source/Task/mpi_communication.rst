.. MPI communication in the task graph

.. _task_mpi_communication:
.. highlight:: c

MPI Communication and Foreign Cells
===================================

This page explains how SWIFT moves data between MPI ranks inside the task
graph. It is written for developers who add or change tasks and who need to
know what a rank can trust when it reads a cell that belongs to another rank.

The short version: every rank owns some cells (the *local* cells) and holds
read-only copies of the neighbouring cells owned by other ranks (the *foreign*
cells). The list of the cells that two ranks exchange is called a *proxy*, and
each rank has one proxy for every rank it interacts with. The owner of a cell
sends its data with *send* tasks. The rank that holds the copy fills it with
matching *recv* tasks. Some fields of a foreign
copy are refreshed at different moments than the same fields of the local
cell, so they can be out of date while a task reads them.

All names below refer to the source files in ``src/``. Look up the functions
there for the exact details.


Local and foreign cells
~~~~~~~~~~~~~~~~~~~~~~~

The top-level grid is known by every rank. ``space.cells_top`` holds all the
top-level cells, and each one has a ``nodeID`` that says which rank owns it.
From the point of view of one rank, a cell is *local* when
``c->nodeID == engine_rank`` and *foreign* otherwise.

What is a proxy?
^^^^^^^^^^^^^^^^

A *proxy* (``struct proxy`` in ``proxy.h``) describes the relationship between
this rank and one other rank. Each rank holds one proxy for every rank that
owns at least one top-level cell this rank needs to interact with. It holds no
proxy for the ranks it never interacts with. The proxies of a rank are stored
in ``engine.proxies``, and ``engine.proxy_ind[r]`` gives the proxy for rank
``r``.

A proxy contains two lists of top-level cells:

* ``cells_in``: the foreign cells owned by the other rank that this rank
  receives. This rank keeps a foreign copy of each of them.
* ``cells_out``: the local cells of this rank that the other rank needs. This
  rank sends them.

For each cell in these lists, ``cells_in_type`` and ``cells_out_type`` store why
the cell is there, as the bit flags ``proxy_cell_type_hydro`` and
``proxy_cell_type_gravity``. A cell that is needed for both keeps both flags:
``proxy_addcell_in()`` and ``proxy_addcell_out()`` add the new flag to the
existing entry instead of listing the cell twice. The flags decide which
communication tasks are created for the cell.

Proxies are built by ``engine_makeproxies()`` (``engine_proxy.c``). It walks
over all pairs of top-level cells that are close enough to interact (direct
neighbours for hydro, and the cells selected by the opening angle for
gravity). Pairs where both cells are local, or where both cells are foreign,
are skipped. For every remaining pair, ``engine_add_proxy()`` puts the foreign
cell in ``cells_in`` and the local cell in ``cells_out`` of the proxy of the
foreign cell's rank. It creates that proxy if it does not exist yet.

Every rank runs the same walk over the same pairs. A pair made of a cell ``A``
on rank 0 and a cell ``B`` on rank 1 is therefore seen by both ranks, with the
roles swapped:

.. code-block:: none

   rank 0                                   rank 1
   proxy for rank 1                         proxy for rank 0
     cells_out: A   (local, sent)             cells_out: B   (local, sent)
     cells_in:  B'  (foreign copy)            cells_in:  A'  (foreign copy)

What one rank has in ``cells_out`` is what the other rank has in ``cells_in``.
This symmetry allows the send tasks of one rank to match the receive tasks of
the other one (see "Who talks to whom" below).

A proxy also holds the buffers and the MPI requests for two exchanges with its
rank:

* The exchange of the cell trees at each rebuild (``proxy_cells_exchange()``,
  called from ``engine_exchange_cells()``). It sends the description of the
  ``cells_out`` cells and builds the ``cells_in`` cells from what arrives.
* The exchange of stray particles (``proxy_parts_exchange_first()`` and
  ``proxy_parts_exchange_second()``, called from ``engine_exchange_strays()``).
  The ``parts_in``, ``parts_out``, ``gparts_in`` and similar arrays hold the
  particles that left the domain of one rank and now belong to the other one.
  These are transfers of ownership, not foreign copies.

The proxies are not rebuilt at every rebuild of the tree. They are created
again when the grid or the domain decomposition changes: at start-up, after a
regrid (``space_regrid()``) and after a repartition (``engine_repartition()``).
The cell trees are exchanged again at every rebuild.

The lists ``cells_in`` and ``cells_out`` contain top-level cells only. The
sub-cells travel inside the cell tree.

The cell tree of a foreign cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The tree below a top-level cell is not the same everywhere:

* A rank builds the sub-tree of its own top-level cells in ``space_split()``.
* A rank gets the sub-tree of a foreign top-level cell only if that cell is in
  one of its proxies. The owner packs the cell tree in ``cell_pack()`` and the
  other rank rebuilds it in ``cell_unpack()`` (both in ``cell_pack.c``, called
  from ``proxy_cells_exchange()`` in ``proxy.c``).

What a foreign cell contains
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A foreign cell is a copy with a limited content.

* **Bookkeeping copied at rebuild time**, from the ``pcell`` structure
  (``cell.h``, ``cell_unpack()``): particle counts per type, ``h_max`` per
  type, ``ti_end_min``, ``ti_old_part``, ``maxdepth`` and, for gravity, the
  multipole.
* **Particle buffers.** The particles of foreign cells do not live in the
  arrays of the local cells. They are stored in the ``*_foreign`` buffers of
  the space (``parts_foreign``, ``gparts_foreign``, ``sparts_foreign``,
  ``bparts_foreign``), which ``engine_allocate_foreign_particles()`` allocates
  and ``cell_link_foreign_parts()``, ``cell_link_foreign_gparts()`` and
  ``cell_link_sparts()`` point the foreign cells into.
* **No extended particle data.** The linking code sets only ``hydro.parts``.
  No message carries ``xpart`` structures, so a foreign cell has no valid
  ``hydro.xparts``. Do not dereference it.
* **No drift tasks.** The drift, force-finishing and time-integration tasks
  are created for local cells only (``engine_make_hierarchical_tasks_hydro()``
  and its siblings check ``c->nodeID == e->nodeID``). Under
  ``SWIFT_DEBUG_CHECKS``, the drift functions in ``cell_drift.c`` stop with an
  error if they are called on a foreign cell.
* **Sink particles** have no communication channels in this version of the
  code. ``engine_addtasks_send_sinks()`` and ``engine_addtasks_recv_sinks()``
  are placeholders.


Who talks to whom
~~~~~~~~~~~~~~~~~

A common first impression is that a message goes from the cell ``ci`` of a pair
to the cell ``cj`` of the same pair. This is not how it works.

**A message is about one cell.** The rank that owns the cell sends the
particles of that cell. The rank that holds the foreign copy of the *same cell*
receives them into the buffer that the copy points to. The other cell of the
pair takes no part in that message.

A pair of neighbouring cells ``A`` (owned by rank 0) and ``B`` (owned by
rank 1) therefore produces two independent messages per data channel:

.. code-block:: text

      rank 0                                          rank 1

      local A          ===== send task attached to A =====>   foreign copy A'
                                        (recv task attached to A')

      foreign copy B'  <===== recv task attached to B' =====   local B
                                        (send task attached to B)

Each rank sends its own cell and receives the foreign one.

Send and recv tasks are attached to cells
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Every cell has four lists of links to tasks: ``c->mpi.send``, ``c->mpi.recv``,
``c->mpi.pack`` and ``c->mpi.unpack`` (``cell.h``).

* A send task is attached to the *local* cell that it sends.
* A recv task is attached to the *foreign* cell that it fills.

The task itself stores the cell in ``t->ci``. Look at
``engine_addtasks_send_hydro()`` in ``engine_maketasks.c``::

    t_xv = scheduler_addtask(s, task_type_send, task_subtype_xv, ci->mpi.tag,
                             0, ci, cj);

and at ``engine_addtasks_recv_hydro()``::

    t_xv = scheduler_addtask(s, task_type_recv, task_subtype_xv, c->mpi.tag, 0,
                             c, NULL);

In the recv task there is only one cell. In the send task there are two, and
this is the source of the confusion:

* ``ci`` is the cell whose data is sent. The message size and the buffer are
  taken from ``t->ci`` in ``scheduler_enqueue()`` (``scheduler.c``).
* ``cj`` is *not* a cell that receives data. It is a foreign top-level cell of
  the destination rank, and it is used for one thing: ``cj->nodeID`` is the
  rank to send to (it is the ``dest`` argument of ``MPI_Isend`` in
  ``scheduler_enqueue()``). The send-creation functions in
  ``engine_maketasks.c`` document it as a *dummy cell containing the nodeID of
  the receiving node*. ``engine_maketasks()`` fills it with the first foreign
  cell of the proxy (``p->cells_in[0]``). The same ``cj->nodeID`` is used to
  decide whether the send task is needed at all, by checking whether any task
  attached to ``ci`` involves a cell of that rank.

A recv task uses ``t->ci->nodeID`` as the source rank, in ``MPI_Irecv``.

How a send is matched with a recv
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The match is made by the MPI tag and by a communicator that depends on the
task subtype (``subtaskMPI_comms`` in ``task.c``, one duplicate of
``MPI_COMM_WORLD`` per subtype). The tag is stored in ``c->mpi.tag`` and is a
property of the *cell*.

The owning rank assigns the tag in ``cell_ensure_tagged()`` when it creates a
send task. The tags of a whole tree are then sent to the other rank by
``proxy_tags_exchange()``, which is called after the send tasks are created and
before the recv tasks are created (``engine_maketasks()``). The recv tasks then
read ``c->mpi.tag`` from the foreign copy. Both ends of a message use the tag of
the same cell, on the two ranks.

Channels
^^^^^^^^

Each channel is a task subtype. All of them send or receive the data of one
cell, and the type of data is fixed by the subtype (see ``scheduler_enqueue()``
in ``scheduler.c``).

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Subtype
     - Content of the message
   * - ``xv``, ``rho``, ``gradient``, ``part_prep1``, ``rt_gradient``,
       ``rt_transport``
     - The full array of ``struct part`` of the cell, sent at different points
       of the step. The receiver's buffer holds whatever was sent last.
   * - ``limiter``
     - Only the ``time_bin`` of every particle, put in a buffer by a *pack*
       task and read back by an *unpack* task.
   * - ``gpart``, ``fof``
     - The array of ``gpart_foreign`` (or ``gpart_fof_foreign``) of the cell,
       through a pack task.
   * - ``spart_density``, ``spart_prep2``
     - The array of ``struct spart`` of the cell.
   * - ``bpart_rho``, ``bpart_feedback``, ``part_swallow``, ``bpart_merger``
     - The black hole channels.
   * - ``sf_counts``, ``grav_counts``
     - The particle counts (and offsets) of a whole top-level tree, as
       ``pcell_sf_stars`` and ``pcell_sf_grav``. See the section on counts
       below.
   * - ``tend``
     - The time-step information of a whole top-level tree, as
       ``pcell_step``: ``ti_end_min`` and, for each particle type,
       ``dx_max_part``. See the section on ``dx_max_part`` below.

A message whose size depends on the cell is sized by the *receiver's* copy:
``scheduler_enqueue()`` computes the count of a recv from ``t->ci->hydro.count``
(or the equivalent for the other types) on the foreign cell. The receiver must
therefore already know the count that the sender will use. This is why some
counts have their own channel.

Tasks are created once and shared by descendants
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The send functions recurse down the tree of a local cell. At each cell they
look at the tasks attached to the cell (for example ``ci->hydro.density``) and
check whether one of them involves a cell of the destination rank. The first
cell where this is true, going down the tree, creates the send tasks. The
tasks are passed to the recursive calls, and every descendant that also has a
task for that rank gets a *link* to the *same* task in its own ``mpi.send``
list (``engine_addlink()``). The recv functions do the same on the foreign
copy, at the first cell that has a task (for example ``c->hydro.density !=
NULL``).

One send task can therefore serve several cells at different depths, and it
always sends the whole particle array of the cell it was created for
(``t->ci``), including the particles of descendants that have no task of
their own.

When are the tasks active
^^^^^^^^^^^^^^^^^^^^^^^^^

Send and recv tasks are activated in ``cell_unskip.c`` (for example
``cell_unskip_hydro_tasks()``), while it walks over the pair tasks. Only a pair
task with one local and one foreign cell activates them. Pair tasks with two
foreign cells are never created (see ``engine_make_hydroloop_tasks_mapper()``).

For a pair with a foreign ``ci`` and a local ``cj``, the code activates:

* the recv tasks of ``ci->mpi.recv``: the local rank needs the foreign cell
  when its own cell is active;
* the send tasks of ``cj->mpi.send``, toward ``ci->nodeID``: the neighbouring
  rank needs the local cell when the foreign cell is active.

The case ``ci`` local and ``cj`` foreign is the mirror image. Note that the
send is looked up on the *local* cell and the recv on the *foreign* cell. This
is again the same rule: one message, one cell.

A worked example
^^^^^^^^^^^^^^^^

Take two neighbouring top-level cells, ``A`` owned by rank 0 and ``B`` owned by
rank 1, with hydro only and both cells active in the step.

On rank 0, ``A`` is local and ``B`` is a foreign copy (called ``B'``). The
task graph of rank 0 contains:

* the pair density and force tasks of ``(A, B')``, run by rank 0;
* ``send_xv(A)`` and ``send_rho(A)`` attached to ``A``, with destination
  rank 1;
* ``recv_xv(B')`` and ``recv_rho(B')`` attached to ``B'``, with source rank 1.

On rank 1 it is the mirror image: the pair tasks of ``(A', B)``, ``send_xv(B)``
and ``send_rho(B)``, ``recv_xv(A')`` and ``recv_rho(A')``.

The messages of the ``xv`` channel are:

.. code-block:: text

      rank 0                                   rank 1
      A.parts   ======== tag(A) ========>      A'.parts   (in parts_foreign)
      B'.parts  <======= tag(B) =========      B.parts
      (in parts_foreign)

and the same two messages again for the ``rho`` channel, later in the step.
Rank 0 reads ``B'`` in its copy of the pair tasks only after ``recv_xv(B')``
has completed, because ``recv_xv`` unlocks the pair tasks and the sort of ``B'``
(``engine_addtasks_recv_hydro()``). The ``rho`` receive is unlocked by the
pair density tasks, so that the buffer of ``B'`` is not overwritten while a
task is reading it.


Both ranks run the pair task
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The task graph is built by every rank over the full top-level grid, and only
pairs with two foreign cells are skipped. A pair of cells owned by two
different ranks therefore exists as **two separate task instances**, one on
each rank, and each rank runs its own.

Each instance has one local side and one foreign side. From the point of view
of one instance:

* The **local** cell is up to date. The task can read and write it, subject to
  the usual locks and dependencies.
* The **foreign** cell is a read-only snapshot. Its particle data are those of
  the last message that arrived for the channels that were active. Its
  bookkeeping fields (see the next section) may be older than the local ones.

The two instances do not share memory. Each receive overwrites the buffer of
the foreign cell with the data that the owner sent.


Foreign data that lags by design
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Several fields of a foreign cell are refreshed by a message that is sent at a
different time than the corresponding fields of the local cell change. Code
that runs on a foreign cell must not assume that these fields are current.

The dx_max_part of a foreign cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``dx_max_part`` is the largest displacement of any particle in the cell since
the last tree rebuild. The rebuild criterion
(``cell_need_rebuild_for_hydro_pair()`` in ``cell.h``) and several debug checks
use it.

On the **local** cell, ``dx_max_part`` is updated by the drift
(``cell_drift_part()`` and its siblings in ``cell_drift.c``), so it is always
current.

On the **foreign** copy it is set at two moments only:

1. At a rebuild. ``space_rebuild_recycle_mapper()`` resets it to zero on every
   top-level cell, and ``cell_unpack()`` creates the sub-cells of a foreign
   tree with a value of zero.
2. At the end of a step, by the ``tend`` channel. The owner packs the value in
   ``cell_pack_end_step()`` when the send task is enqueued in
   ``scheduler_enqueue()``, and the other rank stores it in
   ``cell_unpack_end_step()`` (``runner_main.c``). Only the top-level cells that
   the time-step tasks marked as updated are exchanged
   (``space_mark_cell_as_updated()``, ``engine_synchronize_times()``).

The receive of the particle data does not touch it: ``runner_do_recv_part()``
updates ``h_max``, ``h_max_active`` and ``ti_old_part`` of the foreign cell, but
not ``dx_max_part``.

The order of operations in ``engine_step()`` (``engine.c``) is:

.. code-block:: text

      engine_prepare()
          engine_unskip()             activate tasks, evaluate rebuild criteria
          (allreduce of forcerebuild over all ranks)
          engine_rebuild()            only if a rebuild is due; resets dx_max_part
                                      and unskips again once the tasks are new
      engine_launch("tasks")          drifts, pair tasks, recv/send, timesteps
      engine_synchronize_times()      tend exchange: foreign dx_max_part refreshed
      engine_collect_end_of_step()    decides forcerebuild for the NEXT step

A rebuild always happens at the start of a step, before the tasks of that step
are launched.

This gives the following timeline for a rebuild at step ``N``:

.. code-block:: text

      step N     rebuild: dx_max_part = 0 on local and foreign copies.
                 Particles were drifted to the current time before the
                 rebuild, so the drift tasks of this step have nothing left
                 to move. The tend at the end of N sends 0.
      step N+1   the drift tasks move the particles: the local dx_max_part
                 grows during the step. The foreign copy still holds the value
                 of the tend of step N, which is 0.
      step N+2   the foreign copy holds the local value of the end of step
                 N+1. The local value has grown again during N+2.

The lag is always there: during any step, a foreign ``dx_max_part`` is the
value at the end of the previous step. It is most visible on the first step
after a rebuild, because the foreign value is exactly zero while the local
value already has a non-zero size.

**Consequence for checks.** The ``SWIFT_DEBUG_CHECKS`` blocks in
``runner_doiact_functions_hydro.h``, ``runner_doiact_functions_limiter.h``,
``runner_doiact_functions_stars.h`` and ``cache.h`` verify that particles are in
the right frame. They compare the position of a particle with respect to the
cell corner against a threshold built from the cell width and ``dx_max_part``
(``shift_threshold_x``, ``shift_threshold_y``, ``shift_threshold_z``). A
threshold that uses the ``dx_max_part`` of a foreign cell is too small for a
particle that has legitimately moved. Such a check must use a bound for the
foreign side that does not rely on the lagged value. A cell width is a
possible choice: the rebuild criterion above keeps ``dx_max_part`` below the
cell size (``dmin``) at the start of each step, so a cell width is a generous
bound. The check on the foreign side is then loose, but the rank that owns the
cell runs the same pair task with the current value and applies the tight
check there.

The same fields are read by ``cell_need_rebuild_for_hydro_pair()``, so the
rebuild criterion sees the lagged value on the foreign side too.

The time-bins of a foreign cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When the time-step limiter is on, the owner sends the ``time_bin`` of every
particle in the ``limiter`` channel. The receiver stores them in
``cell_unpack_timebin()`` (``cell_pack.c``), called by
``runner_do_unpack_limiter()`` (``runner_pack.c``).

The pack task is not ordered with the tasks that write the bins:

* ``engine_addtasks_send_hydro()`` makes the limiter pack task depend on
  ``ci->super->timestep`` only. There is no dependency on
  ``timestep_limiter`` or ``timestep_sync``.
* ``task_lock()`` (``task.c``) has no case for ``task_type_pack``, so the pack
  task takes no cell lock. The tasks ``timestep_limiter`` and
  ``timestep_sync`` lock the cell.

``timestep_sync_part()`` (``timestep_sync.h``) writes a *temporary negative*
bin, ``p->time_bin = -min_active_bin``, and ``runner_do_sync()``
(``runner_time_integration.c``) replaces it by the new bin a few statements
later. A pack task that runs at that moment copies the negative value. The
receiver can therefore see a ``time_bin`` that is zero or negative for a
particle that is being re-timestepped in this step.

Such a particle counts as active: the rule used everywhere is
``time_bin <= max_active_bin``. Code that handles unpacked bins must apply the
same comparison and must not call functions that check the bin against the
current time:

* ``part_is_active()`` (``active.h``) computes, under ``SWIFT_DEBUG_CHECKS``,
  the end of the step from the bin before it compares. ``get_integer_timestep()``
  returns zero for a bin that is zero or negative, so the end of the step is
  zero, and ``part_is_active()`` stops with the error *particle in an impossible
  time-zone*. ``part_is_starting()`` has a similar debug check that does not
  fire for such a bin, but it also derives a time from the bin: do not rely on
  it for unpacked data.
* ``part_is_active_no_debug()`` and a direct comparison with
  ``e->max_active_bin`` are safe. ``runner_do_recv_part()`` (``runner_recv.c``)
  uses the direct comparison.

Counts and the layout of the foreign buffers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Particle creation (star formation, sink star formation) changes the number of
particles of a cell between two rebuilds, and the receiver sizes its receive
from its own count (see the section on channels). The ``sf_counts`` and
``grav_counts`` channels carry the new counts of a whole top-level tree
(``cell_pack_sf_counts()``, ``cell_pack_grav_counts()``), and the receiver
stores them in ``cell_unpack_sf_counts()`` and ``cell_unpack_grav_counts()``.
These tasks exist at the top level only.

The order between a count message and the data message is set in
``engine_maketasks.c``, and it is not the same for all types:

* For stars, the ``sf_counts`` receive unlocks the ``spart_density`` receive
  (``engine_addtasks_recv_stars()``): the counts arrive first, and the sender
  makes the send of ``spart_density`` wait for ``sf_counts`` in the same way.
* For gravity, when star formation is on, the ``gpart`` send unlocks star
  formation, which unlocks the
  ``grav_counts`` send (``engine_addtasks_send_gravity()``), and the ``gpart``
  receive unlocks the ``grav_counts`` receive
  (``engine_addtasks_recv_gravity()``). The particle data therefore describe
  the state before star formation, and the counts describe the state after it.

When you add a channel whose message size can change during a step, decide
explicitly which count the receiver uses for each message, and make the task
dependencies enforce it.

The receiving rank does not reproduce the layout of the sender's arrays.
``cell_link_foreign_parts()`` and ``cell_link_foreign_gparts()`` walk down the
foreign tree and attach the buffer at the level where the recv task exists
(``cell_get_recv()``). A part of the tree that has no recv task takes no space
in the foreign buffer, whereas the sender's array of a top-level cell contains
all the particles of that cell. An index or an offset that is valid in the
sender's array is therefore not, in general, valid in the receiver's buffer.
Derive the position of a foreign cell in the buffer from the receiver's own
tree.

Drifts and sorts of foreign cells
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* **Foreign cells are never drifted.** Only the owner drifts. The send tasks
  depend on the drift task of the super cell (``engine_addtasks_send_hydro()``,
  ``engine_addtasks_send_stars()``, ``engine_addtasks_send_gravity()``), so the
  data leave the rank already drifted. After a receive of ``xv``,
  ``runner_do_recv_part()`` sets ``c->hydro.ti_old_part = ti_current``, which
  is what makes ``cell_are_part_drifted()`` true for the foreign cell.
* **Foreign cells are sorted by the receiver.** The sort task of a foreign super
  cell exists, and ``recv_xv`` unlocks it. ``runner_do_recv_part()`` clears the
  ``sorted`` flags when ``xv`` arrives (``clear_sorts``), because the
  positions have changed. The sort arrays of a foreign cell are therefore
  allocated, resized and freed by tasks of the receiving rank, with the usual
  locks (``hydro.extra_sort_lock``).
* **The drift must target the cell of the send task.** A send task is shared
  by descendants at several depths, and it sends the whole particle array of
  ``t->ci``. When a pair task activates a send, the drift of the sent data
  must be activated for ``t->ci``, not for the cell of the pair task that
  happened to trigger the activation. ``scheduler_activate_send()`` returns the
  link to the task for this reason. ``cell_unskip_stars_tasks()`` shows the
  pattern::

      struct link *l_send_spart = scheduler_activate_send(
          s, cj->mpi.send, task_subtype_spart_density, ci_nodeID);
      ...
      cell_activate_drift_spart(l_send_spart->t->ci, s);

  If only the triggering cell is drifted, particles of the other descendants
  that share the same send task are sent undrifted.


Writing MPI-safe code in SWIFT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Never assume that a field of a foreign cell is current. Particle data are as
  fresh as the last message of the channel. ``dx_max_part`` is at best one step
  old. Time-bins can hold a transient value.
* Do not call helpers that check the state against the current time on data that
  just arrived. Compare the bins directly (``time_bin <= max_active_bin``).
* Before dereferencing an array of a foreign cell, check that it exists.
  ``hydro.xparts`` does not.
* A debug check that bounds motion with ``dx_max_part`` needs a foreign-side
  bound that does not use the lagged value.
* Create a send or recv task once, at the highest level that needs it, and link
  it to the descendants. Decide whether a task exists from the existence of the
  shared task, not from a local particle count of the cell that you happen to
  be visiting.
* When a send task is shared, activate drifts (and anything else that must
  precede the send) for ``t->ci`` of the task, using the link that
  ``scheduler_activate_send()`` returns.
* Give every recv the same count that the sender uses, and enforce the order of
  count and data messages with task dependencies.
* Remember that the pair task of a cross-rank pair runs on both ranks, with the
  local and foreign roles swapped. A result that depends on reading the foreign
  cell must be valid on both.
