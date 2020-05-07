Various hints on settings needed to get various MPIs running with SWIFT.

Last update 5th May 2020.


## Intel MPI

_Intel MPI 2018_ usually runs without any needs for special settings.

_Intel MPI 2019 and 2020_ can run for small tests, but without flags will
generally deadlock in the MPI exchanges of the engine, or worse. In that case
try the following settings.

```
  FI_OFI_RXM_RX_SIZE=4096
  FI_OFI_RXM_TX_SIZE=4096
  FI_UNIVERSE_SIZE=2048
```

If you want use the `release_mt` library, then you also need to use:

```
  source $I_MPI_ROOT/intel64/bin/mpivars.sh release_mt
```

when initializing the library environment. Some success has also been seen
using the asynchronous progression settings:

```
  I_MPI_ASYNC_PROGRESS=1
  I_MPI_ASYNC_PROGRESS_THREADS=1

```
(note these are tested with `2019 update-4` and `2020 update-1`).

## OpenMPI

_Open MPI_ comes in many flavours with many combinations of underlying
transport libraries and running on many different fabrics. A complete
description of all combinations is beyond the scope of this guide.

On Mellanox hardware, we have had success running version 4.0 with the
UCX layer version 1.6 and using the following settings:

```
-mca coll_hcoll_enable 0
UCX_TLS=ud_x,shm,self
UCX_RC_MLX5_TM_ENABLE=n
UCX_DC_MLX5_TM_ENABLE=n
```
