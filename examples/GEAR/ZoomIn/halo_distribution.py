#!/usr/bin/env

import yt
from yt.extensions.astro_analysis.halo_analysis.api import HaloCatalog


def doPlot(f, name, i):
    hc = HaloCatalog(data_ds=f, finder_method="fof",
                     finder_kwargs={"padding": 0.02})
    hc.create()
    print(hc)
    print(dir(hc))
    exit()
