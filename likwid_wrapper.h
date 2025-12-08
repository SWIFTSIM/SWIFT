#ifndef SWIFT_LIKWID_WRAPPER_H
#define SWIFT_LIKWID_WRAPPER_H

#include <config.h>
#include "inline.h"

#ifdef WITH_LIKWID
#define LIKWID_PERFMON
#include <likwid-marker.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#endif

__attribute__((always_inline)) INLINE void swift_likwid_marker_init(void) {
	    LIKWID_MARKER_INIT;
}

static __attribute__((always_inline)) INLINE void swift_likwid_marker_register(const char* regionTag) {
	    LIKWID_MARKER_REGISTER(regionTag);
}

static __attribute__((always_inline)) INLINE void swift_likwid_marker_start_region(const char* regionTag) {
	    LIKWID_MARKER_START(regionTag);
}

static __attribute__((always_inline)) INLINE void swift_likwid_marker_stop_region(const char* regionTag) {
	    LIKWID_MARKER_STOP(regionTag);
}

static __attribute__((always_inline)) INLINE void swift_likwid_marker_close(void) {
	    LIKWID_MARKER_CLOSE;
}

#endif /* SWIFT_LIKWID_WRAPPER_H */

