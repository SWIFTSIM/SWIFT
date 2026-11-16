# This file is part of SWIFT.
# Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
#                    Matthieu Schaller (schaller@strw.leidenuniv.nl).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Non-standard files that should be part of the distribution.
# Every file tracked by git under examples/ is copied into the tarball,
# so new examples are included without having to be listed here.
# Outside a git checkout (e.g. 'make distcheck' re-running 'make dist' from
# the unpacked tarball), all the files present are copied instead, except
# build products.
# Files automake has already copied (Makefile.am, Makefile.in) are skipped.
# Symbolic links are copied as links (-P); they all point inside examples/.
dist-hook:
	@dest=`cd "$(distdir)" && pwd`; \
	cd "$(srcdir)" && \
	if git ls-files --error-unmatch Makefile.am >/dev/null 2>&1; then \
	  git ls-files .; \
	else \
	  find . \( -name .deps -o -name .libs \) -prune -o \
	    \( -type f -o -type l \) ! -name Makefile ! -name '*.o' \
	    ! -name '*.lo' ! -name '*.la' ! -name cooling_rates -print; \
	fi | while read -r f; do \
	  test -e "$$dest/$$f" || test -h "$$dest/$$f" || \
	  { mkdir -p "$$dest/`dirname "$$f"`" && cp -pP "$$f" "$$dest/$$f"; } || exit 1; \
	done
