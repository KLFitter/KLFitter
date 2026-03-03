#!/bin/bash
# Copyright (c) 2009--2026, the KLFitter developer team
#
# This file is part of KLFitter.
#
# KLFitter is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# KLFitter is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with KLFitter. If not, see <http://www.gnu.org/licenses/>.
#
# Helper script to apply patches to the BAT source directory.
# Usage: patch_bat.sh [source directory] [patch directory]

set -e

SOURCE_DIR=$1
PATCH_DIR=$2

if [ ! -d "$SOURCE_DIR" ]; then
    echo "Error: Source directory $SOURCE_DIR does not exist."
    exit 1
fi

if [ ! -d "$PATCH_DIR" ]; then
    echo "Error: Patch directory $PATCH_DIR does not exist."
    exit 1
fi

cd "$SOURCE_DIR"

for p in "$PATCH_DIR"/*.patch; do
    if [ -f "$p" ]; then
        echo "Applying patch: $p"
        patch -p1 < "$p"
    fi
done
