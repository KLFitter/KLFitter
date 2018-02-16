#!/bin/sh
# Copyright (c) 2009--2017, the KLFitter developer team
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
# ==========================================================
#
# This is a short script to download and build the BAT library. The
# individual steps include:
#  - Download a tar archive of BAT v0.9.4.1.
#  - Verify the integrity of the tar file.
#  - Extract the source.
#  - Configure and build the source.
#  - Copy the source into the installation path ($1).
#  - Clean up.
#
# Usage: CompileBAT.sh [absolute install path]

set -e

if [ $# -eq 1 ]
then
    echo "Downloading and installing BAT into $1."
else
    echo "Usage: $0 /path/to/install/BAT"
    exit 1
fi

# Define the directories.
base_dir=$PWD
build_dir=$base_dir/BATBuild
target_dir=$1
mkdir -p $build_dir $target_dir

# This function verifies the SHA256 ID of a given file.
# Call with: verify_file_hash [file name] [valid hash]
verify_file_hash() {
    # Extract the SHA256 hash from the file.
    file_hash=`sha256sum "$1"`
    file_hash=${file_hash%$1}
    file_hash="$(echo $file_hash | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"

    if [ "$file_hash" != "$2" ]; then
        echo "SHA256 hash of $1 is invalid. Aborting"
        return 1
    else
        echo "SHA256 hash of $1 matches."
        return 0
    fi
}

# Download the BAT tar archive and verify it.
cd $build_dir
tar_name="BAT-0.9.4.1.tar.gz"
valid_hash="d46c6f834cb5888bbf4db393887190380132fa48816e0804f79c4a3cc344ef87"
wget https://github.com/bat/bat/releases/download/v0.9.4.1/$tar_name
if !verify_file_hash "$tar_name" "$valid_hash"; then
    exit 1;
fi

# Perform the actual configure and make commands.
tar xzf "$tar_name"
cd BAT-0.9.4.1
./configure --with-rootsys=`root-config --prefix` --prefix=$target_dir
make -j || make -j || make -j

# Copy built BAT into the installation path ($1).
make install

# Go back to the start and clean up.
cd $base_dir
rm -rf $build_dir
