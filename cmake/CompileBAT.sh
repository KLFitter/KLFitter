#!/bin/sh

set -e

# Define the directories.
base_dir=$PWD
build_dir=$base_dir/BATBuild
target_dir=$base_dir/external/BAT

# Make the directories..
mkdir -p $build_dir
mkdir -p $target_dir

# This function verifies the SHA256 ID of a given file.
verify_file_hash() {
    # Extract the SHA256 hash from the file.
    file_hash=`sha256sum "$1"`
    file_hash=${file_hash%$1}
    file_hash="$(echo $file_hash | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"

    if [ "$file_hash" != "$2" ]; then
        echo "SHA256 hash of $1 is invalid. Aborting"
        return 0
    else
        echo "SHA256 hash of $1 matches."
        return 1
    fi
}

# Download the archive and verify it.
cd $build_dir
tar_name="BAT-0.9.4.1.tar.gz"
valid_hash="d46c6f834cb5888bbf4db393887190380132fa48816e0804f79c4a3cc344ef87"
wget https://github.com/bat/bat/releases/download/v0.9.4.1/$tar_name
if verify_file_hash "$tar_name" "$valid_hash"; then
    exit 1;
fi

# Perform the actual configure and make commands.
tar xzf "$tar_name"
cd BAT-0.9.4.1
./configure --with-rootsys=`root-config --prefix` --prefix=$target_dir
make -j || make -j || make -j
make install
cd $base_dir
rm -rf $build_dir
