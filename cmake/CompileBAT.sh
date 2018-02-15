#!/bin/sh

set -e

# Save the build directory.
build_dir=$PWD

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
tar_name="BAT-0.9.4.1.tar.gz"
valid_hash="d46c6f834cb5888bbf4db393887190380132fa48816e0804f79c4a3cc344ef87"
wget https://github.com/bat/bat/releases/download/v0.9.4.1/$tar_name
if verify_file_hash "$tar_name" "$valid_hash"; then
    exit 1;
fi

# Extract the tar archive, move its contents into "BATBuild" folder.
tar xzf BAT-0.9.4.1.tar.gz
mv BAT-0.9.4.1 BATBuild

# Perform the actual configure and make commands.
cd BATBuild
./configure --with-rootsys=$ROOTSYS --prefix=$build_dir
make -j || make -j || make -j
make install
cd $build_dir
