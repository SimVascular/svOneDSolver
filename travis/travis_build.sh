#!/bin/bash

set -e

if $WITH_CMAKE; then
  MAKE="make --jobs=$NUM_THREADS --keep-going"
  mkdir -p $BUILD_DIR
  cd $BUILD_DIR
  source $SCRIPTS/travis_cmake_config.sh
  pushd $BUILD_DIR
  $MAKE
  popd
fi

