#!/bin/bash
cd ../../ParMooN_CMG/BUILD
if `make -j8` | grep -iF 'error'; then
  echo "error"
  exit
fi
