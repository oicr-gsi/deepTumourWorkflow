#!/bin/bash
set -o nounset
set -o pipefail

cd $1


ls | sort

find -name *.json -xtype f | md5sum" \; 

