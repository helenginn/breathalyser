#!/bin/bash

check_version=`grep 'CHECK_VERSION_NUMBER' $1 | cut -d' ' -f3`
check_commit="notcompiledfromgit"

command -v git > /dev/null 2>&1
if [ $? -eq 0 ]; then
	check_commit="."`git rev-parse HEAD`
fi

sed "s/_version_commit_/${check_version}${check_commit}/g" $1 > $2


