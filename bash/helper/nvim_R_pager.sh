#!/bin/bash
if command -v nvim &> /dev/null; then
  col -b | nvim +"set ft=help" -
  exit $?  # Exit with nvim's status code
elif command -v vim &> /dev/null; then
  col -b | vim +"set ft=help" -
  exit $?  # Exit with vim's status code
else
  echo "Error: Neither nvim nor vim found in $PATH" >&2
  exit 1
fi
