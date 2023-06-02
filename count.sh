#!/bin/zsh

echo -n "lines: "
find . -not -path "./node_modules" -type f \( -name "*.c" -o -name "*.h" \) -exec wc -l {} + | awk '{total += $1} END{print total}'
