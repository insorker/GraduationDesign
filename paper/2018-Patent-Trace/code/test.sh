#!/bin/zsh

for i in {1..10}; do
  ./a.out | tail -n1
  sleep 1s
done
