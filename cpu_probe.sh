#!/bin/bash

echo "Probing /proc/cpuinfo to update cpu.mk"
echo -n "Have:"

if grep sse2 /proc/cpuinfo > /dev/null
then
    echo -n " SSE2"
    echo "SSE2= 1" > cpu.mk
fi

if grep avx /proc/cpuinfo > /dev/null
then
    echo -n " AVX"
    echo "AVX= 1" > cpu.mk
fi

if grep avx2 /proc/cpuinfo > /dev/null
then
    echo -n " AVX2"
    echo "AVX2= 1" > cpu.mk
fi

echo
