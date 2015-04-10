#!/bin/bash

arr=($(seq 0 89))

for i in ${arr[@]}; do
	mkdir "sim${i}"
done

echo sim* | xargs -n 1 cp stat.m emle_option_con.m  emle_option_Lagranian.m  emle_optionscore.m
