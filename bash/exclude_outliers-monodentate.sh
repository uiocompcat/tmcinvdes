#!/bin/bash
python -m tmcinvdes.analysis.exclude_outliers -d monodentate \
                                              -i datasets/07_uncond-labeled \
                                              -o datasets/08a_uncond-excluded \
                                              -r datasets/08b_uncond-included \
                                              -x test
