#!/bin/bash
fn="$1"
r="$2"
python delineate/joint_from_all_scenarios.py "${fn}" "${r}" | sed '/^[^P]/d' | sed '/[}][{]/d' | sed '/[(] [{]/d'
