#!/bin/bash
fn="$1"
r="$2"
d="$(dirname $0)"
tfn=.from_joint.txt
# generate test output from joint script
bash "${d}/marginal_from_joint.sh" "${fn}" "${r}" > "${tfn}"
# compare each line to the calc from the marginal script...
linenum=1
nl=$(wc -l "${tfn}" | awk '{print $1}')
while test $linenum -lt $nl ; do
    line="$(head -n$linenum $tfn | tail -n1)" || exit
    taxa=$(echo $line | sed 's/.*[{]//' | sed 's/[}].*//' | sed 's/,/ /g')
    jprob=$(echo $line | sed 's/.*[=] //')
    linenum=`expr 1 + $linenum`
    margout=$(python delineate/marginal_prob.py "${fn}" "${r}" $taxa) || exit
    mprob=$(echo $margout | sed 's/.*[=] //')
    if python -c "assert(abs(${mprob} - ${jprob})/${jprob} < 1e-08)" 2>/dev/null ; then
        echo "  correct         $taxa     $mprob   $jprob"
    else
        echo "INCORRECT         $taxa     $mprob   $jprob"
    fi
done
