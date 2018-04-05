#!/bin/bash
nfails=0
ntests=0
fn="$1"
r="$2"
d="$(dirname $0)"
ffn=.full_joint.txt
# generate test output from joint script
echo "generating test output. overwriting ${ffn} if it exists..."
python delineate/joint_from_all_scenarios.py "${fn}" "${r}" > "${ffn}"

sjfn=.sorted_joint_from_joint.txt
echo "generating sorted test output. overwriting ${sjfn} if it exists..."
cat "${ffn}" | grep '^Pr[(] ' | sort > ${sjfn}

osjfn=".eff-${sjfn}"
echo "testing efficient_joint_prob.py . overwriting ${osjfn} if it exists..."
python delineate/efficient_joint_prob.py "${fn}" "${r}" | sort > "${osjfn}" || exit 1
linenum=1
nl=$(wc -l "${sjfn}" | awk '{print $1}')
nl=$(expr 1 + $nl)
while test $linenum -lt $nl ; do
    tline="$(head -n$linenum $sjfn | tail -n1)" || exit
    oline="$(head -n$linenum $osjfn | tail -n1)" || exit
    ttaxa=$(echo $tline | sed 's/.*[(]//' | sed 's/[)].*//')
    otaxa=$(echo $oline | sed 's/.*[(]//' | sed 's/[)].*//')
    if ! test $ttaxa = $otaxa ; then
        echo "line $linenum of test output is for \"$ttaxa\" but in efficient_joint_prob output it was for \"$otaxa\""
        exit 1
    fi
    tprob=$(echo $tline | sed 's/.*[=] //')
    oprob=$(echo $oline | sed 's/.*[=] //')
    linenum=`expr 1 + $linenum`
    if python -c "assert(abs(${oprob} - ${tprob})/${tprob} < 1e-08)" 2>/dev/null ; then
        echo "  correct         $ttaxa     $oprob   $tprob"
        ntests=$(expr 1 + $ntests)
    else
        echo "INCORRECT         $ttaxa     $oprob   $tprob"
        nfails=$(expr 1 + $nfails)
    fi
done


tfn=.marg_from_joint.txt
echo "generating marginal test output. overwriting ${tfn} if it exists..."
cat "${ffn}" | sed '/^[^P]/d' | sed '/[}][{]/d' | sed '/[(] [{]/d' > "${tfn}"

# compare each line to the calc from the marginal script...
linenum=1
nl=$(wc -l "${tfn}" | awk '{print $1}')
nl=$(expr 1 + $nl)
while test $linenum -lt $nl ; do
    line="$(head -n$linenum $tfn | tail -n1)" || exit
    taxa=$(echo $line | sed 's/.*[{]//' | sed 's/[}].*//' | sed 's/,/ /g')
    jprob=$(echo $line | sed 's/.*[=] //')
    linenum=`expr 1 + $linenum`
    margout=$(python delineate/marginal_prob.py "${fn}" "${r}" $taxa) || exit
    mprob=$(echo $margout | sed 's/.*[=] //')
    if python -c "assert(abs(${mprob} - ${jprob})/${jprob} < 1e-08)" 2>/dev/null ; then
        echo "  correct         $taxa     $mprob   $jprob"
        ntests=$(expr 1 + $ntests)
    else
        echo "INCORRECT         $taxa     $mprob   $jprob"
        nfails=$(expr 1 + $nfails)
    fi
done

echo "${ntests} tests ${nfails} errors"
exit $nfails
