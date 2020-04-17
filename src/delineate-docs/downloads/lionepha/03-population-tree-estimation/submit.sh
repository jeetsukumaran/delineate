#! /bin/bash
set -e -o pipefail
set -x
qsub 001-rep001-hkyg_strictclock_bd_poplineages_p095.job >> submitted.log
sleep 1.0
qsub 001-rep002-hkyg_strictclock_bd_poplineages_p095.job >> submitted.log
sleep 1.0
qsub 001-rep003-hkyg_strictclock_bd_poplineages_p095.job >> submitted.log
sleep 1.0
qsub 001-rep004-hkyg_strictclock_bd_poplineages_p095.job >> submitted.log
sleep 1.0
qsub 001-rep005-hkyg_strictclock_bd_poplineages_p095.job >> submitted.log
sleep 1.0
qsub 001-rep006-hkyg_strictclock_bd_poplineages_p095.job >> submitted.log
sleep 1.0
