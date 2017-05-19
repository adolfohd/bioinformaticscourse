#!/bin/bash
rm log/out.*
rm log/log.*
rm log/err.*
echo "log files deleted"
echo "submitting to condor..."
condor_submit submit
echo "submission to condor complete"
