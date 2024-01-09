#!/bin/bash

echo -n "last test run:" > ../tests/test_log.txt
date >> ../tests/test_log.txt
py.test -q -s ../tests/test_*.py > /tmp/test_output
cat /tmp/test_output | tail -1 >> ../tests/test_log.txt
cat /tmp/test_output >> ../tests/test_log.txt
rm /tmp/test_output