#!/bin/bash
t=$1
read P T0 < <(.venv/bin/python -c "import json;d=json.load(open('res/$t.json'));print(d['P'],d['T0'])")
.venv/bin/python vet.py $t >/dev/null 2>&1
D=$(.venv/bin/python -c "import json;print(json.load(open('vet_$t.json'))['dur_h'])")
.venv/bin/python xcheck.py $t $P $T0 $D 2>&1 | grep -v Warn > xcheck_$t.txt
.venv/bin/python shape.py $t $P $T0 > shape_$t.txt 2>&1
.venv/bin/python shape.py $t $P $T0 QLP > shape_${t}_qlp.txt 2>&1
