"""Portable launcher. Optional project-local .deps is used only when installed."""
from pathlib import Path
import runpy
import sys

root=Path(__file__).resolve().parent
sys.path.insert(0,str(root))
if (root/'.deps').is_dir():sys.path.insert(0,str(root/'.deps'))
commands={'science':'experiments.science','reproduce':'experiments.reproduce_version1',
          'benchmark':'experiments.benchmark','manuscript':'experiments.build_manuscript',
          'fetch':'experiments.fetch_population_data'}
if len(sys.argv)<2 or sys.argv[1] not in {*commands,'test'}:
    raise SystemExit('Usage: python run.py {test|reproduce|science|benchmark|manuscript|fetch} [arguments]')
command=sys.argv.pop(1)
if command=='test':
    import pytest
    raise SystemExit(pytest.main([str(root/'tests'),*sys.argv[1:]]))
runpy.run_module(commands[command],run_name='__main__')
