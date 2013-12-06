#!/usr/bin/python
import os.path, time, subprocess

last_built = '.updated'

if os.path.exists(last_built):
  build_time = os.path.getmtime(last_built)
else:
  build_time = 0

if os.path.getmtime('make.py') > build_time:
  build_all = True
  print 'make file changed, rebuilding all'
else:
  build_all = False


py_cmds = []
slice_cmds = ['B y', 'B z', 'E x', 'S y', 'S z']
for cmd in slice_cmds:
  py_cmds.append('plot-vortices.py slice ' + cmd)
print py_cmds

latex_docs = ['presentation.tex']

for f in py_cmds:
  if os.path.getmtime(f.split(' ')[0]) > build_time or build_all:
    print 'Running: python', f
    subprocess.call(['python'] + f.split(' '))
for f in latex_docs:
  if os.path.getmtime(f.split(' ')[0]) > build_time or build_all:
    print 'Running: latex', f
    subprocess.call(['pdflatex'] + f.split(' '))


print 'all done, updating time'
with file(last_built, 'a'):
  os.utime(last_built, None)
