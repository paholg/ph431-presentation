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
mayavi_cmds = []

slice_args = ['B y', 'B z', 'E x', 'S y', 'S z']
mayavi_args = ['3d E', '3d B', '3d S']

latex_docs = ['presentation.tex']


for arg in slice_args:
  py_cmds.append('plot-vortices.py slice ' + arg)

for arg in mayavi_args:
  mayavi_cmds.append('plot-vortices.py ' + arg)

print py_cmds
print mayavi_cmds

for f in py_cmds:
  if os.path.getmtime(f.split(' ')[0]) > build_time or build_all:
    print 'Running: python', f
    subprocess.call(['python'] + f.split(' '))

xcmd = 'xvfb-run -n 99 --server-args="-screen 0 1024x768x24" ./'

for f in mayavi_cmds:
  if os.path.getmtime(f.split(' ')[0]) > build_time or build_all:
    print 'Running:', xcmd + f
    subprocess.call(xcmd + f, shell=True)



for f in latex_docs:
  if os.path.getmtime(f.split(' ')[0]) > build_time or build_all:
    print 'Running: latex', f
    subprocess.call(['pdflatex'] + f.split(' '))


print 'all done, updating time'
with file(last_built, 'a'):
  os.utime(last_built, None)
