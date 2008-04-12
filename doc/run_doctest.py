import os
import os.path

path = "source"

import doctest
for f in os.listdir(path):
    if f.endswith(".txt"):
        doctest.testfile(os.path.join(path, f), module_relative=False)
