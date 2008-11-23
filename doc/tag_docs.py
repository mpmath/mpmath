script = \
"""
<!-- Generate pageview statistics when this document is viewed on the mpmath website -->
<script src="http://www.google-analytics.com/urchin.js" type="text/javascript"></script>
<script type="text/javascript">

if ((""+document.location).match("google"))
{
    _uacct = "UA-2697185-2";
    urchinTracker();
}
</script>
"""

tag = "Generate pageview statistics"

import os
import os.path

paths = ["build", "build/functions", "build/calculus"]

for path in paths:
    for fname in os.listdir(path):
        if fname.endswith(".html"):
            f = open(os.path.join(path, fname), "r+w")
            if script not in f.read():
                f.seek(0)
                lines = f.readlines()
                lines.insert(lines.index("  </body>\n"), script)
                f.seek(0)
                f.write("".join(lines))
                print "modified", fname
