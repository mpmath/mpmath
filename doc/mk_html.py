from docutils.core import publish_file

publish_file(
  source_path='manual.rst',
  destination_path='manual.html',
  writer_name='html')

lines = open("manual.html", "r").readlines()

lines.insert(lines.index("</style>\n")-1, \
"""
body {
    margin:0px;
    padding:25px;
    background-color: #ccc;
    font-size: 13px; font-family: arial, sans-serif;
    line-height:1.5em;
}

div.document {
    max-width: 700px;
    color: #000;
    background-color: #fff;
    padding:25px;
    border:5px solid #ddd;
}

h1 {
    margin-top:0;
    padding-top:20px;
}

table {
    border-collapse: collapse;
}

pre, tt {
  font-family: consolas, lucida console, courier new, monospace;
}

pre.literal-block, pre.doctest-block {
  line-height:1.3em;
  border-top:1px solid #ccc;
  border-bottom:1px solid #ccc;
  background-color:#f0f0f0;
}

""")

lines.insert(lines.index("</body>\n"), \
"""
<!-- Generate pageview statistics -->
<script src="http://www.google-analytics.com/urchin.js" type="text/javascript"></script>
<script type="text/javascript">
_uacct = "UA-2697185-2";
urchinTracker();
</script>
""")

open("manual.html", "w").write("".join(lines))