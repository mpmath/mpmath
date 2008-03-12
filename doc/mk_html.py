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
    padding:5px;
    border-left: 25px solid #ccc;
    background-color: #ccc;
    font-size: 13px; font-family: arial, sans-serif;
    line-height:1.5em;
}

div.document {
    max-width: 700px;
    color: black;
    background-color: white;
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

""")

open("manual.html", "w").write("".join(lines))