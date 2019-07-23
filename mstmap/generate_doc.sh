doxygen Doxyfile
# npm install moxygen -g
mkdir doc/md
moxygen --html-anchors --output doc/md/tmap.md doc/xml

cd python-doc
make html
make markdown
make latexpdf